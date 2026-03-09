function [result, dbg] = ofdm_decoder(rx_audio, p)
% OFDM_DECODER  Decode one or more short acoustic OFDM packets.
%
%   [result, dbg] = ofdm_decoder(rx_audio, p)

    if nargin < 2 || isempty(p)
        p = default_params();
    end

    rx_audio = rx_audio(:);
    packets = {};
    pos = 1;
    last_dbg = struct();

    while pos < length(rx_audio)
        if isfield(p, 'oracle_wake_start') && ~isempty(p.oracle_wake_start) && pos == 1
            found = true;
            start_idx = max(1, round(double(p.oracle_wake_start)));
            score = NaN;
        else
            [found, start_idx, score] = find_wake_tone(rx_audio, pos, p);
        end
        if ~found
            break;
        end

        segment = rx_audio(start_idx:end);
        [ok, pkt_info, consumed, dbg_one] = rx_one_packet(segment, p);
        dbg_one.wake_score = score;
        dbg_one.wake_start = start_idx;
        last_dbg = dbg_one;

        if ok
            packets{end+1} = pkt_info; %#ok<AGROW>
            pos = start_idx + consumed;
        else
            pos = start_idx + round(0.03 * p.fs);
        end
    end

    result = reconstruct_message_from_packets(packets);
    dbg = last_dbg;
end

function [ok, pkt_info, consumed, dbg] = rx_one_packet(rx_segment, p)
    ok = false;
    pkt_info = struct();
    dbg = struct();
    consumed = 0;

    wake_len = round(p.wake_ms * 1e-3 * p.fs);
    guard_len = round(p.wake_guard_ms * 1e-3 * p.fs);
    coarse_start = wake_len + guard_len + 1;
    if length(rx_segment) < coarse_start + 200
        return;
    end

    sym_len = p.Nfft + p.Ncp;
    bps = bits_per_modulation(p.modulation);
    max_payload_bytes = p.packet_payload_bytes + 16;
    max_bits = max_payload_bytes * 8;
    max_data_ofdm = ceil(max_bits / (numel(p.used_bins) * bps)) + 4;
    needed_after_sync = 2*p.sync_half_len + sym_len + max_data_ofdm*sym_len;
    search_end = min(length(rx_segment), coarse_start + p.sync_search_len + needed_after_sync);
    chunk = rx_segment(coarse_start:search_end);
    rbb = iq_downconvert_and_lpf(chunk, p);
    dbg.rbb = rbb;

    [metric, peak_idx_raw] = sync_metric_known_preamble(rbb, p);
    peak_idx = select_sync_peak(metric, p);
    cfo_est_hz = estimate_cfo_from_repeat(rbb, peak_idx, p.sync_half_len, p.fs);
    dbg.sync_metric = metric;
    dbg.peak_idx = peak_idx;
    dbg.peak_idx_raw = peak_idx_raw;
    dbg.cfo_est_hz = cfo_est_hz;

    if isfield(p, 'oracle_sync_start') && ~isempty(p.oracle_sync_start)
        sync_start = max(1, round(double(p.oracle_sync_start)));
        dbg.peak_idx = sync_start;
    else
        sync_start = peak_idx;
    end
    if isfield(p, 'oracle_cfo_est_hz') && ~isempty(p.oracle_cfo_est_hz)
        cfo_est_hz = double(p.oracle_cfo_est_hz);
        dbg.cfo_est_hz = cfo_est_hz;
    end
    if sync_start < 1
        return;
    end

    n = (0:length(rbb)-1).';
    rbb = rbb .* exp(-1j * 2*pi * cfo_est_hz * n / p.fs);
    dbg.rbb_cfo = rbb;

    xsync_len = 2 * p.sync_half_len;
    train_len = p.Nfft + p.Ncp;
    train_start = sync_start + xsync_len;
    train_end = train_start + train_len - 1;
    if train_end > length(rbb)
        return;
    end

    rtrain = rbb(train_start:train_end);
    rtrain = rtrain(p.Ncp+1:end);
    Ytrain = fft(rtrain, p.Nfft);
    train_known = known_training_symbols(numel(p.used_bins), p.modulation);
    Hest = Ytrain(p.used_bins) ./ train_known;

    dbg.Hest = Hest;
    dbg.train_rx_raw = Ytrain(p.used_bins);
    dbg.train_rx_eq = Ytrain(p.used_bins) ./ Hest;

    data_start = train_end + 1;
    max_payload_bytes = p.packet_payload_bytes + 16;
    max_bits = max_payload_bytes * 8;
    bps = bits_per_modulation(p.modulation);
    max_data_ofdm = ceil(max_bits / (numel(p.used_bins) * bps)) + 2;

    sym_len = p.Nfft + p.Ncp;
    total_needed = data_start + max_data_ofdm * sym_len - 1;
    if total_needed > length(rbb)
        max_data_ofdm = floor((length(rbb) - data_start + 1) / sym_len);
    end
    if max_data_ofdm <= 0
        return;
    end

    rx_syms = [];
    dbg.rx_syms_raw = [];
    dbg.rx_syms_eq = [];

    for i = 1:max_data_ofdm
        s0 = data_start + (i-1)*sym_len;
        s1 = s0 + sym_len - 1;
        if s1 > length(rbb)
            break;
        end
        rt = rbb(s0:s1);
        rt = rt(p.Ncp+1:end);
        Y = fft(rt, p.Nfft);
        Xraw = Y(p.used_bins);
        Xeq = Xraw ./ Hest;
        rx_syms = [rx_syms; Xeq(:)]; %#ok<AGROW>
        dbg.rx_syms_raw = [dbg.rx_syms_raw; Xraw(:)]; %#ok<AGROW>
        dbg.rx_syms_eq = [dbg.rx_syms_eq; Xeq(:)]; %#ok<AGROW>
    end

    rx_bits = demap_bits(rx_syms, p.modulation);
    rx_bytes = bits_to_bytes(rx_bits);
    [valid, pkt_info, pkt_total_bytes] = parse_packet_bytes(rx_bytes);
    if ~valid
        return;
    end

    num_data_bits = pkt_total_bytes * 8;
    used_ofdm = ceil(num_data_bits / (numel(p.used_bins) * bits_per_modulation(p.modulation)));
    num_used_syms = used_ofdm * numel(p.used_bins);
    dbg.rx_syms_raw = dbg.rx_syms_raw(1:min(num_used_syms, numel(dbg.rx_syms_raw)));
    dbg.rx_syms_eq = dbg.rx_syms_eq(1:min(num_used_syms, numel(dbg.rx_syms_eq)));

    consumed_bb = (sync_start - 1) + xsync_len + train_len + used_ofdm * sym_len;
    consumed = coarse_start - 1 + consumed_bb;
    ok = true;
end

function [found, best_idx, best_score] = find_wake_tone(x, pos, p)
    found = false; best_idx = -1; best_score = 0;
    [ref, preN] = make_wake_tone_ref(p);
    N = length(ref);
    ref_energy = sum(abs(ref).^2) + 1e-12;
    hop = max(1, floor(N / 8));
    i_end = min(length(x)-N+1, pos + p.wake_search_len);
    for i = pos:hop:i_end
        seg = x(i:i+N-1);
        corrv = sum(conj(ref) .* seg);
        norm_score = abs(corrv)^2 / ((sum(abs(seg).^2) + 1e-12) * ref_energy);
        if norm_score > best_score
            best_score = norm_score;
            best_idx = i;
        end
    end
    if best_score > p.detect_threshold
        found = true;
        best_idx = best_idx + preN;
    end
end

function [ref, preN] = make_wake_tone_ref(p)
    N = round(p.wake_ms * 1e-3 * p.fs);
    preN = round(p.wake_ref_pre_ms * 1e-3 * p.fs);
    postN = round(p.wake_guard_ms * 1e-3 * p.fs);
    n = (0:N-1).';
    w = sin(2*pi * p.wake_freq * n / p.fs);
    ramp = min(round(0.001 * p.fs), floor(N/4));
    env = ones(N,1);
    if ramp > 1
        r = (0:ramp-1).' / ramp;
        env(1:ramp) = r;
        env(end-ramp+1:end) = flipud(r);
    end
    w = 0.7 * w .* env;
    ref = [zeros(preN, 1); w; zeros(postN, 1)];
end

function [metric, peak_idx_raw] = sync_metric_known_preamble(r, p)
    L = p.sync_half_len;
    ref_half = known_sync_half(L);
    ref = [ref_half; ref_half];
    ref_energy = sum(abs(ref).^2) + 1e-12;

    N = length(r);
    M = max(1, N - 2*L + 1);
    metric = zeros(M, 1);
    for d = 1:M
        seg = r(d:d+2*L-1);
        c = sum(seg .* conj(ref));
        e = sum(abs(seg).^2) + 1e-12;
        metric(d) = abs(c)^2 / (e * ref_energy);
    end
    [~, peak_idx_raw] = max(metric);
end

function cfo_est_hz = estimate_cfo_from_repeat(r, sync_start, L, fs)
    if sync_start < 1 || (sync_start + 2*L - 1) > length(r)
        cfo_est_hz = 0;
        return;
    end
    a = r(sync_start:sync_start+L-1);
    b = r(sync_start+L:sync_start+2*L-1);
    P = sum(a .* conj(b));
    cfo_est_hz = -angle(P) * fs / (2*pi*L);
end

function peak_idx = select_sync_peak(metric, p)
    if isempty(metric)
        peak_idx = 1;
        return;
    end

    max_m = max(metric);
    thr_rel = p.sync_peak_rel_threshold * max_m;
    thr_abs = p.sync_metric_abs_threshold;
    thr = max(thr_rel, thr_abs);
    cand = find(metric >= thr, 1, 'first');
    if isempty(cand)
        [~, peak_idx] = max(metric);
    else
        peak_idx = cand;
    end
end

function s = known_sync_half(L)
    idx = (0:L-1).';
    a = exp(1j * pi/2 * mod(3*idx + 1, 4));
    s = a / sqrt(mean(abs(a).^2));
end

function rbb = iq_downconvert_and_lpf(y, p)
    fs = p.fs;
    fc = p.fc;
    n = (0:length(y)-1).';
    if isfield(p, 'disable_modulation') && p.disable_modulation
        rmix = y;
    else
        rmix = y .* exp(-1j * 2*pi * fc * n / fs);
    end

    if isfield(p, 'disable_lpf') && p.disable_lpf
        rbb = rmix;
        return;
    end

    L = 65;
    % Keep all configured OFDM bins plus margin after downconversion.
    df = fs / p.Nfft;
    max_sc_hz = max(p.used_bins) * df;
    bw = min(0.45 * fs, max(5000, 1.35 * max_sc_hz));
    h = fir_lowpass(L, bw, fs);
    rbb = filter(h, 1, rmix);
    gd = floor((L-1)/2);
    rbb = [rbb(gd+1:end); zeros(gd,1)];
end

function h = fir_lowpass(L, fc, fs)
    n = -(L-1)/2:(L-1)/2;
    h = 2*fc/fs * sinc(2*fc*n/fs);
    w = 0.54 - 0.46*cos(2*pi*(0:L-1)/(L-1));
    h = h(:) .* w(:);
    h = h / sum(h);
end

function result = reconstruct_message_from_packets(packets)
    result = struct();
    result.success = false;
    result.payload = uint8([]);
    result.num_packets_ok = numel(packets);
    result.num_packets_total = 0;
    if isempty(packets)
        return;
    end
    frag_count = double(packets{1}.frag_count);
    result.num_packets_total = frag_count;
    by_index = cell(frag_count, 1);
    for i = 1:numel(packets)
        idx = double(packets{i}.frag_index) + 1;
        if idx >= 1 && idx <= frag_count
            by_index{idx} = packets{i}.payload;
        end
    end
    for i = 1:frag_count
        if isempty(by_index{i})
            return;
        end
    end
    payload = uint8([]);
    for i = 1:frag_count
        payload = [payload; by_index{i}(:)]; %#ok<AGROW>
    end
    result.payload = payload;
    result.success = true;
end

function [valid, info, total_bytes] = parse_packet_bytes(rx_bytes)
    valid = false;
    info = struct();
    total_bytes = 0;
    if numel(rx_bytes) < 11
        return;
    end
    if ~(rx_bytes(1) == hex2dec('A5') && rx_bytes(2) == hex2dec('5A'))
        return;
    end
    version = rx_bytes(3);
    sid = bytes_to_uint16(rx_bytes(5:6));
    frag_index = rx_bytes(7);
    frag_count = rx_bytes(8);
    plen = rx_bytes(9);
    total_bytes = 9 + double(plen) + 2;
    if numel(rx_bytes) < total_bytes
        return;
    end
    body = rx_bytes(1:9+double(plen));
    rx_crc = bytes_to_uint16(rx_bytes(10+double(plen):11+double(plen)));
    calc_crc = crc16_ccitt(body);
    if rx_crc ~= calc_crc
        return;
    end
    info.version = version;
    info.session_id = sid;
    info.frag_index = frag_index;
    info.frag_count = frag_count;
    info.payload = rx_bytes(10:9+double(plen));
    valid = true;
end

function bits = demap_bits(syms, modulation)
    switch upper(modulation)
        case 'BPSK'
            bits = uint8(real(syms) < 0);
        case 'QPSK'
            bits = zeros(2*length(syms), 1, 'uint8');
            bits(1:2:end) = uint8(real(syms) < 0);
            bits(2:2:end) = uint8(imag(syms) < 0);
        otherwise
            error('Unsupported modulation');
    end
end

function bps = bits_per_modulation(modulation)
    switch upper(modulation)
        case 'BPSK'
            bps = 1;
        case 'QPSK'
            bps = 2;
        otherwise
            error('Unsupported modulation');
    end
end

function s = known_training_symbols(N, modulation)
    switch upper(modulation)
        case 'BPSK'
            s = ones(N, 1); s(2:2:end) = -1;
        case 'QPSK'
            base = [1+1j; 1-1j; -1+1j; -1-1j] / sqrt(2);
            s = zeros(N,1);
            for i = 1:N
                s(i) = base(mod(i-1,4)+1);
            end
        otherwise
            error('Unsupported modulation');
    end
end

function crc = crc16_ccitt(data)
    crc = uint16(hex2dec('FFFF'));
    poly = uint16(hex2dec('1021'));
    for i = 1:numel(data)
        crc = bitxor(crc, bitshift(uint16(data(i)), 8));
        for b = 1:8
            if bitand(crc, uint16(hex2dec('8000')))
                crc = bitxor(bitshift(crc, 1), poly);
            else
                crc = bitshift(crc, 1);
            end
            crc = bitand(crc, uint16(hex2dec('FFFF')));
        end
    end
end

function bytes = bits_to_bytes(bits)
    bits = uint8(bits(:));
    nbytes = floor(numel(bits) / 8);
    bytes = zeros(nbytes, 1, 'uint8');
    for i = 1:nbytes
        v = uint8(0);
        for b = 1:8
            v = bitshift(v, 1);
            v = bitor(v, bits((i-1)*8 + b));
        end
        bytes(i) = v;
    end
end

function x = bytes_to_uint16(b)
    b = uint8(b(:));
    x = bitor(bitshift(uint16(b(1)), 8), uint16(b(2)));
end

function p = default_params()
    p.fs = 48000;
    p.fc = 17000;
    p.Nfft = 96;
    p.Ncp = 72;
    p.used_bins = [2 3 4 5];
    p.modulation = 'BPSK';
    p.wake_ms = 12;
    p.wake_freq = 16500;
    p.wake_guard_ms = 4;
    p.sync_half_len = 64;
    p.packet_payload_bytes = 24;
    p.detect_threshold = 0.12;
    p.wake_search_len = 6000;
    p.wake_ref_pre_ms = 15;
    p.sync_search_len = 500;
    p.disable_lpf = false;
    p.disable_modulation = false;
    p.oracle_wake_start = [];
    p.oracle_sync_start = [];
    p.oracle_cfo_est_hz = [];
    p.sync_peak_rel_threshold = 0.80;
    p.sync_metric_abs_threshold = 0.02;
end
