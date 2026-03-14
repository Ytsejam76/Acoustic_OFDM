% Copyright (c) 2026 Elias S. G. Carotti

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
            pos = pos + p.wake_miss_hop;
            continue;
        end

        segment = rx_audio(start_idx:end);
        [ok, pkt_info, consumed, dbg_one] = rx_one_packet(segment, p);
        dbg_one.wake_score = score;
        dbg_one.wake_start = start_idx;
        last_dbg = dbg_one;

        if ok
            packets{end+1} = pkt_info; %#ok<AGROW>
            pos = max(start_idx + consumed, start_idx + p.wake_rearm_hop);
        else
            pos = start_idx + p.wake_retry_hop;
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
    [data_bins, ~, ~] = ofdm_bin_plan(p);
    num_data_carriers = numel(data_bins);
    if num_data_carriers <= 0
        return;
    end
    bps = bits_per_modulation(p.modulation);
    max_payload_bytes = p.packet_payload_bytes + 16;
    max_bits = max_payload_bytes * 8;
    max_data_ofdm = ceil(max_bits / (num_data_carriers * bps)) + 4;
    needed_after_sync = 2*p.sync_half_len + sym_len + max_data_ofdm*sym_len;
    search_end = min(length(rx_segment), coarse_start + p.sync_search_len + needed_after_sync);
    pre = min(coarse_start - 1, p.rx_preroll);
    chunk = rx_segment(coarse_start-pre:search_end);
    rbb_full = iq_downconvert_and_lpf(chunk, p);
    rbb = rbb_full(pre+1:end);
    dbg.rbb = rbb;

    [metric, peak_idx_raw] = sync_metric_known_preamble(rbb, p);
    [sync_cands, sync_scores, peak_ratio] = select_sync_candidates(metric, p);
    if isempty(sync_cands)
        return;
    end
    peak_idx = sync_cands(1);
    cfo_est_hz = estimate_cfo_from_repeat(rbb, peak_idx, p.sync_half_len, p.fs);
    dbg.sync_metric = metric;
    dbg.peak_idx = peak_idx;
    dbg.peak_idx_raw = peak_idx_raw;
    dbg.cfo_est_hz = cfo_est_hz;
    dbg.sync_candidates = sync_cands(:);
    dbg.sync_candidate_scores = sync_scores(:);
    dbg.sync_peak_ratio = peak_ratio;

    if isfield(p, 'oracle_sync_start') && ~isempty(p.oracle_sync_start)
        sync_cands = max(1, round(double(p.oracle_sync_start)));
        dbg.peak_idx = sync_cands(1);
        use_forced_cfo = isfield(p, 'oracle_cfo_est_hz') && ~isempty(p.oracle_cfo_est_hz);
        if use_forced_cfo
            cfo_force = double(p.oracle_cfo_est_hz);
        else
            cfo_force = [];
        end
    else
        cfo_force = [];
    end

    max_cfo_hz = 300;
    if isfield(p, 'max_cfo_hz')
        max_cfo_hz = p.max_cfo_hz;
    end

    best_found = false;
    best_score = -Inf;
    best_pkt = struct();
    best_consumed_bb = 0;
    best_dbg_try = struct();
    best_sync = 0;
    best_ci = 0;
    best_cfo = 0;

    for ci = 1:numel(sync_cands)
        sync_start = sync_cands(ci);
        if sync_start < 1
            continue;
        end

        if isempty(cfo_force)
            cfo_est = estimate_cfo_from_repeat(rbb, sync_start, p.sync_half_len, p.fs);
            cfo_list = p.cfo_grid_hz(:).';
            if abs(cfo_est) <= max_cfo_hz
                cfo_list = [cfo_est, cfo_list];
            end
            cfo_list = unique(cfo_list, 'stable');
        else
            cfo_list = cfo_force;
        end

        for cfi = 1:numel(cfo_list)
            cfo_try = cfo_list(cfi);
            [ok_try, pkt_try, consumed_bb, dbg_try] = try_decode_from_sync(rbb, sync_start, cfo_try, p);
            if ok_try
                score = sync_scores(min(ci, numel(sync_scores))) / (1 + abs(cfo_try) / max(max_cfo_hz, 1));
                if score > best_score
                    best_found = true;
                    best_score = score;
                    best_pkt = pkt_try;
                    best_consumed_bb = consumed_bb;
                    best_dbg_try = dbg_try;
                    best_sync = sync_start;
                    best_ci = ci;
                    best_cfo = cfo_try;
                end
            end
        end
    end

    if best_found
        ok = true;
        pkt_info = best_pkt;
        consumed = coarse_start - 1 + best_consumed_bb;
        dbg.cfo_est_hz = best_cfo;
        dbg.sync_candidate_used = best_ci;
        dbg.peak_idx = best_sync;
        dbg.rbb_cfo = best_dbg_try.rbb_cfo;
        dbg.Hest = best_dbg_try.Hest;
        dbg.train_rx_raw = best_dbg_try.train_rx_raw;
        dbg.train_rx_eq = best_dbg_try.train_rx_eq;
        dbg.rx_syms_raw = best_dbg_try.rx_syms_raw;
        dbg.rx_syms_eq = best_dbg_try.rx_syms_eq;
        dbg.pll_phase_hist = best_dbg_try.pll_phase_hist;
        dbg.pll_freq_hist = best_dbg_try.pll_freq_hist;
        dbg.pll_err_hist = best_dbg_try.pll_err_hist;
    end
end

function [found, best_idx, best_score] = find_wake_tone(x, pos, p)
    found = false; best_idx = -1; best_score = 0;
    [ref, ~] = make_wake_tone_ref(p);
    N = length(ref);
    ref_energy = sum(abs(ref).^2) + 1e-12;
    hop = max(1, floor(N / 8));
    i_end = min(length(x)-N+1, pos + p.wake_search_len);
    min_score = 0.005;
    if isfield(p, 'wake_min_score')
        min_score = p.wake_min_score;
    end
    thr = max(min_score, p.detect_threshold);

    first_hit = -1;
    for i = pos:hop:i_end
        seg = x(i:i+N-1);
        corrv = sum(conj(ref) .* seg);
        norm_score = abs(corrv)^2 / ((sum(abs(seg).^2) + 1e-12) * ref_energy);
        if norm_score > best_score
            best_score = norm_score;
            best_idx = i;
        end
        if first_hit < 0 && norm_score >= thr
            first_hit = i;
            break;
        end
    end

    if first_hit > 0
        % Refine around the first crossing to avoid jumping to later false wakes.
        j0 = max(pos, first_hit - 2*hop);
        j1 = min(i_end, first_hit + 2*hop);
        best_score = -Inf;
        best_idx = first_hit;
        for j = j0:j1
            seg = x(j:j+N-1);
            corrv = sum(conj(ref) .* seg);
            norm_score = abs(corrv)^2 / ((sum(abs(seg).^2) + 1e-12) * ref_energy);
            if norm_score > best_score
                best_score = norm_score;
                best_idx = j;
            end
        end
        found = (best_score >= thr);
    end
end

function [ref, preN] = make_wake_tone_ref(p)
    N = round(p.wake_ms * 1e-3 * p.fs);
    preN = round(p.wake_ref_pre_ms * 1e-3 * p.fs);
    postN = round(p.wake_guard_ms * 1e-3 * p.fs);
    n = (0:N-1).';
    if isfield(p, 'use_chirp_sync') && p.use_chirp_sync
        t = n / p.fs;
        T = max(t(end), 1/p.fs);
        k = (p.sync_chirp_f1 - p.sync_chirp_f0) / T;
        phase = 2*pi * (p.sync_chirp_f0 * t + 0.5 * k * t.^2);
        w = sin(phase);
    else
        w = sin(2*pi * p.wake_freq * n / p.fs);
    end
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

function [cands, scores, peak_ratio] = select_sync_candidates(metric, p)
    if isempty(metric)
        cands = 1;
        scores = 1;
        peak_ratio = Inf;
        return;
    end

    max_search = numel(metric);
    if isfield(p, 'sync_search_len')
        max_search = min(max_search, p.sync_search_len);
    end
    mm = metric(1:max_search);
    max_m = max(mm);
    thr_rel = p.sync_peak_rel_threshold * max_m;
    thr_abs = p.sync_metric_abs_threshold;
    thr = max(thr_rel, thr_abs);

    num_cands = 6;
    if isfield(p, 'sync_num_candidates')
        num_cands = p.sync_num_candidates;
    end
    min_sep = max(1, round(p.sync_half_len / 4));
    if isfield(p, 'sync_candidate_min_sep')
        min_sep = max(1, round(p.sync_candidate_min_sep));
    end

    [vals, idx] = sort(mm, 'descend');
    if numel(vals) >= 2
        peak_ratio = vals(1) / (vals(2) + 1e-12);
    else
        peak_ratio = Inf;
    end

    expected_idx = 1;
    if isfield(p, 'sync_expected_idx')
        expected_idx = p.sync_expected_idx;
    end
    expected_span = 220;
    if isfield(p, 'sync_expected_span')
        expected_span = max(1, p.sync_expected_span);
    end
    expected_weight = 1.0;
    if isfield(p, 'sync_expected_weight')
        expected_weight = max(0, p.sync_expected_weight);
    end
    dist = abs(idx - expected_idx);
    rank_val = vals .* exp(-expected_weight * (dist / expected_span));

    [~, order] = sort(rank_val, 'descend');
    cands = [];
    scores = [];

    % Always keep one strong early candidate near the expected sync position.
    exp_lo = max(1, round(expected_idx - expected_span));
    exp_hi = min(numel(mm), round(expected_idx + expected_span));
    mm_exp = mm(exp_lo:exp_hi);
    [exp_best_val, exp_off] = max(mm_exp);
    exp_best_idx = exp_lo + exp_off - 1;
    if exp_best_val >= thr
        cands(end+1) = exp_best_idx; %#ok<AGROW>
        scores(end+1) = exp_best_val; %#ok<AGROW>
    end

    for oi = 1:numel(order)
        i = order(oi);
        if vals(i) < thr
            break;
        end
        if isempty(cands) || min(abs(idx(i) - cands)) >= min_sep
            cands(end+1) = idx(i); %#ok<AGROW>
            scores(end+1) = rank_val(i); %#ok<AGROW>
            if numel(cands) >= num_cands
                break;
            end
        end
    end
    if isempty(cands)
        [~, best] = max(mm);
        cands = best;
        scores = mm(best);
    end
end

function [ok, pkt_info, consumed_bb, dbg] = try_decode_from_sync(rbb, sync_start, cfo_est_hz, p)
    ok = false;
    pkt_info = struct();
    consumed_bb = 0;
    dbg = struct();

    n = (0:length(rbb)-1).';
    rbb_cfo = rbb .* exp(-1j * 2*pi * cfo_est_hz * n / p.fs);
    dbg.rbb_cfo = rbb_cfo;

    xsync_len = 2 * p.sync_half_len;
    train_len = p.Nfft + p.Ncp;
    train_start = sync_start + xsync_len;
    train_end = train_start + train_len - 1;
    if train_end > length(rbb_cfo)
        return;
    end

    rtrain = rbb_cfo(train_start:train_end);
    rtrain = rtrain(p.Ncp+1:end);
    Ytrain = fft(rtrain, p.Nfft);
    [data_bins, pilot_bins, used_bins] = ofdm_bin_plan(p);
    num_data_carriers = numel(data_bins);
    if num_data_carriers <= 0
        return;
    end
    [has_pilot, pilot_pos, data_pos] = bin_positions(used_bins, pilot_bins, data_bins);
    train_known = known_training_symbols(numel(used_bins), p.modulation);
    Hest = Ytrain(used_bins) ./ train_known;
    dbg.Hest = Hest;
    dbg.train_rx_raw = Ytrain(used_bins);
    dbg.train_rx_eq = Ytrain(used_bins) ./ Hest;

    data_start = train_end + 1;
    max_payload_bytes = p.packet_payload_bytes + 16;
    max_bits = max_payload_bytes * 8;
    bps = bits_per_modulation(p.modulation);
    max_data_ofdm = ceil(max_bits / (num_data_carriers * bps)) + 2;

    sym_len = p.Nfft + p.Ncp;
    total_needed = data_start + max_data_ofdm * sym_len - 1;
    if total_needed > length(rbb_cfo)
        max_data_ofdm = floor((length(rbb_cfo) - data_start + 1) / sym_len);
    end
    if max_data_ofdm <= 0
        return;
    end

    rx_syms = [];
    dbg.rx_syms_raw = [];
    dbg.rx_syms_eq = [];
    pll_on = isfield(p, 'pll_enable') && p.pll_enable;
    pll_kp = 0.05;
    pll_ki = 0.002;
    if isfield(p, 'pll_kp')
        pll_kp = p.pll_kp;
    end
    if isfield(p, 'pll_ki')
        pll_ki = p.pll_ki;
    end
    phase_hat = 0;
    freq_hat = 0;
    dbg.pll_phase_hist = [];
    dbg.pll_freq_hist = [];
    dbg.pll_err_hist = [];
    dbg.pilot_gain_hist = [];
    for i = 1:max_data_ofdm
        s0 = data_start + (i-1)*sym_len;
        s1 = s0 + sym_len - 1;
        if s1 > length(rbb_cfo)
            break;
        end
        rt = rbb_cfo(s0:s1);
        rt = rt(p.Ncp+1:end);
        Y = fft(rt, p.Nfft);
        Xraw_used = Y(used_bins);
        Xeq_used = Xraw_used ./ Hest;
        if has_pilot
            pref = known_pilot_symbols(numel(pilot_bins), i);
            prx = Xeq_used(pilot_pos);
            g = sum(prx .* conj(pref)) / (sum(abs(pref).^2) + 1e-12);
            if isfinite(g) && abs(g) > 1e-12
                Xeq_used = Xeq_used / g;
            else
                g = 1;
            end
            dbg.pilot_gain_hist(end+1, 1) = g; %#ok<AGROW>
        end
        Xraw = Xraw_used(data_pos);
        Xeq = Xeq_used(data_pos);
        if pll_on
            Xeq = Xeq .* exp(-1j * phase_hat);
            Xdec = hard_slice_symbols(Xeq, p.modulation);
            ph_err = angle(sum(Xeq .* conj(Xdec)));
            if ~isfinite(ph_err)
                ph_err = 0;
            end
            freq_hat = freq_hat + pll_ki * ph_err;
            phase_hat = phase_hat + freq_hat + pll_kp * ph_err;
            dbg.pll_phase_hist(end+1, 1) = phase_hat; %#ok<AGROW>
            dbg.pll_freq_hist(end+1, 1) = freq_hat; %#ok<AGROW>
            dbg.pll_err_hist(end+1, 1) = ph_err; %#ok<AGROW>
        end
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
    used_ofdm = ceil(num_data_bits / (num_data_carriers * bits_per_modulation(p.modulation)));
    num_used_syms = used_ofdm * num_data_carriers;
    dbg.rx_syms_raw = dbg.rx_syms_raw(1:min(num_used_syms, numel(dbg.rx_syms_raw)));
    dbg.rx_syms_eq = dbg.rx_syms_eq(1:min(num_used_syms, numel(dbg.rx_syms_eq)));
    consumed_bb = (sync_start - 1) + xsync_len + train_len + used_ofdm * sym_len;
    ok = true;
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

function out = hard_slice_symbols(syms, modulation)
    switch upper(modulation)
        case 'BPSK'
            out = ones(size(syms));
            out(real(syms) < 0) = -1;
        case 'QPSK'
            re = ones(size(syms));
            im = ones(size(syms));
            re(real(syms) < 0) = -1;
            im(imag(syms) < 0) = -1;
            out = (re + 1j * im) / sqrt(2);
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

function p = known_pilot_symbols(N, sym_idx)
    idx = (0:N-1).';
    p = exp(1j * pi/2 * mod((sym_idx - 1) + idx, 4));
end

function [data_bins, pilot_bins, used_bins] = ofdm_bin_plan(p)
    used_bins = p.used_bins(:).';
    pilot_bins = select_pilot_bins(p, used_bins);
    data_bins = setdiff(used_bins, pilot_bins, 'stable');
end

function pilot_bins = select_pilot_bins(p, used_bins)
    pilot_bins = [];
    if ~pilots_enabled(p)
        return;
    end

    if isfield(p, 'pilot_bins') && ~isempty(p.pilot_bins)
        candidates = p.pilot_bins(:).';
    else
        candidates = used_bins;
    end
    pilot_bins = intersect(used_bins, candidates, 'stable');

    np = requested_num_pilots(p);
    if ~isempty(np)
        np = max(0, min(np, numel(pilot_bins)));
        pilot_bins = pilot_bins(1:np);
    end
end

function np = requested_num_pilots(p)
    np = [];
    if isfield(p, 'num_pilots') && ~isempty(p.num_pilots)
        np = round(double(p.num_pilots(1)));
        if ~isfinite(np)
            np = [];
        end
    end
end

function on = pilots_enabled(p)
    if isfield(p, 'use_pilots') && ~isempty(p.use_pilots)
        on = logical(p.use_pilots);
        return;
    end
    on = strcmpi(p.modulation, 'QPSK');
end

function [has_pilot, pilot_pos, data_pos] = bin_positions(used_bins, pilot_bins, data_bins)
    if isempty(pilot_bins)
        has_pilot = false;
        pilot_pos = [];
    else
        [tfp, pilot_pos] = ismember(pilot_bins, used_bins);
        has_pilot = all(tfp);
        if ~has_pilot
            pilot_pos = [];
        end
    end
    [tfd, data_pos] = ismember(data_bins, used_bins);
    if ~all(tfd)
        data_pos = [];
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
    p.pilot_bins = [2 4];
    p.num_pilots = [];
    p.modulation = 'BPSK';
    p.use_pilots = [];
    p.wake_ms = 12;
    p.wake_freq = 16500;
    p.wake_guard_ms = 4;
    p.use_chirp_sync = true;
    p.sync_chirp_f0 = 4000;
    p.sync_chirp_f1 = 8000;
    p.sync_half_len = 64;
    p.packet_payload_bytes = 24;
    p.detect_threshold = 0.005;
    p.wake_min_score = 0.005;
    p.wake_search_len = 8000;
    p.wake_miss_hop = round(0.012 * p.fs);
    p.wake_retry_hop = round(0.002 * p.fs);
    p.wake_rearm_hop = round(0.008 * p.fs);
    p.wake_ref_pre_ms = 15;
    p.sync_search_len = 1500;
    p.disable_lpf = false;
    p.disable_modulation = false;
    p.oracle_wake_start = [];
    p.oracle_sync_start = [];
    p.oracle_cfo_est_hz = [];
    p.sync_peak_rel_threshold = 0.80;
    p.sync_metric_abs_threshold = 0.02;
    p.sync_num_candidates = 20;
    p.sync_candidate_min_sep = 16;
    p.max_cfo_hz = 300;
    p.cfo_grid_hz = -12:2:12;
    p.sync_peak_ratio_min = 1.20;
    p.sync_expected_idx = 1;
    p.sync_expected_span = 220;
    p.sync_expected_weight = 1.0;
    p.rx_preroll = 128;
    p.pll_enable = true;
    p.pll_kp = 0.05;
    p.pll_ki = 0.002;
end
