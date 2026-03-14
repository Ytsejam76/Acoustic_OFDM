% Copyright (c) 2026 Elias S. G. Carotti

function [tx_audio, meta] = ofdm_encoder(payload, p)
% OFDM_ENCODER  Build a short packetized acoustic OFDM burst.
%
%   [tx_audio, meta] = ofdm_encoder(payload, p)
%
% payload : uint8 vector
% p       : parameter struct; if omitted, defaults are used
%
% tx_audio: real passband waveform suitable for playback or simulation
% meta    : debug metadata

    if nargin < 2 || isempty(p)
        p = default_params();
    end

    payload = uint8(payload(:));
    if ~isfield(p, 'modulation')
        p.modulation = 'BPSK';
    end
    if ~isfield(p, 'disable_modulation')
        p.disable_modulation = false;
    end

    chunks = split_payload(payload, p.packet_payload_bytes);
    num_packets = numel(chunks);

    tx_audio = [];
    packet_meta = cell(num_packets, 1);

    for i = 1:num_packets
        pkt_bytes = build_packet_bytes(chunks{i}, uint8(i-1), uint8(num_packets), p);
        [pkt_audio, pkt_dbg] = tx_one_packet(pkt_bytes, p);

        pre_sil  = zeros(round(0.015 * p.fs), 1);
        post_sil = zeros(round(0.020 * p.fs), 1);

        tx_audio = [tx_audio; pre_sil; pkt_audio; post_sil]; %#ok<AGROW>
        packet_meta{i} = pkt_dbg;
    end

    m = max(abs(tx_audio));
    if m > 0
        tx_audio = 0.85 * tx_audio / m;
    end

    meta = struct();
    meta.num_packets = num_packets;
    meta.packet_meta = packet_meta;
    meta.params = p;
end

function [tx_audio, dbg] = tx_one_packet(pkt_bytes, p)
    bits = bytes_to_bits(pkt_bytes);
    [data_bins, pilot_bins, used_bins] = ofdm_bin_plan(p);

    [payload_syms, bps] = map_bits(bits, p.modulation); %#ok<NASGU>
    num_data_carriers = numel(data_bins);
    if num_data_carriers <= 0
        error('No data carriers left after pilot allocation');
    end
    syms_per_ofdm = num_data_carriers;
    num_data_ofdm = ceil(numel(payload_syms) / syms_per_ofdm);

    pad_len = num_data_ofdm * syms_per_ofdm - numel(payload_syms);
    if pad_len > 0
        payload_syms = [payload_syms; zeros(pad_len, 1)]; %#ok<AGROW>
    end
    payload_grid = reshape(payload_syms, syms_per_ofdm, num_data_ofdm);

    Xtrain = zeros(p.Nfft, 1);
    train_known = known_training_symbols(numel(used_bins), p.modulation);
    Xtrain(used_bins) = train_known;
    xtrain = ifft(Xtrain, p.Nfft);
    xtrain_cp = [xtrain(end-p.Ncp+1:end); xtrain];

    xdata = [];
    for i = 1:num_data_ofdm
        X = zeros(p.Nfft, 1);
        X(data_bins) = payload_grid(:, i);
        if ~isempty(pilot_bins)
            X(pilot_bins) = known_pilot_symbols(numel(pilot_bins), i);
        end
        xt = ifft(X, p.Nfft);
        xt_cp = [xt(end-p.Ncp+1:end); xt];
        xdata = [xdata; xt_cp]; %#ok<AGROW>
    end

    syncA = known_sync_half(p.sync_half_len);
    xsync = [syncA; syncA];

    xbb = [xsync; xtrain_cp; xdata];
    if p.disable_modulation
        passband = xbb;
    else
        passband = iq_upconvert(xbb, p.fs, p.fc);
    end
    wake = make_wake_tone(p);
    guard = zeros(round(p.wake_guard_ms * 1e-3 * p.fs), 1);

    tx_audio = [wake; guard; passband];
    m = max(abs(tx_audio));
    if m > 0
        tx_audio = 0.85 * tx_audio / m;
    end

    dbg = struct();
    dbg.packet_bytes = pkt_bytes;
    dbg.xbb = xbb;
    dbg.train_known = train_known;
    dbg.num_data_ofdm = num_data_ofdm;
    dbg.payload_grid = payload_grid;
    dbg.data_bins = data_bins;
    dbg.pilot_bins = pilot_bins;
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
    p.session_id = uint16(1234);
    p.disable_modulation = false;
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

function p = known_pilot_symbols(N, sym_idx)
    idx = (0:N-1).';
    p = exp(1j * pi/2 * mod((sym_idx - 1) + idx, 4));
end

function chunks = split_payload(payload, chunk_size)
    n = numel(payload);
    m = ceil(n / chunk_size);
    chunks = cell(m,1);
    for i = 1:m
        s0 = (i-1)*chunk_size + 1;
        s1 = min(i*chunk_size, n);
        chunks{i} = payload(s0:s1);
    end
end

function pkt = build_packet_bytes(payload, frag_index, frag_count, p)
    magic = uint8([hex2dec('A5'); hex2dec('5A')]);
    version = uint8(1);
    switch upper(p.modulation)
        case 'BPSK'
            mod_id = uint8(1);
        case 'QPSK'
            mod_id = uint8(2);
        otherwise
            error('Unsupported modulation');
    end
    sid = uint16_to_bytes(p.session_id);
    plen = uint8(numel(payload));

    header = [magic; version; mod_id; sid(:); frag_index; frag_count; plen];
    body = [header; payload(:)];
    crc = crc16_ccitt(body);
    pkt = [body; uint16_to_bytes(crc)];
end

function [syms, bps] = map_bits(bits, modulation)
    switch upper(modulation)
        case 'BPSK'
            bps = 1;
            bits = bits(:);
            syms = 1 - 2*double(bits);
        case 'QPSK'
            bps = 2;
            bits = bits(:);
            if mod(length(bits), 2) ~= 0
                bits = [bits; 0];
            end
            b = reshape(bits, 2, []).';
            syms = (1 - 2*double(b(:,1))) + 1j*(1 - 2*double(b(:,2)));
            syms = syms / sqrt(2);
        otherwise
            error('Unsupported modulation');
    end
end

function s = known_training_symbols(N, modulation)
    switch upper(modulation)
        case 'BPSK'
            s = ones(N, 1);
            s(2:2:end) = -1;
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

function s = known_sync_half(L)
    idx = (0:L-1).';
    a = exp(1j * pi/2 * mod(3*idx + 1, 4));
    s = a / sqrt(mean(abs(a).^2));
end

function y = iq_upconvert(xbb, fs, fc)
    n = (0:length(xbb)-1).';
    y = real(xbb .* exp(1j * 2*pi * fc * n / fs));
end

function w = make_wake_tone(p)
    N = round(p.wake_ms * 1e-3 * p.fs);
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

function bits = bytes_to_bits(bytes)
    bytes = uint8(bytes(:));
    bits = zeros(numel(bytes)*8, 1, 'uint8');
    k = 1;
    for i = 1:numel(bytes)
        for b = 7:-1:0
            bits(k) = bitget(bytes(i), b+1);
            k = k + 1;
        end
    end
end

function b = uint16_to_bytes(x)
    x = uint16(x);
    b = uint8([bitshift(x, -8); bitand(x, 255)]);
end
