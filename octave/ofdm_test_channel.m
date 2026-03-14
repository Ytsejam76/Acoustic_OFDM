% Copyright (c) 2026 Elias S. G. Carotti

function [result, dbg, tx_audio, rx_audio, payload, p] = ofdm_test_channel(varargin)
% OFDM_TEST_CHANNEL  Simulated loop for the acoustic OFDM toy modem.
%
%   [result, dbg, tx_audio, rx_audio, payload, p] = ofdm_test_channel()
%
% This function:
%   1. generates a random payload
%   2. encodes it to passband audio
%   3. applies a simple test channel
%   4. decodes the received audio
%   5. plots waveform, sync metric, and constellation

    p = default_params();
    if nargin >= 1 && ~isempty(varargin{1})
        user_p = varargin{1};
        f = fieldnames(user_p);
        for i = 1:numel(f)
            p.(f{i}) = user_p.(f{i});
        end
    end

    if nargin >= 2 && ~isempty(varargin{2})
        payload = uint8(varargin{2}(:));
    else
        payload = uint8(randi([0 255], p.payload_bytes, 1));
    end

    [tx_audio, enc_meta] = ofdm_encoder(payload, p); %#ok<NASGU>
    rx_audio = simulate_audio_channel(tx_audio, p);
    [result, dbg] = ofdm_decoder(rx_audio, p);
    [ber, bit_errors, bit_total] = compute_payload_ber(payload, result.payload);
    result.ber = ber;
    result.bit_errors = bit_errors;
    result.bit_total_compared = bit_total;

    if p.verbose
        fprintf('\n==== TEST CHANNEL RESULT ====\n');
        fprintf('Modulation: %s\n', p.modulation);
        fprintf('Payload bytes: %d\n', numel(payload));
        fprintf('Decoded packets: %d / %d\n', result.num_packets_ok, result.num_packets_total);
        fprintf('Success: %d\n', result.success);
        if result.success
            fprintf('Payload matches: %d\n', isequal(result.payload(:), payload(:)));
        end
        fprintf('Channel SNR: %.1f dB\n', p.snr_db);
        fprintf('Channel CFO: %.2f Hz\n', p.cfo_hz);
        fprintf('Timing offset: %d samples\n', p.timing_offset);
        if isnan(result.ber)
            fprintf('BER: n/a (no overlapping decoded bytes)\n');
        else
            fprintf('BER: %.6f (%d / %d bit errors)\n', result.ber, result.bit_errors, result.bit_total_compared);
        end
    end

    fig_handles = [];
    if p.make_plots
        fig_handles = make_plots(tx_audio, rx_audio, dbg, p);
    end

    if p.save_images || p.save_decoder_constellation
        ensure_out_dir(p.out_dir);
    end

    if p.save_images
        save_plot_images(fig_handles, p.out_dir);
    end

    if p.save_decoder_constellation
        save_decoder_constellation(dbg, p.out_dir);
    end

    if p.pause_before_exit
        if p.pause_seconds < 0
            fprintf('\nPress any key to continue...\n');
            pause;
        else
            fprintf('\nPausing for %.2f seconds...\n', p.pause_seconds);
            pause(p.pause_seconds);
        end
    end
end

function [ber, bit_errors, bit_total] = compute_payload_ber(tx_payload, rx_payload)
    tx_payload = uint8(tx_payload(:));
    rx_payload = uint8(rx_payload(:));
    nbytes = min(numel(tx_payload), numel(rx_payload));
    if nbytes == 0
        ber = NaN;
        bit_errors = 0;
        bit_total = 0;
        return;
    end

    tx_use = tx_payload(1:nbytes);
    rx_use = rx_payload(1:nbytes);
    xor_bytes = bitxor(tx_use, rx_use);
    bit_errors = sum(arrayfun(@(x) sum(bitget(x, 1:8)), xor_bytes));
    bit_total = 8 * nbytes;
    ber = double(bit_errors) / double(bit_total);
end

function y = simulate_audio_channel(x, p)
    y = [zeros(p.timing_offset, 1); x];

    % Real multipath on the passband waveform.
    if ~isempty(p.channel_taps)
        h = real(p.channel_taps(:));
        y = filter(h, 1, y);
    end

    % Frequency error.
    n = (0:length(y)-1).';
    if p.disable_modulation
        y = y .* exp(1j * 2*pi * p.cfo_hz * n / p.fs);
    else
        ya = analytic_signal(y);
        y = real(ya .* exp(1j * 2*pi * p.cfo_hz * n / p.fs));
    end

    % Optional amplitude ripple.
    if p.apply_am_ripple
        y = y .* (1 + p.am_ripple_depth * sin(2*pi*p.am_ripple_hz*n/p.fs));
    end

    % AWGN.
    Ps = mean(abs(y).^2);
    Pn = Ps / (10^(p.snr_db/10));
    if isreal(y)
        y = y + sqrt(Pn) * randn(size(y));
    else
        y = y + sqrt(Pn/2) * (randn(size(y)) + 1j*randn(size(y)));
    end

    m = max(abs(y));
    if m > 0
        y = 0.85 * y / m;
    end
end

function xa = analytic_signal(x)
    x = x(:);
    N = length(x);
    X = fft(x);
    H = zeros(N,1);
    if mod(N,2) == 0
        H(1) = 1;
        H(N/2+1) = 1;
        H(2:N/2) = 2;
    else
        H(1) = 1;
        H(2:(N+1)/2) = 2;
    end
    xa = ifft(X .* H);
end

function fig_handles = make_plots(tx_audio, rx_audio, dbg, p)
    fig_handles = [];

    fig_handles(end+1) = figure('name', 'tx_waveform');
    plot(tx_audio);
    grid on;
    title('TX audio waveform');
    xlabel('Sample'); ylabel('Amplitude');

    fig_handles(end+1) = figure('name', 'rx_waveform');
    plot(rx_audio);
    grid on;
    title('RX audio waveform');
    xlabel('Sample'); ylabel('Amplitude');

    if isfield(dbg, 'sync_metric') && ~isempty(dbg.sync_metric)
        fig_handles(end+1) = figure('name', 'sync_metric');
        plot(dbg.sync_metric);
        grid on;
        title(sprintf('Sync metric (peak at %d, CFO est %.2f Hz)', dbg.peak_idx, dbg.cfo_est_hz));
        xlabel('Candidate start'); ylabel('Metric');
    end

    if isfield(dbg, 'rx_syms_raw') && ~isempty(dbg.rx_syms_raw)
        fig_handles(end+1) = figure('name', 'constellation_before_eq');
        plot(real(dbg.rx_syms_raw), imag(dbg.rx_syms_raw), '.');
        grid on; axis equal;
        xlabel('In-Phase'); ylabel('Quadrature');
        title('Constellation before equalization');
    end

    if isfield(dbg, 'rx_syms_eq') && ~isempty(dbg.rx_syms_eq)
        fig_handles(end+1) = figure('name', 'constellation_after_eq');
        plot(real(dbg.rx_syms_eq), imag(dbg.rx_syms_eq), 'x');
        grid on; axis equal;
        xlabel('In-Phase'); ylabel('Quadrature');
        title('Constellation after equalization');
        hold on;
        ideal = ideal_constellation(p.modulation);
        plot(real(ideal), imag(ideal), 'ro', 'markersize', 10, 'linewidth', 2);
        hold off;
    end
end

function ensure_out_dir(out_dir)
    if exist(out_dir, 'dir') ~= 7
        mkdir(out_dir);
    end
end

function save_plot_images(fig_handles, out_dir)
    for i = 1:numel(fig_handles)
        h = fig_handles(i);
        if ~ishandle(h)
            continue;
        end
        fig_name = get(h, 'name');
        if isempty(fig_name)
            fig_name = sprintf('figure_%02d', i);
        end
        png_path = fullfile(out_dir, [fig_name '.png']);
        saveas(h, png_path);
    end
    fprintf('Saved %d plot image(s) in: %s\n', numel(fig_handles), out_dir);
end

function save_decoder_constellation(dbg, out_dir)
    if isfield(dbg, 'rx_syms_eq') && ~isempty(dbg.rx_syms_eq)
        sym = dbg.rx_syms_eq(:); %#ok<NASGU>
        basename = 'decoder_constellation_eq';
    elseif isfield(dbg, 'rx_syms_raw') && ~isempty(dbg.rx_syms_raw)
        sym = dbg.rx_syms_raw(:); %#ok<NASGU>
        basename = 'decoder_constellation_raw';
    else
        fprintf('Decoder constellation not available: no symbols in dbg.\n');
        return;
    end

    mat_path = fullfile(out_dir, [basename '.mat']);
    save(mat_path, 'sym');

    csv_data = [real(sym), imag(sym)];
    csv_path = fullfile(out_dir, [basename '.csv']);
    csvwrite(csv_path, csv_data);
    fprintf('Saved decoder constellation in: %s and %s\n', mat_path, csv_path);
end

function ideal = ideal_constellation(modulation)
    switch upper(modulation)
        case 'BPSK'
            ideal = [-1; +1];
        case 'QPSK'
            ideal = [1+1j; 1-1j; -1+1j; -1-1j] / sqrt(2);
        otherwise
            ideal = [];
    end
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
    p.use_chirp_sync = true;
    p.sync_chirp_f0 = 4000;
    p.sync_chirp_f1 = 8000;
    p.sync_half_len = 64;
    p.packet_payload_bytes = 24;
    p.payload_bytes = 64;
    p.session_id = uint16(1234);
    p.detect_threshold = 0.005;
    p.wake_min_score = 0.005;
    p.wake_search_len = 8000;
    p.wake_miss_hop = round(0.012 * p.fs);
    p.wake_retry_hop = round(0.002 * p.fs);
    p.wake_rearm_hop = round(0.008 * p.fs);
    p.wake_ref_pre_ms = 15;
    p.sync_search_len = 1500;
    p.sync_peak_rel_threshold = 0.80;
    p.sync_metric_abs_threshold = 0.02;
    p.cfo_grid_hz = -12:2:12;
    p.rx_preroll = 128;
    p.pll_enable = true;
    p.pll_kp = 0.05;
    p.pll_ki = 0.002;

    p.snr_db = 24;
    p.cfo_hz = 4;
    p.timing_offset = 120;
    p.channel_taps = [1.0; 0.22; -0.08];
    p.apply_am_ripple = true;
    p.am_ripple_depth = 0.03;
    p.am_ripple_hz = 40;
    p.make_plots = true;
    p.pause_before_exit = true;
    p.pause_seconds = -1;
    p.save_images = true;
    p.save_decoder_constellation = true;
    p.out_dir = fullfile(pwd, '..', 'images');
    p.disable_modulation = false;
    p.verbose = true;
end
