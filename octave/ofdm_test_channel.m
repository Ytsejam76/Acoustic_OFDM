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

    if p.make_plots
        make_plots(tx_audio, rx_audio, dbg, p);
    end
end

function y = simulate_audio_channel(x, p)
    y = [zeros(p.timing_offset, 1); x];

    % Real multipath on the passband waveform.
    if ~isempty(p.channel_taps)
        h = real(p.channel_taps(:));
        y = filter(h, 1, y);
    end

    % Small passband frequency error.
    n = (0:length(y)-1).';
    y = real(y .* exp(1j * 2*pi * p.cfo_hz * n / p.fs));

    % Optional amplitude ripple.
    if p.apply_am_ripple
        y = y .* (1 + p.am_ripple_depth * sin(2*pi*p.am_ripple_hz*n/p.fs));
    end

    % AWGN.
    Ps = mean(y.^2);
    Pn = Ps / (10^(p.snr_db/10));
    y = y + sqrt(Pn) * randn(size(y));

    m = max(abs(y));
    if m > 0
        y = 0.85 * y / m;
    end
end

function make_plots(tx_audio, rx_audio, dbg, p)
    figure;
    plot(tx_audio);
    grid on;
    title('TX audio waveform');
    xlabel('Sample'); ylabel('Amplitude');

    figure;
    plot(rx_audio);
    grid on;
    title('RX audio waveform');
    xlabel('Sample'); ylabel('Amplitude');

    if isfield(dbg, 'sync_metric') && ~isempty(dbg.sync_metric)
        figure;
        plot(dbg.sync_metric);
        grid on;
        title(sprintf('Sync metric (peak at %d, CFO est %.2f Hz)', dbg.peak_idx, dbg.cfo_est_hz));
        xlabel('Candidate start'); ylabel('Metric');
    end

    if isfield(dbg, 'rx_syms_raw') && ~isempty(dbg.rx_syms_raw)
        figure;
        plot(real(dbg.rx_syms_raw), imag(dbg.rx_syms_raw), '.');
        grid on; axis equal;
        xlabel('In-Phase'); ylabel('Quadrature');
        title('Constellation before equalization');
    end

    if isfield(dbg, 'rx_syms_eq') && ~isempty(dbg.rx_syms_eq)
        figure;
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
    p.Nfft = 32;
    p.Ncp = 8;
    p.used_bins = [3 4 5 6 7 8 9 10];
    p.modulation = 'BPSK';
    p.wake_ms = 12;
    p.wake_freq = 16500;
    p.wake_guard_ms = 4;
    p.sync_half_len = 32;
    p.packet_payload_bytes = 24;
    p.payload_bytes = 64;
    p.session_id = uint16(1234);
    p.detect_threshold = 0.35;
    p.sync_search_len = 500;

    p.snr_db = 24;
    p.cfo_hz = 4;
    p.timing_offset = 120;
    p.channel_taps = [1.0; 0.22; -0.08];
    p.apply_am_ripple = true;
    p.am_ripple_depth = 0.03;
    p.am_ripple_hz = 40;
    p.make_plots = true;
end
