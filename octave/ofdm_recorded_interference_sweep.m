% Copyright (c) 2026 Elias S. G. Carotti

function stats = ofdm_recorded_interference_sweep(varargin)
% OFDM_RECORDED_INTERFERENCE_SWEEP  Compare equalizers on fixed recording chunks.
%
%   stats = ofdm_recorded_interference_sweep()
%   stats = ofdm_recorded_interference_sweep(cfg)
%
% This helper keeps the experiment offline and reproducible by reusing the
% same recorded-noise chunks across all equalizer modes.

    cfg = default_cfg();
    if nargin >= 1 && ~isempty(varargin{1})
        user_cfg = varargin{1};
        f = fieldnames(user_cfg);
        for i = 1:numel(f)
            if strcmp(f{i}, 'base_params') && isstruct(user_cfg.base_params)
                bf = fieldnames(user_cfg.base_params);
                for bi = 1:numel(bf)
                    cfg.base_params.(bf{bi}) = user_cfg.base_params.(bf{bi});
                end
            else
                cfg.(f{i}) = user_cfg.(f{i});
            end
        end
    end

    stats = run_recording_sweep(cfg);
    print_summary(stats);
end

function stats = run_recording_sweep(cfg)
    modes = cfg.equalizer_modes(:);
    nmodes = numel(modes);
    nchunks = cfg.num_chunks;
    success_count = zeros(nmodes, 1);
    decoded_count = zeros(nmodes, 1);
    bit_errors = zeros(nmodes, 1);
    bit_total = zeros(nmodes, 1);
    erasure_bits = zeros(nmodes, 1);
    total_tx_bits = zeros(nmodes, 1);
    success_by_chunk = false(nchunks, nmodes);

    total_runs = nmodes * nchunks;
    done_runs = 0;
    if cfg.show_progress
        fprintf('[rec-sweep] Starting: %d chunks x %d modes\n', nchunks, nmodes);
    end

    for ci = 1:nchunks
        payload = chunk_payload(cfg, ci);
        offset_samples = chunk_offset_samples(cfg, ci);
        for mi = 1:nmodes
            p = cfg.base_params;
            p.equalizer_mode = char(modes{mi});
            p.interference_file = cfg.interference_file;
            p.interference_offset_samples = offset_samples;
            p.interference_ratio_db = cfg.interference_ratio_db;
            p.pause_before_exit = false;
            p.make_plots = false;
            p.save_images = false;
            p.save_decoder_constellation = false;
            p.verbose = false;

            result = ofdm_test_channel(p, payload);
            tx_bits_this = 8 * double(numel(payload));
            success_count(mi) = success_count(mi) + double(result.success);
            success_by_chunk(ci, mi) = logical(result.success);
            total_tx_bits(mi) = total_tx_bits(mi) + tx_bits_this;
            if result.bit_total_compared > 0
                decoded_count(mi) = decoded_count(mi) + 1;
                bit_errors(mi) = bit_errors(mi) + double(result.bit_errors);
                bit_total(mi) = bit_total(mi) + double(result.bit_total_compared);
            end
            erasure_bits(mi) = erasure_bits(mi) + max(0, tx_bits_this - double(result.bit_total_compared));

            done_runs = done_runs + 1;
            if cfg.show_progress
                print_progress_line(done_runs, total_runs, ci, nchunks, char(modes{mi}));
            end
        end
    end
    if cfg.show_progress
        fprintf('\n');
    end

    decode_rate = decoded_count / nchunks;
    ber_decoded = NaN(nmodes, 1);
    ber_effective = NaN(nmodes, 1);
    for mi = 1:nmodes
        if bit_total(mi) > 0
            ber_decoded(mi) = bit_errors(mi) / bit_total(mi);
        end
        if total_tx_bits(mi) > 0
            ber_effective(mi) = (bit_errors(mi) + erasure_bits(mi)) / total_tx_bits(mi);
        end
    end

    stats = struct();
    stats.modes = modes;
    stats.num_chunks = nchunks;
    stats.success_count = success_count;
    stats.decode_rate = decode_rate;
    stats.ber_decoded = ber_decoded;
    stats.ber_effective = ber_effective;
    stats.bit_errors = bit_errors;
    stats.bit_total = bit_total;
    stats.erasure_bits = erasure_bits;
    stats.total_tx_bits = total_tx_bits;
    stats.success_by_chunk = success_by_chunk;
    stats.chunk_offsets_samples = arrayfun(@(i) chunk_offset_samples(cfg, i), (1:nchunks).');
    stats.cfg = cfg;
end

function payload = chunk_payload(cfg, chunk_idx)
    if isfield(cfg, 'payload_seed') && ~isempty(cfg.payload_seed)
        rand('seed', double(cfg.payload_seed) + chunk_idx - 1);
    end
    payload = uint8(randi([0 255], cfg.base_params.payload_bytes, 1));
end

function offset_samples = chunk_offset_samples(cfg, chunk_idx)
    offset_samples = cfg.start_offset_samples + (chunk_idx - 1) * cfg.chunk_advance_samples;
end

function print_summary(stats)
    fprintf('\n==== RECORDED INTERFERENCE SWEEP ====\n');
    fprintf('Chunks per mode: %d\n', stats.num_chunks);
    fprintf('Mode                     Success   DecodeRate   BER(decoded)   BER(effective)\n');
    for i = 1:numel(stats.modes)
        if isnan(stats.ber_decoded(i))
            ber_str = 'n/a';
        else
            ber_str = sprintf('%.3e', stats.ber_decoded(i));
        end
        if isnan(stats.ber_effective(i))
            ber_eff_str = 'n/a';
        else
            ber_eff_str = sprintf('%.3e', stats.ber_effective(i));
        end
        fprintf('%-24s %3d/%-3d   %10.3f   %s   %s\n', ...
            char(stats.modes{i}), stats.success_count(i), stats.num_chunks, ...
            stats.decode_rate(i), ber_str, ber_eff_str);
    end
end

function print_progress_line(done_runs, total_runs, chunk_idx, total_chunks, mode)
    width = 30;
    ratio = done_runs / max(1, total_runs);
    nfill = max(0, min(width, floor(width * ratio)));
    bar = [repmat('=', 1, nfill), repmat('-', 1, width - nfill)];
    fprintf('\r[rec-sweep] [%s] %5.1f%% (%d/%d) chunk=%d/%d mode=%s', ...
        bar, 100 * ratio, done_runs, total_runs, chunk_idx, total_chunks, mode);
end

function cfg = default_cfg()
    cfg = struct();
    cfg.equalizer_modes = {'pilot-denoise', 'pilot-denoise-temporal', 'pilot-denoise-wiener-psd'};
    cfg.num_chunks = 12;
    cfg.interference_file = fullfile('octave', 'noise_recordings', 'interference.wav');
    cfg.interference_ratio_db = 12;
    cfg.start_offset_samples = 0;
    cfg.chunk_advance_samples = 12000;
    cfg.payload_seed = 12345;
    cfg.show_progress = true;
    cfg.base_params = default_base_params();
end

function p = default_base_params()
    p = struct();
    p.fs = 48000;
    p.fc = 17000;
    p.Nfft = 96;
    p.Ncp = 72;
    p.used_bins = [2 3 4 5 6 7 8 9];
    p.pilot_bins = [2 4 7 9];
    p.base_freq_hz = [];
    p.num_pilots = [];
    p.modulation = 'QPSK';
    p.use_pilots = true;
    p.equalizer_mode = 'training-pilot';
    p.temporal_window = 4;
    p.disturbance_temporal_alpha = 0.75;
    p.disturbance_freq_smooth = 5;
    p.disturbance_psd_gain = 1.0;
    p.residual_tap_order_mode = 'all';
    p.residual_tap_order = 1;
    p.residual_tap_order_max = 7;
    p.wake_ms = 12;
    p.wake_freq = 16500;
    p.wake_guard_ms = 4;
    p.use_chirp_sync = true;
    p.sync_chirp_f0 = 4000;
    p.sync_chirp_f1 = 8000;
    p.sync_half_len = 64;
    p.packet_payload_bytes = 24;
    p.payload_bytes = 24;
    p.session_id = uint16(1234);
    p.detect_threshold = 0.005;
    p.sync_search_len = 1500;
    p.cfo_grid_hz = -12:2:12;
    p.rx_preroll = 128;
    p.pll_enable = true;
    p.pll_kp = 0.05;
    p.pll_ki = 0.002;
    p.snr_db = 80;
    p.cfo_hz = 4;
    p.timing_offset = 120;
    p.channel_taps = [1.0; 0.22; -0.08];
    p.echo_delays_ms = [];
    p.echo_gains = [];
    p.echo_phases_deg = [];
    p.apply_am_ripple = true;
    p.am_ripple_depth = 0.03;
    p.am_ripple_hz = 40;
    p.oracle_wake_start = [];
    p.oracle_sync_start = [];
    p.oracle_cfo_est_hz = [];
end
