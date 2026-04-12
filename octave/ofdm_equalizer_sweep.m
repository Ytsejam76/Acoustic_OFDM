% Copyright (c) 2026 Elias S. G. Carotti

function stats = ofdm_equalizer_sweep(varargin)
% OFDM_EQUALIZER_SWEEP  Monte Carlo comparison of equalizer modes.
%
%   stats = ofdm_equalizer_sweep()
%   stats = ofdm_equalizer_sweep(cfg)
%
% cfg fields (all optional):
%   equalizer_modes      : cell array of mode strings
%   num_trials           : trials per mode
%   oracle_sync          : if true, force oracle wake/sync/CFO
%   show_progress        : if true, print live progress
%   make_octave_plots    : if true, plot mode comparison
%   base_params          : struct merged into each trial
%   save_plot            : if true, save comparison plot and stats
%   out_dir              : output directory
%   plot_filename        : output image name
%
% Returned stats fields:
%   modes, success_count, decode_rate, per, ber_decoded, ber_effective,
%   bit_errors, bit_total, erasure_bits, total_tx_bits, num_trials, cfg

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

    stats = run_mode_sweep(cfg);
    print_summary(stats);
    if cfg.make_octave_plots
        plot_mode_comparison(stats);
    end

    if cfg.save_plot
        ensure_out_dir(cfg.out_dir);
        if cfg.make_octave_plots
            saveas(gcf, fullfile(cfg.out_dir, cfg.plot_filename));
        end
        save('-v7', fullfile(cfg.out_dir, 'equalizer_sweep_stats.mat'), 'stats');
    end
end

function stats = run_mode_sweep(cfg)
    modes = cfg.equalizer_modes(:);
    nmodes = numel(modes);
    success_count = zeros(nmodes, 1);
    decoded_count = zeros(nmodes, 1);
    bit_errors = zeros(nmodes, 1);
    bit_total = zeros(nmodes, 1);
    erasure_bits = zeros(nmodes, 1);
    total_tx_bits = zeros(nmodes, 1);

    total_trials = nmodes * cfg.num_trials;
    done_trials = 0;
    if cfg.show_progress
        fprintf('[eq-sweep] Starting: %d modes x %d trials\n', nmodes, cfg.num_trials);
    end

    for mi = 1:nmodes
        mode = char(modes{mi});
        if cfg.show_progress
            fprintf('[eq-sweep] Mode %s (%d/%d)\n', mode, mi, nmodes);
        end
        for t = 1:cfg.num_trials
            p = cfg.base_params;
            p.equalizer_mode = mode;
            p = apply_oracle_sync_if_needed(p, cfg.oracle_sync);
            p.pause_before_exit = false;
            p.make_plots = false;
            p.save_images = false;
            p.save_decoder_constellation = false;
            p.verbose = false;

            result = ofdm_test_channel(p);
            tx_bits_this = 8 * double(p.payload_bytes);
            success_count(mi) = success_count(mi) + double(result.success);
            total_tx_bits(mi) = total_tx_bits(mi) + tx_bits_this;
            if result.bit_total_compared > 0
                decoded_count(mi) = decoded_count(mi) + 1;
                bit_errors(mi) = bit_errors(mi) + double(result.bit_errors);
                bit_total(mi) = bit_total(mi) + double(result.bit_total_compared);
            end
            erasure_bits(mi) = erasure_bits(mi) + max(0, tx_bits_this - double(result.bit_total_compared));

            done_trials = done_trials + 1;
            if cfg.show_progress
                print_progress_line(done_trials, total_trials);
            end
        end
    end
    if cfg.show_progress
        fprintf('\n');
    end

    decode_rate = decoded_count / cfg.num_trials;
    per = 1 - (success_count / cfg.num_trials);
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
    stats.success_count = success_count;
    stats.decode_rate = decode_rate;
    stats.per = per;
    stats.ber_decoded = ber_decoded;
    stats.ber_effective = ber_effective;
    stats.bit_errors = bit_errors;
    stats.bit_total = bit_total;
    stats.erasure_bits = erasure_bits;
    stats.total_tx_bits = total_tx_bits;
    stats.num_trials = cfg.num_trials;
    stats.cfg = cfg;
end

function print_summary(stats)
    fprintf('\n==== OFDM EQUALIZER SWEEP ====\n');
    fprintf('Trials per mode: %d\n', stats.num_trials);
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
            char(stats.modes{i}), stats.success_count(i), stats.num_trials, ...
            stats.decode_rate(i), ber_str, ber_eff_str);
    end
end

function plot_mode_comparison(stats)
    nmodes = numel(stats.modes);
    x = 1:nmodes;

    figure('name', 'equalizer_sweep_compare');

    subplot(2,1,1);
    bar(x, [stats.success_count / stats.num_trials, stats.decode_rate], 'grouped');
    grid on;
    ylim([0 1]);
    ylabel('Rate');
    title('Equalizer mode success and decode rate');
    set(gca, 'xtick', x, 'xticklabel', stats.modes);
    xtickangle(20);
    legend('Success rate', 'Decode rate', 'location', 'northeast');

    subplot(2,1,2);
    vals = stats.ber_effective;
    vals(isnan(vals) | vals <= 0) = eps;
    semilogy(x, vals, '-o', 'linewidth', 1.5);
    grid on;
    ylabel('BER');
    title('Effective BER by equalizer mode');
    set(gca, 'xtick', x, 'xticklabel', stats.modes);
    xtickangle(20);
end

function p = apply_oracle_sync_if_needed(p, oracle_sync)
    if oracle_sync
        pre_sil = round(0.015 * p.fs);
        p.oracle_wake_start = p.timing_offset + pre_sil + 1;
        p.oracle_sync_start = 1;
        p.oracle_cfo_est_hz = p.cfo_hz;
    else
        p.oracle_wake_start = [];
        p.oracle_sync_start = [];
        p.oracle_cfo_est_hz = [];
    end
end

function cfg = default_cfg()
    cfg = struct();
    cfg.equalizer_modes = {'pilot-denoise', 'pilot-denoise-temporal', 'pilot-denoise-wiener'};
    cfg.num_trials = 20;
    cfg.oracle_sync = true;
    cfg.show_progress = true;
    cfg.make_octave_plots = true;
    cfg.base_params = default_base_params();
    cfg.save_plot = true;
    cfg.out_dir = fullfile(pwd, '..', 'images');
    cfg.plot_filename = 'equalizer_sweep.png';
end

function print_progress_line(done_trials, total_trials)
    width = 30;
    ratio = done_trials / max(1, total_trials);
    nfill = max(0, min(width, floor(width * ratio)));
    bar = [repmat('=', 1, nfill), repmat('-', 1, width - nfill)];
    fprintf('\r[eq-sweep] [%s] %5.1f%% (%d/%d)', bar, 100*ratio, done_trials, total_trials);
end

function ensure_out_dir(out_dir)
    if exist(out_dir, 'dir') ~= 7
        mkdir(out_dir);
    end
end

function p = default_base_params()
    p = struct();
    p.fs = 48000;
    p.fc = 17000;
    p.Nfft = 96;
    p.Ncp = 72;
    p.used_bins = [2 3 4 5];
    p.pilot_bins = [2 4];
    p.base_freq_hz = [];
    p.num_pilots = [];
    p.modulation = 'BPSK';
    p.use_pilots = [];
    p.equalizer_mode = 'training-pilot';
    p.temporal_window = 4;
    p.residual_tap_order_mode = 'all';
    p.residual_tap_order = 1;
    p.residual_tap_order_max = 7;
    p.wake_ms = 12;
    p.wake_freq = 16500;
    p.wake_guard_ms = 4;
    p.wake_miss_hop = round(0.012 * p.fs);
    p.wake_retry_hop = round(0.002 * p.fs);
    p.wake_rearm_hop = round(0.008 * p.fs);
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
    p.snr_db = 24;
    p.cfo_hz = 0;
    p.timing_offset = 0;
    p.channel_taps = 1;
    p.echo_profile = 'none';
    p.echo_delays_ms = [];
    p.echo_gains = [];
    p.echo_phases_deg = [];
    p.apply_am_ripple = false;
    p.am_ripple_depth = 0.03;
    p.am_ripple_hz = 40;
    p.disable_lpf = false;
    p.disable_modulation = false;
    p.pause_before_exit = false;
    p.pause_seconds = -1;
    p.make_plots = false;
    p.save_images = false;
    p.save_decoder_constellation = false;
    p.out_dir = fullfile(pwd, '..', 'images');
    p.verbose = false;
end
