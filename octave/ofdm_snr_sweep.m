% Copyright (c) 2026 Elias S. G. Carotti

function stats = ofdm_snr_sweep(varargin)
% OFDM_SNR_SWEEP  Monte Carlo sweep of OFDM performance vs SNR.
%
%   stats = ofdm_snr_sweep()
%   stats = ofdm_snr_sweep(cfg)
%
% cfg fields (all optional):
%   snr_db_list           : vector of SNR values in dB
%   num_trials            : trials per SNR point
%   oracle_sync           : if true, force oracle wake/sync/CFO
%   compare_oracle        : if true, run both oracle and estimated sync
%   show_progress         : if true, print live progress bar
%   base_params           : struct merged into each test run
%   save_plot             : if true, save summary plot
%   out_dir               : output directory for plot and .mat stats
%   plot_filename         : output image name
%
% Returned stats fields:
%   snr_db, per, decode_rate, ber_decoded, success_count, decoded_count,
%   bit_errors, bit_total, num_trials, cfg

    cfg = default_cfg();
    if nargin >= 1 && ~isempty(varargin{1})
        user_cfg = varargin{1};
        f = fieldnames(user_cfg);
        for i = 1:numel(f)
            cfg.(f{i}) = user_cfg.(f{i});
        end
    end

    if cfg.compare_oracle
        stats_oracle = run_one_sweep(cfg, true, 'Oracle sync');
        stats_no_oracle = run_one_sweep(cfg, false, 'Estimated sync');

        stats = struct();
        stats.oracle = stats_oracle;
        stats.nonoracle = stats_no_oracle;
        stats.cfg = cfg;

        print_summary(stats.oracle, 'Oracle sync');
        print_summary(stats.nonoracle, 'Estimated sync');
        plot_sweep_comparison(stats.oracle, stats.nonoracle);
    else
        stats = run_one_sweep(cfg, cfg.oracle_sync, ternary(cfg.oracle_sync, 'Oracle sync', 'Estimated sync'));
        print_summary(stats, ternary(cfg.oracle_sync, 'Oracle sync', 'Estimated sync'));
        plot_sweep(stats);
    end

    if cfg.save_plot
        ensure_out_dir(cfg.out_dir);
        fig_path = fullfile(cfg.out_dir, cfg.plot_filename);
        saveas(gcf, fig_path);
        mat_path = fullfile(cfg.out_dir, 'snr_sweep_stats.mat');
        save('-v7', mat_path, 'stats');
        fprintf('Saved sweep plot to: %s\n', fig_path);
        fprintf('Saved sweep stats to: %s\n', mat_path);
    end
end

function stats = run_one_sweep(cfg, oracle_sync, mode_label)
    snr_db = cfg.snr_db_list(:).';
    nsnr = numel(snr_db);

    success_count = zeros(nsnr, 1);
    decoded_count = zeros(nsnr, 1);
    bit_errors = zeros(nsnr, 1);
    bit_total = zeros(nsnr, 1);

    total_trials = nsnr * cfg.num_trials;
    done_trials = 0;
    if cfg.show_progress
        fprintf('[%s] Starting sweep: %d SNR points x %d trials\n', ...
            mode_label, nsnr, cfg.num_trials);
    end

    for si = 1:nsnr
        snr = snr_db(si);
        if cfg.show_progress
            fprintf('[%s] SNR %.1f dB (%d/%d)\n', mode_label, snr, si, nsnr);
        end
        for t = 1:cfg.num_trials
            p = cfg.base_params;
            p.snr_db = snr;
            p.pause_before_exit = false;
            p.make_plots = false;
            p.save_images = false;
            p.save_decoder_constellation = false;
            p.verbose = false;

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

            result = ofdm_test_channel(p);
            success_count(si) = success_count(si) + double(result.success);
            if result.bit_total_compared > 0
                decoded_count(si) = decoded_count(si) + 1;
                bit_errors(si) = bit_errors(si) + double(result.bit_errors);
                bit_total(si) = bit_total(si) + double(result.bit_total_compared);
            end

            done_trials = done_trials + 1;
            if cfg.show_progress
                print_progress_line(mode_label, done_trials, total_trials);
            end
        end
    end
    if cfg.show_progress
        fprintf('\n');
    end

    per = 1 - (success_count / cfg.num_trials);
    decode_rate = decoded_count / cfg.num_trials;
    ber_decoded = NaN(nsnr, 1);
    for si = 1:nsnr
        if bit_total(si) > 0
            ber_decoded(si) = bit_errors(si) / bit_total(si);
        end
    end

    stats = struct();
    stats.snr_db = snr_db(:);
    stats.per = per;
    stats.decode_rate = decode_rate;
    stats.ber_decoded = ber_decoded;
    stats.success_count = success_count;
    stats.decoded_count = decoded_count;
    stats.bit_errors = bit_errors;
    stats.bit_total = bit_total;
    stats.num_trials = cfg.num_trials;
    stats.oracle_sync = oracle_sync;
end

function print_summary(stats, label)
    fprintf('\n==== OFDM SNR SWEEP ====\n');
    fprintf('Mode: %s\n', label);
    fprintf('Trials per SNR: %d\n', stats.num_trials);
    fprintf('SNR(dB)   PER      DecodeRate   BER(decoded)\n');
    for i = 1:numel(stats.snr_db)
        if isnan(stats.ber_decoded(i))
            ber_str = 'n/a';
        else
            ber_str = sprintf('%.3e', stats.ber_decoded(i));
        end
        fprintf('%7.1f   %6.3f   %10.3f   %s\n', ...
            stats.snr_db(i), stats.per(i), stats.decode_rate(i), ber_str);
    end
end

function plot_sweep(stats)
    figure('name', 'snr_sweep');

    subplot(2,1,1);
    plot(stats.snr_db, stats.per, '-o', 'linewidth', 1.5);
    hold on;
    plot(stats.snr_db, 1 - stats.decode_rate, '--s', 'linewidth', 1.2);
    hold off;
    grid on;
    xlabel('SNR (dB)');
    ylabel('Rate');
    title('Packet Error and Erasure Rates vs SNR');
    legend('PER', '1 - Decode Rate', 'location', 'northeast');
    ylim([0 1]);

    subplot(2,1,2);
    idx = ~isnan(stats.ber_decoded);
    if any(idx)
        semilogy(stats.snr_db(idx), stats.ber_decoded(idx), '-x', 'linewidth', 1.5);
    else
        plot(stats.snr_db, nan(size(stats.snr_db)));
    end
    grid on;
    xlabel('SNR (dB)');
    ylabel('BER');
    title('BER on Decoded Bits vs SNR');
end

function plot_sweep_comparison(stats_oracle, stats_no_oracle)
    figure('name', 'snr_sweep_compare');

    subplot(2,1,1);
    plot(stats_oracle.snr_db, stats_oracle.per, '-o', 'linewidth', 1.5);
    hold on;
    plot(stats_no_oracle.snr_db, stats_no_oracle.per, '--s', 'linewidth', 1.5);
    hold off;
    grid on;
    xlabel('SNR (dB)');
    ylabel('PER');
    title('Packet Error Rate vs SNR');
    legend('Oracle sync', 'Estimated sync', 'location', 'northeast');
    ylim([0 1]);

    subplot(2,1,2);
    h = [];
    labels = {};
    idx1 = ~isnan(stats_oracle.ber_decoded);
    idx2 = ~isnan(stats_no_oracle.ber_decoded);
    if any(idx1)
        h(end+1) = semilogy(stats_oracle.snr_db(idx1), stats_oracle.ber_decoded(idx1), '-x', 'linewidth', 1.5); %#ok<AGROW>
        labels{end+1} = 'Oracle sync'; %#ok<AGROW>
        hold on;
    else
        hold on;
    end
    if any(idx2)
        h(end+1) = semilogy(stats_no_oracle.snr_db(idx2), stats_no_oracle.ber_decoded(idx2), '--d', 'linewidth', 1.5); %#ok<AGROW>
        labels{end+1} = 'Estimated sync'; %#ok<AGROW>
    end
    hold off;
    grid on;
    xlabel('SNR (dB)');
    ylabel('BER');
    title('BER on Decoded Bits vs SNR');
    if ~isempty(h)
        legend(h, labels, 'location', 'southwest');
    end
end

function ensure_out_dir(out_dir)
    if exist(out_dir, 'dir') ~= 7
        mkdir(out_dir);
    end
end

function cfg = default_cfg()
    cfg = struct();
    cfg.snr_db_list = 0:3:30;
    cfg.num_trials = 50;
    cfg.oracle_sync = true;
    cfg.compare_oracle = false;
    cfg.show_progress = true;
    cfg.base_params = default_base_params();
    cfg.save_plot = true;
    cfg.out_dir = fullfile(pwd, '..', 'images');
    cfg.plot_filename = 'snr_sweep.png';
end

function print_progress_line(mode_label, done_trials, total_trials)
    width = 30;
    ratio = done_trials / max(1, total_trials);
    nfill = max(0, min(width, floor(width * ratio)));
    bar = [repmat('=', 1, nfill), repmat('-', 1, width - nfill)];
    fprintf('\r[%s] [%s] %5.1f%% (%d/%d)', ...
        mode_label, bar, 100*ratio, done_trials, total_trials);
end

function y = ternary(cond, a, b)
    if cond
        y = a;
    else
        y = b;
    end
end

function p = default_base_params()
    p = struct();
    p.fs = 48000;
    p.fc = 17000;
    p.Nfft = 96;
    p.Ncp = 72;
    p.used_bins = [2 3 4 5];
    p.modulation = 'BPSK';
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
