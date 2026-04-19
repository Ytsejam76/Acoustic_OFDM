% Copyright (c) 2026 Elias S. G. Carotti

function stats = ofdm_recorded_interference_compare(varargin)
% OFDM_RECORDED_INTERFERENCE_COMPARE  Compare equalizers on shared noise chunks.
%
%   stats = ofdm_recorded_interference_compare()
%   stats = ofdm_recorded_interference_compare(cfg)
%
% This helper is intended for exploratory visual comparison. Each segment
% reuses the same payload and the same interference chunk across all modes,
% and saves one overlaid constellation plot per segment plus one compact
% Wiener PSD diagnostic plot when the Wiener PSD mode is present.

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

    ensure_out_dir(cfg.out_dir);
    stats = run_compare(cfg);
    write_summary(stats);
    print_summary(stats);
end

function stats = run_compare(cfg)
    modes = cfg.equalizer_modes(:);
    nmodes = numel(modes);
    nsegments = cfg.num_segments;
    total_runs = nmodes * nsegments;
    done_runs = 0;

    rows = struct([]);
    row_idx = 0;

    if cfg.show_progress
        fprintf('[rec-compare] Starting: %d segments x %d modes\n', nsegments, nmodes);
    end

    for si = 1:nsegments
        payload = segment_payload(cfg, si);
        offset_samples = segment_offset_samples(cfg, si);
        segment_rows = struct([]);
        segment_dbg = cell(nmodes, 1);

        for mi = 1:nmodes
            mode = char(modes{mi});
            p = cfg.base_params;
            p.equalizer_mode = mode;
            p.interference_file = cfg.interference_file;
            p.interference_offset_samples = offset_samples;
            p.interference_ratio_db = cfg.interference_ratio_db;
            p.pause_before_exit = false;
            p.make_plots = false;
            p.save_images = false;
            p.save_decoder_constellation = false;
            p.verbose = false;

            [result, dbg] = ofdm_test_channel(p, payload);
            segment_dbg{mi} = dbg;

            row_idx = row_idx + 1;
            rows(row_idx).segment = si;
            rows(row_idx).offset_samples = offset_samples;
            rows(row_idx).mode = mode;
            rows(row_idx).success = logical(result.success);
            rows(row_idx).ber = result.ber;
            rows(row_idx).bit_errors = double(result.bit_errors);
            rows(row_idx).bit_total = double(result.bit_total_compared);
            rows(row_idx).plot_path = segment_plot_path(cfg.out_dir, si);
            rows(row_idx).diag_path = segment_diag_path(cfg.out_dir, si);
            segment_rows(mi) = rows(row_idx); %#ok<AGROW>

            done_runs = done_runs + 1;
            if cfg.show_progress
                print_progress_line(done_runs, total_runs, si, nsegments, mode, result.success);
            end
        end

        if cfg.save_images || cfg.make_plots
            save_segment_compare_plot(cfg, si, offset_samples, segment_rows, segment_dbg);
            save_segment_diagnostic_plot(cfg, si, offset_samples, segment_rows, segment_dbg);
        end
    end
    if cfg.show_progress
        fprintf('\n');
    end

    stats = struct();
    stats.rows = rows;
    stats.cfg = cfg;
end

function payload = segment_payload(cfg, segment_idx)
    if isfield(cfg, 'payload_seed') && ~isempty(cfg.payload_seed)
        rand('seed', double(cfg.payload_seed) + segment_idx - 1);
    end
    payload = uint8(randi([0 255], cfg.base_params.payload_bytes, 1));
end

function offset_samples = segment_offset_samples(cfg, segment_idx)
    offset_samples = cfg.start_offset_samples + (segment_idx - 1) * cfg.chunk_advance_samples;
end

function path = segment_plot_path(root_dir, segment_idx)
    path = fullfile(root_dir, sprintf('segment_%02d_constellation_compare.png', segment_idx));
end

function path = segment_diag_path(root_dir, segment_idx)
    path = fullfile(root_dir, sprintf('segment_%02d_wiener_psd_diag.png', segment_idx));
end

function save_segment_compare_plot(cfg, segment_idx, offset_samples, segment_rows, segment_dbg)
    colors = lines(numel(segment_rows));
    h = figure('visible', ternary_visible(cfg.make_plots), ...
        'name', sprintf('segment_%02d_compare', segment_idx), ...
        'position', [100 100 960 720]);
    hold on;

    labels = cell(numel(segment_rows), 1);
    all_pts = [];
    for i = 1:numel(segment_rows)
        dbg = segment_dbg{i};
        if isfield(dbg, 'rx_syms_eq') && ~isempty(dbg.rx_syms_eq)
            plot(real(dbg.rx_syms_eq), imag(dbg.rx_syms_eq), '.', ...
                'color', colors(i, :), 'markersize', 10);
            all_pts = [all_pts; dbg.rx_syms_eq(:)]; %#ok<AGROW>
        end
        labels{i} = sprintf('%s [ok=%d]', short_mode_name(segment_rows(i).mode), segment_rows(i).success);
    end

    ideal = ideal_constellation_local(cfg.base_params.modulation);
    if ~isempty(ideal)
        plot(real(ideal), imag(ideal), 'ko', 'markersize', 10, 'linewidth', 1.5);
        labels{end+1} = 'ideal';
    end

    hold off;
    grid on;
    axis equal;
    apply_constellation_limits(all_pts, ideal);
    xlabel('In-Phase');
    ylabel('Quadrature');
    title(sprintf('Segment %d | offset %d | interference %.1f dB', ...
        segment_idx, offset_samples, cfg.interference_ratio_db), 'fontsize', 13);
    legend(labels, 'location', 'southoutside', 'orientation', 'horizontal');

    if cfg.save_images
        save_figure_png(h, segment_plot_path(cfg.out_dir, segment_idx));
    end
    if ~cfg.make_plots
        close(h);
    end
end

function save_segment_diagnostic_plot(cfg, segment_idx, offset_samples, segment_rows, segment_dbg)
    widx = find_wiener_psd_mode(segment_rows);
    if isempty(widx)
        return;
    end

    dbg = segment_dbg{widx};
    if ~isfield(dbg, 'disturbance_psd') || isempty(fieldnames(dbg.disturbance_psd))
        return;
    end

    h = figure('visible', ternary_visible(cfg.make_plots), ...
        'name', sprintf('segment_%02d_wiener_diag', segment_idx), ...
        'position', [120 120 960 900]);

    subplot(3,1,1);
    if isfield(dbg.disturbance_psd, 'noise_bins') && ~isempty(dbg.disturbance_psd.noise_bins) ...
            && isfield(dbg.disturbance_psd, 'noise_rx_raw') && ~isempty(dbg.disturbance_psd.noise_rx_raw)
        plot(dbg.disturbance_psd.noise_bins, dbg.disturbance_psd.noise_rx_raw, '.-');
        grid on;
        xlabel('FFT bin');
        ylabel('Power');
        title(sprintf('Segment %d, offset %d: unused-bin disturbance power', segment_idx, offset_samples));
    else
        axis off;
        text(0.1, 0.5, 'No unused-bin disturbance debug available');
    end

    subplot(3,1,2);
    if isfield(dbg.disturbance_psd, 'noise_rx_interp_used') && ~isempty(dbg.disturbance_psd.noise_rx_interp_used)
        plot(cfg.base_params.used_bins, dbg.disturbance_psd.noise_rx_interp_used, 'o-');
        grid on;
        xlabel('Used FFT bin');
        ylabel('Power');
        title('Interpolated disturbance PSD on used bins');
    else
        axis off;
        text(0.1, 0.5, 'No interpolated PSD debug available');
    end

    subplot(3,1,3);
    if isfield(dbg, 'wiener_gain_used') && ~isempty(dbg.wiener_gain_used)
        stem(cfg.base_params.used_bins, dbg.wiener_gain_used, 'filled');
        grid on;
        ylim([0 1.05]);
        xlabel('Used FFT bin');
        ylabel('Gain');
        title(sprintf('Wiener gain for %s', segment_rows(widx).mode));
    else
        axis off;
        text(0.1, 0.5, 'No Wiener gain debug available');
    end

    if cfg.save_images
        save_figure_png(h, segment_diag_path(cfg.out_dir, segment_idx));
    end
    if ~cfg.make_plots
        close(h);
    end
end

function idx = find_wiener_psd_mode(segment_rows)
    idx = [];
    for i = 1:numel(segment_rows)
        if strcmpi(segment_rows(i).mode, 'pilot-denoise-wiener-psd')
            idx = i;
            return;
        end
    end
end

function visible = ternary_visible(make_plots)
    if make_plots
        visible = 'on';
    else
        visible = 'off';
    end
end

function name = short_mode_name(mode)
    switch lower(strtrim(mode))
        case 'pilot-denoise'
            name = 'PD';
        case 'pilot-denoise-temporal'
            name = 'PDT';
        case 'pilot-denoise-wiener-psd'
            name = 'WPSD';
        otherwise
            name = mode;
    end
end

function apply_constellation_limits(all_pts, ideal)
    pts = all_pts(:);
    if nargin >= 2 && ~isempty(ideal)
        pts = [pts; ideal(:)];
    end
    if isempty(pts)
        return;
    end
    xr = [min(real(pts)), max(real(pts))];
    yr = [min(imag(pts)), max(imag(pts))];
    span = max([xr(2) - xr(1), yr(2) - yr(1), 1.5]);
    cx = mean(xr);
    cy = mean(yr);
    half = 0.6 * span;
    xlim([cx - half, cx + half]);
    ylim([cy - half, cy + half]);
end

function save_figure_png(h, path)
    set(h, 'paperpositionmode', 'auto');
    print(h, path, '-dpng', '-r140');
end

function ideal = ideal_constellation_local(modulation)
    switch upper(modulation)
        case 'BPSK'
            ideal = [-1; +1];
        case 'QPSK'
            ideal = [1+1j; 1-1j; -1+1j; -1-1j] / sqrt(2);
        otherwise
            ideal = [];
    end
end

function write_summary(stats)
    summary_path = fullfile(stats.cfg.out_dir, 'summary.tsv');
    fid = fopen(summary_path, 'w');
    if fid < 0
        error('Failed to open summary file: %s', summary_path);
    end
    cleaner = onCleanup(@() fclose(fid));
    fprintf(fid, 'segment\toffset_samples\tmode\tsuccess\tber\tbit_errors\tbit_total\tplot_path\tdiag_path\n');
    for i = 1:numel(stats.rows)
        row = stats.rows(i);
        if isnan(row.ber)
            ber_str = 'nan';
        else
            ber_str = sprintf('%.6g', row.ber);
        end
        fprintf(fid, '%d\t%d\t%s\t%d\t%s\t%d\t%d\t%s\t%s\n', ...
            row.segment, row.offset_samples, row.mode, row.success, ...
            ber_str, row.bit_errors, row.bit_total, row.plot_path, row.diag_path);
    end
end

function print_summary(stats)
    fprintf('\n==== RECORDED INTERFERENCE COMPARE ====\n');
    fprintf('Output root: %s\n', stats.cfg.out_dir);
    fprintf('segment\toffset_samples\tmode\tsuccess\tber\tplot_path\tdiag_path\n');
    for i = 1:numel(stats.rows)
        row = stats.rows(i);
        if isnan(row.ber)
            ber_str = 'n/a';
        else
            ber_str = sprintf('%.3e', row.ber);
        end
        fprintf('%d\t%d\t%s\t%d\t%s\t%s\t%s\n', ...
            row.segment, row.offset_samples, row.mode, row.success, ...
            ber_str, row.plot_path, row.diag_path);
    end
end

function print_progress_line(done_runs, total_runs, segment_idx, total_segments, mode, success)
    width = 30;
    ratio = done_runs / max(1, total_runs);
    nfill = max(0, min(width, floor(width * ratio)));
    bar = [repmat('=', 1, nfill), repmat('-', 1, width - nfill)];
    fprintf('\r[rec-compare] [%s] %5.1f%% (%d/%d) segment=%d/%d mode=%s success=%d', ...
        bar, 100 * ratio, done_runs, total_runs, segment_idx, total_segments, mode, success);
end

function cfg = default_cfg()
    cfg = struct();
    cfg.equalizer_modes = {'pilot-denoise', 'pilot-denoise-temporal', 'pilot-denoise-wiener-psd'};
    cfg.num_segments = 3;
    cfg.interference_file = fullfile('octave', 'noise_recordings', 'interference.wav');
    cfg.interference_ratio_db = 12;
    cfg.start_offset_samples = 0;
    cfg.chunk_advance_samples = 12000;
    cfg.payload_seed = 12345;
    cfg.show_progress = true;
    cfg.make_plots = false;
    cfg.save_images = true;
    cfg.out_dir = fullfile('output', 'octave_recorded_compare');
    cfg.base_params = default_base_params();
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
