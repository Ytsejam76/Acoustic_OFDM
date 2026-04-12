% Copyright (c) 2026 Elias S. G. Carotti

function stats = ofdm_wiener_probe(varargin)
% OFDM_WIENER_PROBE  Run the Wiener equalizer diagnostics across multiple taps.
%
%   stats = ofdm_wiener_probe()
%   stats = ofdm_wiener_probe(cfg)
%
% cfg fields:
%   base_params : parameter struct passed to ofdm_test_channel
%   tap_list    : vector of tap indices to probe

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

    tap_list = cfg.tap_list(:);
    ntaps = numel(tap_list);
    stats = struct();
    stats.tap = tap_list;
    stats.success = zeros(ntaps, 1);
    stats.has_dbg = zeros(ntaps, 1);
    stats.window_used = zeros(ntaps, 1);
    stats.sigma2 = NaN(ntaps, 1);
    stats.ryy0 = NaN(ntaps, 1);
    stats.rhh0 = NaN(ntaps, 1);
    stats.weight0 = NaN(ntaps, 1);
    stats.weight_sum = NaN(ntaps, 1);

    for i = 1:ntaps
        p = cfg.base_params;
        p.equalizer_mode = 'pilot-denoise-wiener';
        p.wiener_debug_tap = tap_list(i);
        p.pause_before_exit = false;
        p.make_plots = false;
        p.save_images = false;
        p.save_decoder_constellation = false;
        p.verbose = false;

        [result, dbg] = ofdm_test_channel(p);
        stats.success(i) = double(result.success);
        if isfield(dbg, 'wiener_dbg') && ~isempty(fieldnames(dbg.wiener_dbg))
            stats.has_dbg(i) = 1;
            if isfield(dbg, 'wiener_window_used')
                stats.window_used(i) = dbg.wiener_window_used;
            end
            if isfield(dbg.wiener_dbg, 'sigma2')
                stats.sigma2(i) = dbg.wiener_dbg.sigma2;
            end
            if isfield(dbg.wiener_dbg, 'ryy') && ~isempty(dbg.wiener_dbg.ryy)
                stats.ryy0(i) = real(dbg.wiener_dbg.ryy(1));
            end
            if isfield(dbg.wiener_dbg, 'rhh') && ~isempty(dbg.wiener_dbg.rhh)
                stats.rhh0(i) = real(dbg.wiener_dbg.rhh(1));
            end
            if isfield(dbg.wiener_dbg, 'weights') && ~isempty(dbg.wiener_dbg.weights)
                stats.weight0(i) = real(dbg.wiener_dbg.weights(1));
                stats.weight_sum(i) = sum(real(dbg.wiener_dbg.weights));
            end
        end
    end

    print_summary(stats);
end

function print_summary(stats)
    fprintf('\n==== WIENER TAP PROBE ====\n');
    fprintf('tap  success  dbg  win  sigma2    Ryy(0)    Rhh(0)    w0       sum(w)\n');
    for i = 1:numel(stats.tap)
        fprintf('%3d  %7d  %3d  %3d  %8.4f  %8.4f  %8.4f  %8.4f  %8.4f\n', ...
            stats.tap(i), stats.success(i), stats.has_dbg(i), stats.window_used(i), ...
            stats.sigma2(i), stats.ryy0(i), stats.rhh0(i), stats.weight0(i), stats.weight_sum(i));
    end
end

function cfg = default_cfg()
    cfg = struct();
    cfg.base_params = default_base_params();
    cfg.base_params.modulation = 'QPSK';
    cfg.base_params.use_pilots = true;
    cfg.base_params.used_bins = [2 3 4 5 6 7 8 9];
    cfg.base_params.pilot_bins = [2 4 7 9];
    cfg.tap_list = (1:8).';
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
    p.payload_bytes = 64;
    p.session_id = uint16(1234);
    p.detect_threshold = 0.005;
    p.wake_min_score = 0.005;
    p.wake_search_len = 8000;
    p.wake_ref_pre_ms = 15;
    p.sync_search_len = 1500;
    p.cfo_grid_hz = -12:2:12;
    p.rx_preroll = 128;
    p.pll_enable = true;
    p.pll_kp = 0.05;
    p.pll_ki = 0.002;
    p.snr_db = 24;
    p.cfo_hz = 4;
    p.timing_offset = 120;
    p.channel_taps = [1.0; 0.22; -0.08];
    p.echo_delays_ms = [];
    p.echo_gains = [];
    p.echo_phases_deg = [];
    p.apply_am_ripple = true;
    p.am_ripple_depth = 0.03;
    p.am_ripple_hz = 40;
    p.disable_lpf = false;
    p.disable_modulation = false;
    p.oracle_wake_start = [];
    p.oracle_sync_start = [];
    p.oracle_cfo_est_hz = [];
    p.verbose = false;
end
