% Copyright (c) 2026 Elias S. G. Carotti

% RUN_BER_SNR_PLOT  Generate BER/PER vs SNR plots for OFDM.
%
% Run:
%   octave --quiet run_ber_snr_plot.m
%   octave --quiet run_ber_snr_plot.m --mod BPSK
%   octave --quiet run_ber_snr_plot.m --mod QPSK
%   octave --quiet run_ber_snr_plot.m --mod both
%   octave --quiet run_ber_snr_plot.m --mod both --echo cp_mix
%   octave --quiet run_ber_snr_plot.m --output my_plot.png
%   octave --quiet run_ber_snr_plot.m --with-oracle
%   octave --quiet run_ber_snr_plot.m --help

args = argv();
mods = {'BPSK', 'QPSK'};
echo_profile = 'none';
output_name = 'ber_per_snr.png';
with_oracle = false;
i = 1;
while i <= numel(args)
    if strcmp(args{i}, '--help') || strcmp(args{i}, '-h')
        fprintf(['Usage: octave --quiet run_ber_snr_plot.m [options]\n\n' ...
            'Options:\n' ...
            '  --help, -h            Show this help and exit.\n' ...
            '  --mod MOD             Select modulation: BPSK | QPSK | both (default: both).\n' ...
            '  --echo PROFILE        Echo profile: none | room_mild | cp_mix | all (default: none).\n' ...
            '  --output FILE.png     Output BER/PER plot filename (default: ber_per_snr.png).\n' ...
            '  --with-oracle         Include oracle-sync curves (default: disabled).\n\n' ...
            'Examples:\n' ...
            '  octave --quiet run_ber_snr_plot.m --mod both --echo all\n' ...
            '  octave --quiet run_ber_snr_plot.m --mod both --with-oracle --output my_plot.png\n']);
        return;
    elseif strcmp(args{i}, '--mod')
        if i + 1 > numel(args)
            error('--mod requires a value');
        end
        v = upper(strtrim(args{i+1}));
        switch v
            case 'BPSK'
                mods = {'BPSK'};
            case 'QPSK'
                mods = {'QPSK'};
            case 'BOTH'
                mods = {'BPSK', 'QPSK'};
            otherwise
                error('Unsupported --mod value: %s', args{i+1});
        end
        i = i + 2;
    elseif strcmp(args{i}, '--echo')
        if i + 1 > numel(args)
            error('--echo requires a value');
        end
        echo_profile = lower(strtrim(args{i+1}));
        i = i + 2;
    elseif strcmp(args{i}, '--output')
        if i + 1 > numel(args)
            error('--output requires a value');
        end
        output_name = strtrim(args{i+1});
        if isempty(output_name)
            error('--output requires a non-empty filename');
        end
        if isempty(regexp(lower(output_name), '\.png$', 'once'))
            output_name = [output_name '.png'];
        end
        i = i + 2;
    elseif strcmp(args{i}, '--with-oracle')
        with_oracle = true;
        i = i + 1;
    else
        error('Unknown option: %s', args{i});
    end
endwhile

out_dir = fullfile(pwd, '..', 'images');
data_dir = fullfile(pwd, '..', 'output_data');
if exist(data_dir, 'dir') ~= 7
    mkdir(data_dir);
end
stats_paths = struct();
const_dbg = struct();
time_dbg = struct();
mod_cfg = struct();
echo_profiles = {echo_profile};
if strcmpi(echo_profile, 'all')
    echo_profiles = {'none', 'room_mild', 'cp_mix'};
end

for ei = 1:numel(echo_profiles)
    echo_name = echo_profiles{ei};
    for mi = 1:numel(mods)
        mod_name = mods{mi};
        mod_tag = lower(mod_name);
        echo_tag = lower(echo_name);

        cfg = struct();
        cfg.snr_db_list = 0:1:30;
        cfg.num_trials = 20;
        cfg.compare_oracle = with_oracle;
        cfg.make_octave_plots = false;
        cfg.save_plot = false;
        cfg.plot_filename = sprintf('octave_tmp_%s_%s.png', mod_tag, echo_tag);
        cfg.out_dir = out_dir;
        cfg.base_params = struct();
        cfg.base_params.modulation = mod_name;
        cfg.base_params.use_pilots = [];
        cfg.base_params.echo_profile = echo_name;

        stats = ofdm_snr_sweep(cfg); %#ok<NASGU>
        stats_path = fullfile(data_dir, sprintf('snr_sweep_stats_%s_%s.mat', mod_tag, echo_tag));
        save('-v7', stats_path, 'stats');
        stats_paths.(echo_tag).(mod_tag) = stats_path;

        % Save one constellation snapshot (before/after EQ) for this mode/profile.
        ps = struct();
        ps.modulation = mod_name;
        ps.use_pilots = [];
        ps.pause_before_exit = false;
        ps.make_plots = false;
        ps.save_images = false;
        ps.save_decoder_constellation = false;
        ps.out_dir = out_dir;
        ps.verbose = false;
        ps.snr_db = 30;
        ps.cfo_hz = 0;
        ps.timing_offset = 0;
        ps.channel_taps = 1;
        ps.apply_am_ripple = false;
        ps.echo_delays_ms = [];
        ps.echo_gains = [];
        ps.echo_phases_deg = [];
        if strcmpi(echo_name, 'room_mild')
            cp_ms = 1e3 * 72 / 48000;
            ps.echo_delays_ms = [0.20, 0.55, 0.90, 1.35] * cp_ms;
            ps.echo_gains = [0.25, 0.12, 0.06, 0.03];
            ps.echo_phases_deg = [0, 20, -35, 50];
        elseif strcmpi(echo_name, 'cp_mix')
            cp_ms = 1e3 * 72 / 48000;
            inside = [0.30, 0.70, 0.93] * cp_ms;
            outside = [1.35, 2.20] * cp_ms;
            ps.echo_delays_ms = [inside, outside];
            ps.echo_gains = [0.40, 0.28, 0.18, 0.12, 0.08];
            ps.echo_phases_deg = [0, 25, -40, 60, -90];
        end
        [~, dbg_const, tx_const, rx_const, ~, p_used] = ofdm_test_channel(ps);
        mod_cfg.(mod_tag) = p_used;
        if isempty(ps.echo_delays_ms)
            echo_desc = 'No echoes';
        else
            max_d_ms = max(ps.echo_delays_ms);
            max_d_m = 0.3 * max_d_ms;
            echo_desc = sprintf('Echoes: N=%d, max %.2f ms (~%.2f m)', ...
                numel(ps.echo_delays_ms), max_d_ms, max_d_m);
        end
        cdbg = struct();
        cdbg.dbg = dbg_const;
        cdbg.snr_db = ps.snr_db;
        cdbg.cfo_hz = ps.cfo_hz;
        cdbg.echo_desc = echo_desc;
        const_dbg.(mod_tag).(echo_tag) = cdbg;

        if ~isempty(tx_const) && ~isempty(rx_const)
            td = struct();
            td.tx = tx_const(:);
            td.rx = rx_const(:);
            td.fs = p_used.fs;
            td.snr_db = ps.snr_db;
            td.cfo_hz = ps.cfo_hz;
            td.echo_desc = echo_desc;
            time_dbg.(mod_tag).(echo_tag) = td;
        end
        close all;

        if numel(mods) == 1 && numel(echo_profiles) == 1
            png_smooth = fullfile(cfg.out_dir, output_name);

            cmd = sprintf('./.venv/bin/python3 plot_snr_sweep_seaborn.py --stats \"%s\" --out \"%s\"', ...
                stats_path, png_smooth);
            [status, out] = system(cmd);
            fprintf('%s', out);
            if status ~= 0
                error('Seaborn rendering failed');
            end
        end

    end
end

% Consolidated constellation image.
if numel(mods) == 2
    if numel(echo_profiles) == 3
        row_echoes = {'none', 'room_mild', 'cp_mix'};
    else
        row_echoes = echo_profiles;
    end
    hconst = figure('name', 'constellation_compare_all', 'visible', 'off', ...
        'position', [100 100 1200 1400]);
    set(hconst, 'paperpositionmode', 'auto');
    nrows = numel(row_echoes);
    ncols = 2;
    lm = 0.055; rm = 0.020; tm = 0.070; bm = 0.045;
    gx = 0.040; gy = 0.050;
    axw = (1 - lm - rm - (ncols - 1) * gx) / ncols;
    axh = (1 - tm - bm - (nrows - 1) * gy) / nrows;
    for ri = 1:numel(row_echoes)
        echo_tag = lower(row_echoes{ri});
        for ci = 1:2
            mod_name = mods{ci};
            mod_tag = lower(mod_name);
            left = lm + (ci - 1) * (axw + gx);
            bottom = 1 - tm - ri * axh - (ri - 1) * gy;
            ax = axes('parent', hconst, 'position', [left bottom axw axh]); %#ok<NASGU>
            hold on;
            if isfield(const_dbg, mod_tag) && isfield(const_dbg.(mod_tag), echo_tag)
                cdbg = const_dbg.(mod_tag).(echo_tag);
                dbg = cdbg.dbg;
                has_raw = isfield(dbg, 'rx_syms_raw') && ~isempty(dbg.rx_syms_raw);
                has_eq = isfield(dbg, 'rx_syms_eq') && ~isempty(dbg.rx_syms_eq);
                if has_raw
                    plot(real(dbg.rx_syms_raw), imag(dbg.rx_syms_raw), '.', ...
                        'color', [0.35 0.35 0.35], 'markersize', 7);
                end
                if has_eq
                    plot(real(dbg.rx_syms_eq), imag(dbg.rx_syms_eq), 'x', ...
                        'color', [0.10 0.45 0.85], 'markersize', 5, 'linewidth', 1.0);
                end
                if strcmpi(mod_name, 'BPSK')
                    ideal = [-1; 1];
                else
                    ideal = [1+1j; 1-1j; -1+1j; -1-1j] / sqrt(2);
                end
                plot(real(ideal), imag(ideal), 'ro', 'markersize', 8, 'linewidth', 1.2);
                lgd = {};
                if has_raw
                    lgd{end+1} = 'Pre-EQ'; %#ok<AGROW>
                end
                if has_eq
                    lgd{end+1} = 'Post-EQ'; %#ok<AGROW>
                end
                lgd{end+1} = 'Ideal'; %#ok<AGROW>
                legend(lgd, 'location', 'northeast');
            end
            grid on;
            axis equal;
            xlabel('In-Phase');
            ylabel('Quadrature');
            if isfield(const_dbg, mod_tag) && isfield(const_dbg.(mod_tag), echo_tag)
                ttl = sprintf('%s | %s | SNR: %.1f dB', upper(mod_name), cdbg.echo_desc, cdbg.snr_db);
            else
                ttl = sprintf('%s | Echo profile: %s', upper(mod_name), echo_tag);
            end
            title(ttl, 'interpreter', 'none');
            hold off;
        end
    end
    pilots_bpsk = NaN;
    pilots_qpsk = NaN;
    bw_khz = NaN;
    if isfield(mod_cfg, 'bpsk')
        cb = mod_cfg.bpsk;
        used = cb.used_bins(:).';
        pilots = [];
        if isfield(cb, 'use_pilots') && ~isempty(cb.use_pilots)
            pon = logical(cb.use_pilots);
        else
            pon = strcmpi(cb.modulation, 'QPSK');
        end
        if pon
            pilots = intersect(used, cb.pilot_bins(:).', 'stable');
            if isfield(cb, 'num_pilots') && ~isempty(cb.num_pilots)
                np = max(0, min(round(double(cb.num_pilots(1))), numel(pilots)));
                pilots = pilots(1:np);
            end
        end
        pilots_bpsk = numel(pilots);
        if ~isempty(used)
            bw_khz = ((max(used) - min(used) + 1) * cb.fs / cb.Nfft) / 1e3;
        end
    end
    if isfield(mod_cfg, 'qpsk')
        cq = mod_cfg.qpsk;
        used = cq.used_bins(:).';
        pilots = [];
        if isfield(cq, 'use_pilots') && ~isempty(cq.use_pilots)
            pon = logical(cq.use_pilots);
        else
            pon = strcmpi(cq.modulation, 'QPSK');
        end
        if pon
            pilots = intersect(used, cq.pilot_bins(:).', 'stable');
            if isfield(cq, 'num_pilots') && ~isempty(cq.num_pilots)
                np = max(0, min(round(double(cq.num_pilots(1))), numel(pilots)));
                pilots = pilots(1:np);
            end
        end
        pilots_qpsk = numel(pilots);
        if isnan(bw_khz) && ~isempty(used)
            bw_khz = ((max(used) - min(used) + 1) * cq.fs / cq.Nfft) / 1e3;
        end
    end
    if isnan(pilots_bpsk)
        pb_str = 'n/a';
    else
        pb_str = sprintf('%d', round(pilots_bpsk));
    end
    if isnan(pilots_qpsk)
        pq_str = 'n/a';
    else
        pq_str = sprintf('%d', round(pilots_qpsk));
    end
    if isnan(bw_khz)
        bw_str = 'n/a';
    else
        bw_str = sprintf('%.2f kHz', bw_khz);
    end
    if exist('sgtitle', 'builtin') || exist('sgtitle', 'file')
        sgtitle({
            'Constellation Comparison (Pre/Post Equalization)'
            sprintf('Channels: %d | Pilots: BPSK=%s, QPSK=%s | BW: %s', ...
                numel(row_echoes), pb_str, pq_str, bw_str)
        }, 'interpreter', 'none');
    end
    print(hconst, fullfile(out_dir, 'constellation_compare.png'), '-dpng', '-r180');
    close(hconst);
end

% Consolidated time-domain images: one figure per modulation.
for mi = 1:numel(mods)
    mod_name = mods{mi};
    mod_tag = lower(mod_name);
    if ~isfield(time_dbg, mod_tag)
        continue;
    end
    if numel(echo_profiles) == 3
        row_echoes = {'none', 'room_mild', 'cp_mix'};
    else
        row_echoes = echo_profiles;
    end
    htd = figure('name', sprintf('time_domain_compare_%s', mod_tag), 'visible', 'off', ...
        'position', [120 120 1200 1300]);
    set(htd, 'paperpositionmode', 'auto');
    nrows = numel(row_echoes);
    lm = 0.070; rm = 0.020; tm = 0.075; bm = 0.055; gy = 0.055;
    axh = (1 - tm - bm - (nrows - 1) * gy) / nrows;
    for ri = 1:nrows
        echo_tag = lower(row_echoes{ri});
        if ~isfield(time_dbg.(mod_tag), echo_tag)
            continue;
        end
        td = time_dbg.(mod_tag).(echo_tag);
        bottom = 1 - tm - ri * axh - (ri - 1) * gy;
        axes('parent', htd, 'position', [lm bottom 1-lm-rm axh]); %#ok<LAXES>
        hold on;
        nt = (0:numel(td.tx)-1).';
        nr = (0:numel(td.rx)-1).';
        plot(nt / td.fs * 1e3, td.tx, '-', 'color', [0.20 0.20 0.20], 'linewidth', 0.9);
        plot(nr / td.fs * 1e3, td.rx, '-', 'color', [0.10 0.45 0.85], 'linewidth', 0.9);
        grid on;
        xlabel('Time (ms)');
        ylabel('Amplitude');
        title(sprintf('%s | SNR: %.1f dB | CFO: %.1f Hz', td.echo_desc, td.snr_db, td.cfo_hz), ...
            'interpreter', 'none');
        legend({'TX (original)', 'RX (after channel)'}, 'location', 'northeast');
        hold off;
    end
    if exist('sgtitle', 'builtin') || exist('sgtitle', 'file')
        sgtitle({
            sprintf('%s Time-Domain Waveforms (TX vs RX)', upper(mod_name))
            'Rows: channel models'
        }, 'interpreter', 'none');
    end
    out_time = fullfile(out_dir, sprintf('time_domain_compare_%s.png', mod_tag));
    print(htd, out_time, '-dpng', '-r180');
    close(htd);
end

if numel(mods) == 2
    if numel(echo_profiles) == 1
        ep = lower(echo_profiles{1});
        if isfield(stats_paths, ep) && isfield(stats_paths.(ep), 'bpsk') && isfield(stats_paths.(ep), 'qpsk')
            combo_smooth = fullfile(out_dir, output_name);

            cmd = sprintf(['./.venv/bin/python3 plot_snr_sweep_seaborn.py --stats-bpsk \"%s\" ', ...
                '--stats-qpsk \"%s\" --out \"%s\"'], ...
                stats_paths.(ep).bpsk, stats_paths.(ep).qpsk, combo_smooth);
            [status, out] = system(cmd);
            fprintf('%s', out);
            if status ~= 0
                error('Seaborn combined rendering failed');
            end
        end
    else
        labels = '';
        bpsk_list = '';
        qpsk_list = '';
        for ei = 1:numel(echo_profiles)
            ep = lower(echo_profiles{ei});
            if ei > 1
                labels = [labels ',']; %#ok<AGROW>
                bpsk_list = [bpsk_list ',']; %#ok<AGROW>
                qpsk_list = [qpsk_list ',']; %#ok<AGROW>
            end
            labels = [labels ep]; %#ok<AGROW>
            bpsk_list = [bpsk_list stats_paths.(ep).bpsk]; %#ok<AGROW>
            qpsk_list = [qpsk_list stats_paths.(ep).qpsk]; %#ok<AGROW>
        end

        combo_smooth = fullfile(out_dir, output_name);

        cmd = sprintf(['./.venv/bin/python3 plot_snr_sweep_seaborn.py --stats-bpsk-list \"%s\" ', ...
            '--stats-qpsk-list \"%s\" --labels \"%s\" --out \"%s\"'], ...
            bpsk_list, qpsk_list, labels, combo_smooth);
        [status, out] = system(cmd);
        fprintf('%s', out);
        if status ~= 0
            error('Seaborn combined rendering failed');
        end
    end
end
