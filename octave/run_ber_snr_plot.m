% Copyright (c) 2026 Elias S. G. Carotti

% RUN_BER_SNR_PLOT  Generate BER/PER vs SNR plots for OFDM.
%
% Run:
%   octave --quiet run_ber_snr_plot.m
%   octave --quiet run_ber_snr_plot.m --mod BPSK
%   octave --quiet run_ber_snr_plot.m --mod QPSK
%   octave --quiet run_ber_snr_plot.m --mod both

args = argv();
mods = {'BPSK', 'QPSK'};
for i = 1:max(0, numel(args)-1)
    if strcmp(args{i}, '--mod')
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
        break;
    end
end

out_dir = fullfile(pwd, '..', 'images');

for mi = 1:numel(mods)
    mod_name = mods{mi};
    mod_tag = lower(mod_name);

    cfg = struct();
    cfg.snr_db_list = 0:1:30;
    cfg.num_trials = 20;
    cfg.compare_oracle = true;
    cfg.save_plot = true;
    cfg.plot_filename = sprintf('octave_tmp_%s.png', mod_tag);
    cfg.out_dir = out_dir;
    cfg.base_params = struct();
    cfg.base_params.modulation = mod_name;
    cfg.base_params.use_pilots = [];

    stats = ofdm_snr_sweep(cfg); %#ok<NASGU>
    stats_path = fullfile(cfg.out_dir, sprintf('snr_sweep_stats_%s.mat', mod_tag));
    save('-v7', stats_path, 'stats');

    png_smooth = fullfile(cfg.out_dir, sprintf('ber_per_vs_snr_%s.png', mod_tag));
    png_raw = fullfile(cfg.out_dir, sprintf('ber_per_vs_snr_%s_raw.png', mod_tag));

    cmd = sprintf('./.venv/bin/python3 plot_snr_sweep_seaborn.py --stats \"%s\" --out \"%s\"', ...
        stats_path, png_smooth);
    [status, out] = system(cmd);
    fprintf('%s', out);
    if status ~= 0
        error('Seaborn rendering failed');
    end

    cmd = sprintf('./.venv/bin/python3 plot_snr_sweep_seaborn.py --stats \"%s\" --out \"%s\" --no-smooth', ...
        stats_path, png_raw);
    [status, out] = system(cmd);
    fprintf('%s', out);
    if status ~= 0
        error('Seaborn rendering failed');
    end

    if strcmpi(mod_name, 'BPSK')
        copyfile(png_smooth, fullfile(cfg.out_dir, 'ber_per_vs_snr.png'));
        copyfile(png_raw, fullfile(cfg.out_dir, 'ber_per_vs_snr_raw.png'));
    end

    octave_tmp = fullfile(cfg.out_dir, cfg.plot_filename);
    if exist(octave_tmp, 'file') == 2
        delete(octave_tmp);
    end
end
