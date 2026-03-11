% RUN_BER_SNR_PLOT  Generate BER/PER vs SNR plot for the OFDM link.
%
% Run:
%   octave --quiet run_ber_snr_plot.m

cfg = struct();
cfg.snr_db_list = 0:1:30;
cfg.num_trials = 20;
cfg.compare_oracle = true;
cfg.save_plot = true;
cfg.plot_filename = 'octave_tmp.png';
cfg.out_dir = fullfile(pwd, '..', 'images');

stats = ofdm_snr_sweep(cfg); %#ok<NASGU>

stats_path = fullfile(cfg.out_dir, 'snr_sweep_stats.mat');
png_path = fullfile(cfg.out_dir, 'ber_per_vs_snr.png');
octave_tmp = fullfile(cfg.out_dir, cfg.plot_filename);
cmd = sprintf('./.venv/bin/python3 plot_snr_sweep_seaborn.py --stats \"%s\" --out \"%s\"', ...
    stats_path, png_path);
[status, out] = system(cmd);
fprintf('%s', out);
if status ~= 0
    error('Seaborn rendering failed');
end
if exist(octave_tmp, 'file') == 2
    delete(octave_tmp);
end
