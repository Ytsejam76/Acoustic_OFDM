% RUN_BER_SNR_PLOT  Generate BER/PER vs SNR plot for the OFDM link.
%
% Run:
%   octave --quiet run_ber_snr_plot.m

cfg = struct();
cfg.snr_db_list = 0:2:30;
cfg.num_trials = 100;
cfg.compare_oracle = true;
cfg.save_plot = true;
cfg.plot_filename = 'ber_per_vs_snr.png';
cfg.out_dir = fullfile(pwd, 'output_plots');

stats = ofdm_snr_sweep(cfg); %#ok<NASGU>

