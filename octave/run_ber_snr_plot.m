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

args = argv();
mods = {'BPSK', 'QPSK'};
echo_profile = 'none';
output_name = 'ber_per_snr.png';
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
    elseif strcmp(args{i}, '--echo')
        echo_profile = lower(strtrim(args{i+1}));
    elseif strcmp(args{i}, '--output')
        output_name = strtrim(args{i+1});
        if isempty(output_name)
            error('--output requires a non-empty filename');
        end
        if isempty(regexp(lower(output_name), '\.png$', 'once'))
            output_name = [output_name '.png'];
        end
    end
end

out_dir = fullfile(pwd, '..', 'images');
stats_paths = struct();
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
        cfg.num_trials = 50;
        cfg.compare_oracle = true;
        cfg.save_plot = true;
        cfg.plot_filename = sprintf('octave_tmp_%s_%s.png', mod_tag, echo_tag);
        cfg.out_dir = out_dir;
        cfg.base_params = struct();
        cfg.base_params.modulation = mod_name;
        cfg.base_params.use_pilots = [];
        cfg.base_params.echo_profile = echo_name;

        stats = ofdm_snr_sweep(cfg); %#ok<NASGU>
        stats_path = fullfile(cfg.out_dir, sprintf('snr_sweep_stats_%s_%s.mat', mod_tag, echo_tag));
        save('-v7', stats_path, 'stats');
        stats_paths.(echo_tag).(mod_tag) = stats_path;

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

        octave_tmp = fullfile(cfg.out_dir, cfg.plot_filename);
        if exist(octave_tmp, 'file') == 2
            delete(octave_tmp);
        end
    end
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
