% Copyright (c) 2026 Elias S. G. Carotti

function [noise_eq_used, state, dbg] = ofdm_estimate_disturbance_psd(Yfull, Hest, used_bins, p, state)
% OFDM_ESTIMATE_DISTURBANCE_PSD  Estimate post-EQ disturbance power on used bins.
%
% This estimator treats all unused FFT bins as measurements of residual
% non-signal energy. The resulting spectrum is interpolated onto the active
% bins and converted into the equalized domain by dividing by |Hest|^2.

    if nargin < 5 || isempty(state)
        state = ofdm_equalizer_init_state();
    end

    dbg = struct();
    nbins = numel(used_bins);
    noise_eq_used = zeros(nbins, 1);
    if isempty(Yfull) || isempty(Hest) || isempty(used_bins)
        return;
    end

    all_bins = fft_observation_bins(numel(Yfull));
    noise_bins = setdiff(all_bins, used_bins(:).', 'stable');
    if isempty(noise_bins)
        return;
    end

    noise_rx = abs(Yfull(noise_bins)).^2 / max(1, numel(Yfull));
    noise_rx = smooth_frequency_samples(noise_rx, p);
    interp_rx = interpolate_noise_psd(noise_bins, noise_rx, used_bins(:));

    if isfield(p, 'disturbance_psd_gain') && ~isempty(p.disturbance_psd_gain)
        interp_rx = double(p.disturbance_psd_gain) * interp_rx;
    end

    temporal_alpha = selected_disturbance_temporal_alpha(p);
    if ~isfield(state, 'disturbance_psd_est') || isempty(state.disturbance_psd_est)
        state.disturbance_psd_est = interp_rx(:);
    else
        state.disturbance_psd_est = temporal_alpha * state.disturbance_psd_est(:) ...
            + (1 - temporal_alpha) * interp_rx(:);
    end

    hpow = max(1.0e-6, abs(Hest(:)).^2);
    noise_eq_used = state.disturbance_psd_est(:) ./ hpow;

    dbg.noise_bins = noise_bins(:);
    dbg.noise_rx_raw = abs(Yfull(noise_bins)).^2 / max(1, numel(Yfull));
    dbg.noise_rx_interp_used = interp_rx(:);
    dbg.noise_eq_used = noise_eq_used(:);
    dbg.temporal_alpha = temporal_alpha;
end

function bins = fft_observation_bins(Nfft)
    bins = (2:Nfft).';
end

function vals = smooth_frequency_samples(vals, p)
    width = 3;
    if isfield(p, 'disturbance_freq_smooth') && ~isempty(p.disturbance_freq_smooth)
        width = max(1, round(double(p.disturbance_freq_smooth)));
    end
    if width <= 1 || numel(vals) <= 2
        vals = vals(:);
        return;
    end

    radius = floor(width / 2);
    tmp = vals(:);
    out = tmp;
    for i = 1:numel(tmp)
        lo = max(1, i - radius);
        hi = min(numel(tmp), i + radius);
        out(i) = mean(tmp(lo:hi));
    end
    vals = out;
end

function interp_vals = interpolate_noise_psd(src_bins, src_vals, dst_bins)
    if isempty(src_bins) || isempty(src_vals)
        interp_vals = zeros(numel(dst_bins), 1);
        return;
    end
    if numel(src_bins) == 1
        interp_vals = repmat(src_vals(1), numel(dst_bins), 1);
        return;
    end

    interp_vals = interp1(double(src_bins(:)), double(src_vals(:)), ...
        double(dst_bins(:)), 'linear', 'extrap');
    interp_vals = max(0, interp_vals(:));
end

function alpha = selected_disturbance_temporal_alpha(p)
    alpha = 0.75;
    if isfield(p, 'disturbance_temporal_alpha') && ~isempty(p.disturbance_temporal_alpha)
        alpha = min(0.98, max(0.0, double(p.disturbance_temporal_alpha)));
    end
end
