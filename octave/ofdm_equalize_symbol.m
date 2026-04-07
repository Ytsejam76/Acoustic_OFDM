% Copyright (c) 2026 Elias S. G. Carotti

function [Xeq_used, state, dbg] = ofdm_equalize_symbol(Xraw_used, Hest, used_bins, pilot_bins, pref, p, state)
% OFDM_EQUALIZE_SYMBOL  Equalize one OFDM symbol with configurable pilot processing.
%
%   [Xeq_used, state, dbg] = ofdm_equalize_symbol(Xraw_used, Hest, used_bins, pilot_bins, pref, p, state)

    if nargin < 7 || isempty(state)
        state = ofdm_equalizer_init_state();
    end

    dbg = struct();
    Xeq_used = Xraw_used ./ Hest;
    [has_pilot, pilot_pos, ~] = bin_positions_local(used_bins, pilot_bins, setdiff(used_bins, pilot_bins, 'stable'));
    if ~has_pilot || isempty(pref)
        return;
    end

    mode = equalizer_mode_local(p);
    residual_modes = {'pilot-denoise', 'pilot-denoise-mdl', 'pilot-denoise-temporal'};
    if any(strcmpi(mode, residual_modes))
        [residual_curve, residual_var, state, dbg] = residual_curve_for_mode(Xeq_used, used_bins, pilot_pos, pref, p, state, mode);
        if ~isempty(residual_curve)
            valid = abs(residual_curve) > 1e-3 & isfinite(real(residual_curve)) & isfinite(imag(residual_curve));
            Xeq_used(valid) = Xeq_used(valid) ./ residual_curve(valid);
        end
        dbg.residual_var = residual_var;
    end

    [g, pilot_residual] = common_pilot_gain(Xeq_used, pilot_pos, pref);
    dbg.pilot_gain = g;
    dbg.pilot_residual = pilot_residual;
    if isfinite(g) && abs(g) > 1e-12
        Xeq_used = Xeq_used / g;
    end
end

function state = ofdm_equalizer_init_state()
    state = struct();
    state.residual_curve_est = [];
end

function mode = equalizer_mode_local(p)
    mode = 'training-pilot';
    if isfield(p, 'equalizer_mode') && ~isempty(p.equalizer_mode)
        mode = char(p.equalizer_mode);
    end
end

function [residual_curve, residual_var, state, dbg] = residual_curve_for_mode(Xeq_used, used_bins, pilot_pos, pref, p, state, mode)
    dbg = struct();
    residual_curve = [];
    residual_var = [];

    pts = pilot_residual_points(Xeq_used, pilot_pos, pref);
    if rows(pts) < 2
        return;
    end

    residual_var = pilot_residual_noise_var(pts);
    curve_len = numel(used_bins);
    order_mode = 'all';
    fixed_order = 1;
    max_order = min(rows(pts), curve_len);
    switch lower(strtrim(mode))
        case 'pilot-denoise-mdl'
            order_mode = 'mdl';
            if isfield(p, 'residual_tap_order_max') && ~isempty(p.residual_tap_order_max)
                max_order = min(max_order, max(1, round(double(p.residual_tap_order_max))));
            end
        case 'pilot-denoise'
            if isfield(p, 'residual_tap_order_mode') && strcmpi(p.residual_tap_order_mode, 'fixed')
                order_mode = 'fixed';
                fixed_order = max(1, round(double(p.residual_tap_order)));
            end
        case 'pilot-denoise-temporal'
            if isfield(p, 'residual_tap_order_mode') && strcmpi(p.residual_tap_order_mode, 'fixed')
                order_mode = 'fixed';
                fixed_order = max(1, round(double(p.residual_tap_order)));
            end
    end

    selected_order = [];
    switch order_mode
        case 'all'
            curve = interpolate_residual_curve(pts, curve_len);
            taps = unitary_ifft(curve);
        case 'fixed'
            selected_order = min([fixed_order, rows(pts), curve_len]);
            taps = fit_delay_taps_least_squares(pts, curve_len, selected_order);
        case 'mdl'
            selected_order = select_delay_tap_order_mdl(pts, curve_len, max_order);
            taps = fit_delay_taps_least_squares(pts, curve_len, selected_order);
        otherwise
            error('Unsupported residual tap order mode');
    end

    if isempty(taps)
        return;
    end

    if ~isempty(selected_order) && selected_order < numel(taps)
        taps(selected_order+1:end) = 0;
    end

    for i = 1:numel(taps)
        power = abs(taps(i))^2;
        shrink = max(0, min(1, (power - residual_var) / (power + 1e-9)));
        taps(i) = taps(i) * shrink;
    end

    residual_curve = unitary_fft(taps);
    residual_curve = clamp_curve_magnitude(residual_curve);
    dbg.selected_order = selected_order;

    if strcmpi(mode, 'pilot-denoise-temporal')
        [residual_curve, state.residual_curve_est] = temporal_ema_residual_curve(residual_curve, state.residual_curve_est, p);
    end
end

function pts = pilot_residual_points(Xeq_used, pilot_pos, pref)
    pts = zeros(numel(pilot_pos), 2);
    n = 0;
    for k = 1:numel(pilot_pos)
        ref_sym = pref(k);
        if abs(ref_sym)^2 <= 1e-9
            continue;
        end
        n = n + 1;
        pts(n, 1) = pilot_pos(k);
        pts(n, 2) = Xeq_used(pilot_pos(k)) / ref_sym;
    end
    pts = pts(1:n, :);
end

function curve = interpolate_residual_curve(pts, curve_len)
    curve = ones(curve_len, 1);
    for i = 1:curve_len
        if i <= pts(1,1)
            curve(i) = pts(1,2);
            continue;
        end
        if i >= pts(end,1)
            curve(i) = pts(end,2);
            continue;
        end
        seg = find(pts(:,1) >= i, 1) - 1;
        if isempty(seg) || seg < 1
            seg = 1;
        end
        i0 = pts(seg,1);
        i1 = pts(seg+1,1);
        z0 = pts(seg,2);
        z1 = pts(seg+1,2);
        t = max(0, min(1, (i - i0) / max(1, i1 - i0)));
        curve(i) = z0 * (1 - t) + z1 * t;
    end
end

function noise_var = pilot_residual_noise_var(pts)
    if rows(pts) < 3
        noise_var = 5.0e-3;
        return;
    end

    err_sum = 0;
    count = 0;
    for i = 1:(rows(pts) - 2)
        i0 = pts(i,1);
        z0 = pts(i,2);
        i1 = pts(i+1,1);
        z1 = pts(i+1,2);
        i2 = pts(i+2,1);
        z2 = pts(i+2,2);
        t = max(0, min(1, (i1 - i0) / max(1, i2 - i0)));
        pred = z0 * (1 - t) + z2 * t;
        err_sum = err_sum + abs(z1 - pred)^2;
        count = count + 1;
    end

    if count == 0
        noise_var = 5.0e-3;
    else
        noise_var = min(2.5e-1, max(1.0e-4, err_sum / count));
    end
end

function taps = fit_delay_taps_least_squares(pts, curve_len, order)
    order = min([max(1, order), rows(pts), curve_len]);
    gram = zeros(order, order);
    rhs = zeros(order, 1);
    for r = 1:rows(pts)
        basis = delay_basis_row(curve_len, pts(r,1), order);
        obs = pts(r,2);
        rhs = rhs + conj(basis).' * obs;
        gram = gram + conj(basis).' * basis;
    end
    gram = gram + 1e-6 * eye(order);
    taps = gram \ rhs;
    taps = [taps; zeros(curve_len - order, 1)];
end

function order = select_delay_tap_order_mdl(pts, curve_len, max_order)
    sample_count = max(1, 2 * rows(pts));
    max_order = min([max(1, max_order), rows(pts), curve_len]);
    best_score = Inf;
    order = 1;
    for k = 1:max_order
        taps = fit_delay_taps_least_squares(pts, curve_len, k);
        mse = max(1.0e-6, min(1.0e3, delay_tap_fit_mse(pts, curve_len, taps)));
        n_params = 2 * k;
        score = sample_count * log(mse) + n_params * log(sample_count);
        if score < best_score
            best_score = score;
            order = k;
        end
    end
end

function mse = delay_tap_fit_mse(pts, curve_len, taps)
    err = 0;
    for i = 1:rows(pts)
        est = synthesize_delay_response_at(curve_len, taps, pts(i,1));
        err = err + abs(pts(i,2) - est)^2;
    end
    mse = err / rows(pts);
end

function basis = delay_basis_row(curve_len, pos, order)
    tap_idx = 0:(order-1);
    basis = exp(-1j * 2*pi * (pos - 1) * tap_idx / curve_len) / sqrt(curve_len);
end

function value = synthesize_delay_response_at(curve_len, taps, pos)
    tap_idx = (0:numel(taps)-1).';
    basis = exp(-1j * 2*pi * (pos - 1) * tap_idx / curve_len) / sqrt(curve_len);
    value = sum(taps(:) .* basis);
end

function y = unitary_fft(x)
    y = fft(x) / sqrt(numel(x));
end

function y = unitary_ifft(x)
    y = ifft(x) * sqrt(numel(x));
end

function curve = clamp_curve_magnitude(curve)
    mags = abs(curve);
    phases = angle(curve);
    mags = min(1.5, max(0.75, mags));
    curve = mags .* exp(1j * phases);
end

function [fused, state_curve] = temporal_ema_residual_curve(current, prev, p)
    if nargin < 3
        p = struct();
    end
    keep = 4;
    if isfield(p, 'temporal_window') && ~isempty(p.temporal_window)
        keep = max(1, round(double(p.temporal_window)));
    end
    if keep <= 1 || isempty(prev) || numel(prev) ~= numel(current)
        fused = current;
        state_curve = current;
        return;
    end

    alpha = min(1.0, max(0.05, 1.0 / keep));
    phase_align = common_phase_delta(prev, current);
    rot = exp(-1j * phase_align);
    current_aligned = current * rot;
    fused = alpha * current_aligned + (1 - alpha) * prev;
    state_curve = fused;
end

function delta = common_phase_delta(reference, current)
    acc = 0;
    for i = 1:min(numel(reference), numel(current))
        rn = abs(reference(i));
        cn = abs(current(i));
        if rn > 1e-6 && cn > 1e-6
            acc = acc + (current(i) / cn) * conj(reference(i) / rn);
        end
    end
    if abs(acc) > 1e-9
        delta = angle(acc);
    else
        delta = 0;
    end
end

function [g, pilot_residual] = common_pilot_gain(Xeq_used, pilot_pos, pref)
    prx = Xeq_used(pilot_pos);
    g = sum(prx .* conj(pref)) / (sum(abs(pref).^2) + 1e-12);
    pilot_residual = prx;
    if ~isfinite(real(g)) || ~isfinite(imag(g))
        g = 1;
    end
end

function [has_pilot, pilot_pos, data_pos] = bin_positions_local(used_bins, pilot_bins, data_bins)
    if isempty(pilot_bins)
        has_pilot = false;
        pilot_pos = [];
    else
        [tfp, pilot_pos] = ismember(pilot_bins, used_bins);
        has_pilot = all(tfp);
        if ~has_pilot
            pilot_pos = [];
        end
    end
    [tfd, data_pos] = ismember(data_bins, used_bins);
    if ~all(tfd)
        data_pos = [];
    end
end
