// Copyright (c) 2026 Elias S. G. Carotti

use std::collections::VecDeque;
use std::f32::consts::PI;

use rustfft::num_complex::Complex32;
use rustfft::FftPlanner;

use crate::config::{EqualizerFeatures, Modulation, OfdmConfig, ResidualTapOrderMode};

/// Packet-local equalizer tracking state.
///
/// Rationale:
/// The equalizer keeps a short history of pilot-derived phase-line parameters
/// and pilot residual variance. This supports temporal least-squares tracking
/// and MMSE regularization without contaminating the static training baseline.
#[derive(Clone, Debug, Default)]
pub(crate) struct EqualizerTrackingState {
    phase_line_hist: VecDeque<(f32, f32, f32)>,
    residual_curve_est: Option<Vec<Complex32>>,
    next_symbol_time: f32,
    noise_var: Option<f32>,
}

/// Resets packet-local pilot tracking.
pub(crate) fn equalizer_reset_tracking(state: &mut EqualizerTrackingState) {
    *state = EqualizerTrackingState::default();
}

/// Estimates the training-symbol channel on the active carriers.
///
/// Rationale:
/// We use the known training symbol as a least-squares channel probe,
///
/// $$\hat H_{\text{train}}[k] = \frac{Y_{\text{train}}[k]}{X_{\text{train}}[k]}$$
///
/// which gives a coarse per-bin channel estimate before any pilot tracking.
pub(crate) fn estimate_channel_from_training(
    ytrain: &[Complex32],
    used_bins: &[usize],
    train_known: &[Complex32],
) -> Vec<Complex32> {
    let mut hest = vec![Complex32::new(1.0, 0.0); used_bins.len()];
    for (k, &bin) in used_bins.iter().enumerate() {
        hest[k] = ytrain[bin] / train_known[k];
    }
    hest
}

/// Chooses the baseline equalizer channel model.
///
/// Rationale:
/// When training-baseline equalization is enabled, the training estimate is
/// used as the static baseline. Otherwise we start from a flat unit channel,
///
/// $$\hat H_0[k] = 1$$
///
/// and let later pilot corrections carry the residual alignment.
pub(crate) fn equalizer_initial_channel(
    cfg: &OfdmConfig,
    ytrain: &[Complex32],
    used_bins: &[usize],
    train_known: &[Complex32],
) -> Vec<Complex32> {
    if cfg
        .equalizer
        .features
        .contains(EqualizerFeatures::TRAINING_BASELINE)
    {
        estimate_channel_from_training(ytrain, used_bins, train_known)
    } else {
        vec![Complex32::new(1.0, 0.0); used_bins.len()]
    }
}

/// Refreshes the baseline channel when a retraining symbol is present.
///
/// Rationale:
/// A retraining symbol resets the static channel prior,
///
/// $$\hat H[k] \leftarrow \hat H_{\text{train}}[k]$$
///
/// instead of forcing the pilot path to absorb long-term drift.
pub(crate) fn equalizer_refresh_channel(
    cfg: &OfdmConfig,
    hest: &mut Vec<Complex32>,
    y: &[Complex32],
    used_bins: &[usize],
    train_known: &[Complex32],
) {
    if cfg
        .equalizer
        .features
        .contains(EqualizerFeatures::TRAINING_BASELINE)
    {
        *hest = estimate_channel_from_training(y, used_bins, train_known);
    }
}

/// Applies a regularized single-bin equalizer.
///
/// Rationale:
/// A raw zero-forcing inverse,
///
/// $$\hat X[k] = Y[k] / \hat H[k]$$
///
/// blows up noise on weak bins. We instead use a bounded MMSE-like inverse,
///
/// $$G[k] = \frac{\hat H^*[k]}{|\hat H[k]|^2 + \varepsilon}$$
///
/// followed by an explicit gain clamp.
pub(crate) fn regularized_equalize(y: Complex32, h: Complex32) -> Complex32 {
    let h_pow = h.norm_sqr();
    let eps = (2.0e-2f32).max(1.0e-1f32 * h_pow);
    let mut g = h.conj() / (h_pow + eps);
    let gmax = 2.0f32;
    let gnorm = g.norm();
    if gnorm > gmax && gnorm.is_finite() {
        g *= gmax / gnorm;
    }
    y * g
}

/// Applies a noise-aware MMSE single-bin equalizer.
///
/// Rationale:
/// Once pilot residuals provide a noise-power estimate
/// $$
/// \sigma_n^2 \approx \mathbb{E}[|\hat X_{\text{pilot}}-X_{\text{pilot}}|^2],
/// $$
/// the equalizer can use the corresponding MMSE inverse
/// $$
/// G[k] = \frac{\hat H^*[k]}{|\hat H[k]|^2 + \sigma_n^2}
/// $$
/// rather than the legacy heuristic regularizer.
fn regularized_equalize_mmse(y: Complex32, h: Complex32, noise_var: f32) -> Complex32 {
    let h_pow = h.norm_sqr();
    let eps = noise_var.clamp(1.0e-4, 2.5e-1);
    let mut g = h.conj() / (h_pow + eps);
    let gmax = 2.0f32;
    let gnorm = g.norm();
    if gnorm > gmax && gnorm.is_finite() {
        g *= gmax / gnorm;
    }
    y * g
}

/// Equalizes one OFDM symbol against a baseline channel model.
///
/// Rationale:
/// We explicitly separate the training estimate into amplitude and phase,
///
/// $$\hat H_0[k] = \hat A_0[k] e^{j \hat \phi_0[k]}$$
///
/// and invert that baseline before applying any pilot-derived residual.
fn equalize_symbol_with_baseline(
    cfg: &OfdmConfig,
    state: &EqualizerTrackingState,
    y: &[Complex32],
    used_bins: &[usize],
    hest: &[Complex32],
) -> Vec<Complex32> {
    let (amp0, phase0) = baseline_channel_model(hest);
    let mut xeq_used = Vec::with_capacity(used_bins.len());
    for (k, &bin) in used_bins.iter().enumerate() {
        let h0 = Complex32::from_polar(amp0[k], phase0[k]);
        let z = if cfg
            .equalizer
            .features
            .contains(EqualizerFeatures::NOISE_AWARE_MMSE)
        {
            if let Some(noise_var) = state.noise_var {
                regularized_equalize_mmse(y[bin], h0, noise_var)
            } else {
                regularized_equalize(y[bin], h0)
            }
        } else {
            regularized_equalize(y[bin], h0)
        };
        xeq_used.push(z);
    }
    xeq_used
}

/// Equalizes one OFDM symbol with the current instantaneous pilot model.
///
/// Rationale:
/// The training symbol gives the static baseline amplitude/phase model, while
/// pilots capture the symbol-local residual phase trend. The current model fits
/// a linear phase residual over bin index,
///
/// $$\Delta \phi(k) = a + b k$$
///
/// then optionally fits a smooth residual amplitude line,
///
/// $$\Delta A(k) = c + d k$$
///
/// and applies both before a final common pilot normalization.
fn equalize_symbol_with_pilots_instantaneous(
    cfg: &OfdmConfig,
    state: &mut EqualizerTrackingState,
    y: &[Complex32],
    used_bins: &[usize],
    pilot_bins: &[usize],
    pref: &[Complex32],
    hest: &[Complex32],
) -> Vec<Complex32> {
    if used_bins.is_empty() || hest.len() != used_bins.len() {
        return Vec::new();
    }

    let mut xeq_used = equalize_symbol_with_baseline(cfg, state, y, used_bins, hest);

    if cfg
        .equalizer
        .features
        .contains(EqualizerFeatures::PILOT_IFFT_DENOISE)
        && !pilot_bins.is_empty()
        && !pref.is_empty()
    {
        if let Some((residual_curve, residual_var)) =
            pilot_residual_ifft_curve_with_cfg(cfg, &xeq_used, used_bins, pilot_bins, pref)
        {
            let residual_curve =
                temporal_ema_residual_curve(cfg, state, residual_curve, residual_var);
            for (z, corr) in xeq_used.iter_mut().zip(residual_curve.iter()) {
                if corr.norm() > 1.0e-3 && corr.re.is_finite() && corr.im.is_finite() {
                    *z /= *corr;
                }
            }
        }

        let mut num = Complex32::new(0.0, 0.0);
        let mut den = 0.0f32;
        for (k, pbin) in pilot_bins.iter().enumerate() {
            if let Some(pos) = used_bins.iter().position(|b| b == pbin) {
                num += xeq_used[pos] * pref[k].conj();
                den += pref[k].norm_sqr();
            }
        }
        if den > 1.0e-9 {
            let g = num / den;
            if g.norm() > 1.0e-6 && g.re.is_finite() && g.im.is_finite() {
                for z in &mut xeq_used {
                    *z /= g;
                }
            }
        }
        return xeq_used;
    }

    if cfg
        .equalizer
        .features
        .contains(EqualizerFeatures::PILOT_PHASE)
        && !pilot_bins.is_empty()
        && !pref.is_empty()
    {
        let mut pilot_phase_pts = Vec::<(f32, f32, f32)>::new();
        for (k, pbin) in pilot_bins.iter().enumerate() {
            if let Some(pos) = used_bins.iter().position(|b| b == pbin) {
                let ref_sym = pref[k];
                if ref_sym.norm_sqr() > 1.0e-9 {
                    let residual = xeq_used[pos] * ref_sym.conj();
                    pilot_phase_pts.push((
                        *pbin as f32,
                        residual.arg(),
                        pilot_weight(xeq_used[pos], ref_sym, cfg),
                    ));
                }
            }
        }
        if pilot_phase_pts.len() >= 2 {
            unwrap_phase_points(&mut pilot_phase_pts);
            if let Some((intercept, slope)) = fit_phase_line(&pilot_phase_pts) {
                for (k, &bin) in used_bins.iter().enumerate() {
                    let ph = intercept + slope * (bin as f32);
                    xeq_used[k] *= Complex32::from_polar(1.0, -ph);
                }
            }
        }
    }

    if cfg
        .equalizer
        .features
        .contains(EqualizerFeatures::PILOT_AMPLITUDE)
        && !pilot_bins.is_empty()
        && !pref.is_empty()
    {
        let mut pilot_amp_pts = Vec::<(f32, f32, f32)>::new();
        for (k, pbin) in pilot_bins.iter().enumerate() {
            if let Some(pos) = used_bins.iter().position(|b| b == pbin) {
                let ref_sym = pref[k];
                let ref_mag = ref_sym.norm();
                if ref_mag > 1.0e-6 {
                    let gain = (xeq_used[pos].norm() / ref_mag).clamp(0.5, 2.0);
                    pilot_amp_pts.push((
                        *pbin as f32,
                        gain,
                        pilot_weight(xeq_used[pos], ref_sym, cfg),
                    ));
                }
            }
        }
        if let Some((intercept, slope)) = fit_real_line(&pilot_amp_pts) {
            for (k, &bin) in used_bins.iter().enumerate() {
                let amp = (intercept + slope * (bin as f32)).clamp(0.5, 2.0);
                xeq_used[k] /= amp;
            }
        }
    }

    let mut num = Complex32::new(0.0, 0.0);
    let mut den = 0.0f32;
    for (k, pbin) in pilot_bins.iter().enumerate() {
        if let Some(pos) = used_bins.iter().position(|b| b == pbin) {
            num += xeq_used[pos] * pref[k].conj();
            den += pref[k].norm_sqr();
        }
    }
    if den > 1.0e-9 {
        let g = num / den;
        if g.norm() > 1.0e-6 && g.re.is_finite() && g.im.is_finite() {
            for z in &mut xeq_used {
                *z /= g;
            }
        }
    }
    xeq_used
}

/// Equalizes one OFDM symbol and applies pilot-derived residual correction.
///
/// Rationale:
/// The working instantaneous path is kept intact as a fallback because it is a
/// known-good operating mode. When temporal least-squares tracking is enabled,
/// only the low-dimensional pilot phase-line parameters are tracked over time;
/// all other stages remain local to the current symbol.
pub(crate) fn equalize_symbol_with_pilots(
    cfg: &OfdmConfig,
    state: &mut EqualizerTrackingState,
    y: &[Complex32],
    used_bins: &[usize],
    pilot_bins: &[usize],
    pref: &[Complex32],
    hest: &[Complex32],
) -> Vec<Complex32> {
    if !cfg
        .equalizer
        .features
        .contains(EqualizerFeatures::TEMPORAL_LS)
    {
        let xeq_used = equalize_symbol_with_pilots_instantaneous(
            cfg, state, y, used_bins, pilot_bins, pref, hest,
        );
        update_noise_variance(state, &xeq_used, used_bins, pilot_bins, pref);
        return xeq_used;
    }

    if used_bins.is_empty() || hest.len() != used_bins.len() {
        return Vec::new();
    }

    let mut xeq_used = equalize_symbol_with_baseline(cfg, state, y, used_bins, hest);

    if cfg
        .equalizer
        .features
        .contains(EqualizerFeatures::PILOT_PHASE)
        && !pilot_bins.is_empty()
        && !pref.is_empty()
    {
        let mut pilot_phase_pts = Vec::<(f32, f32, f32)>::new();
        for (k, pbin) in pilot_bins.iter().enumerate() {
            if let Some(pos) = used_bins.iter().position(|b| b == pbin) {
                let ref_sym = pref[k];
                if ref_sym.norm_sqr() > 1.0e-9 {
                    let residual = xeq_used[pos] * ref_sym.conj();
                    pilot_phase_pts.push((
                        *pbin as f32,
                        residual.arg(),
                        pilot_weight(xeq_used[pos], ref_sym, cfg),
                    ));
                }
            }
        }
        if pilot_phase_pts.len() >= 2 {
            unwrap_phase_points(&mut pilot_phase_pts);
            if let Some(current_line) = fit_phase_line(&pilot_phase_pts) {
                let phase_line = temporal_phase_line_fit(cfg, state, current_line);
                for (k, &bin) in used_bins.iter().enumerate() {
                    let ph = phase_line.0 + phase_line.1 * (bin as f32);
                    xeq_used[k] *= Complex32::from_polar(1.0, -ph);
                }
            }
        }
    }

    if cfg
        .equalizer
        .features
        .contains(EqualizerFeatures::PILOT_AMPLITUDE)
        && !pilot_bins.is_empty()
        && !pref.is_empty()
    {
        let mut pilot_amp_pts = Vec::<(f32, f32, f32)>::new();
        for (k, pbin) in pilot_bins.iter().enumerate() {
            if let Some(pos) = used_bins.iter().position(|b| b == pbin) {
                let ref_sym = pref[k];
                let ref_mag = ref_sym.norm();
                if ref_mag > 1.0e-6 {
                    let gain = (xeq_used[pos].norm() / ref_mag).clamp(0.5, 2.0);
                    pilot_amp_pts.push((
                        *pbin as f32,
                        gain,
                        pilot_weight(xeq_used[pos], ref_sym, cfg),
                    ));
                }
            }
        }
        if let Some((intercept, slope)) = fit_real_line(&pilot_amp_pts) {
            for (k, &bin) in used_bins.iter().enumerate() {
                let amp = (intercept + slope * (bin as f32)).clamp(0.5, 2.0);
                xeq_used[k] /= amp;
            }
        }
    }

    let mut num = Complex32::new(0.0, 0.0);
    let mut den = 0.0f32;
    for (k, pbin) in pilot_bins.iter().enumerate() {
        if let Some(pos) = used_bins.iter().position(|b| b == pbin) {
            num += xeq_used[pos] * pref[k].conj();
            den += pref[k].norm_sqr();
        }
    }
    if den > 1.0e-9 {
        let g = num / den;
        if g.norm() > 1.0e-6 && g.re.is_finite() && g.im.is_finite() {
            for z in &mut xeq_used {
                *z /= g;
            }
        }
    }

    update_noise_variance(state, &xeq_used, used_bins, pilot_bins, pref);
    xeq_used
}

fn temporal_phase_line_fit(
    cfg: &OfdmConfig,
    state: &mut EqualizerTrackingState,
    current: (f32, f32),
) -> (f32, f32) {
    let current_time = state.next_symbol_time;
    let intercept = if let Some((_, prev_i, _)) = state.phase_line_hist.back() {
        normalize_angle_near(current.0, *prev_i)
    } else {
        current.0
    };
    state
        .phase_line_hist
        .push_back((current_time, intercept, current.1));
    state.next_symbol_time += 1.0;

    let keep = cfg.equalizer.temporal_window.max(1);
    while state.phase_line_hist.len() > keep {
        state.phase_line_hist.pop_front();
    }

    if state.phase_line_hist.len() < 2 {
        return (intercept, current.1);
    }

    let mut intercept_pts = Vec::with_capacity(state.phase_line_hist.len());
    let mut slope_pts = Vec::with_capacity(state.phase_line_hist.len());
    for &(t, a, b) in &state.phase_line_hist {
        intercept_pts.push((t, a, 1.0));
        slope_pts.push((t, b, 1.0));
    }
    let fit_a = fit_real_line(&intercept_pts).map(|(c0, c1)| c0 + c1 * current_time);
    let fit_b = fit_real_line(&slope_pts).map(|(d0, d1)| d0 + d1 * current_time);
    match (fit_a, fit_b) {
        (Some(a), Some(b)) => (a, b),
        _ => (intercept, current.1),
    }
}

fn update_noise_variance(
    state: &mut EqualizerTrackingState,
    xeq_used: &[Complex32],
    used_bins: &[usize],
    pilot_bins: &[usize],
    pref: &[Complex32],
) {
    if pilot_bins.is_empty() || pref.is_empty() {
        return;
    }
    let mut err_sum = 0.0f32;
    let mut count = 0usize;
    for (k, pbin) in pilot_bins.iter().enumerate() {
        if let Some(pos) = used_bins.iter().position(|b| b == pbin) {
            err_sum += (xeq_used[pos] - pref[k]).norm_sqr();
            count += 1;
        }
    }
    if count > 0 {
        state.noise_var = Some((err_sum / count as f32).clamp(1.0e-4, 2.5e-1));
    }
}

fn normalize_angle_near(mut value: f32, reference: f32) -> f32 {
    while value - reference > std::f32::consts::PI {
        value -= 2.0 * std::f32::consts::PI;
    }
    while value - reference < -std::f32::consts::PI {
        value += 2.0 * std::f32::consts::PI;
    }
    value
}

/// Builds a smooth baseline channel model from the training estimate.
///
/// Rationale:
/// The raw training LS estimate is too noisy to trust literally on magnitude.
/// We therefore keep phase as a per-bin baseline, but smooth the magnitude
/// across neighboring bins:
///
/// $$\hat A_0[k] \leftarrow \mathrm{smooth}(|\hat H_{\text{train}}[k]|)$$
///
/// which reduces weak-bin overreaction before pilots apply per-symbol phase
/// correction.
fn baseline_channel_model(hest: &[Complex32]) -> (Vec<f32>, Vec<f32>) {
    let amps = smooth_real_line(&hest.iter().map(|h| h.norm()).collect::<Vec<_>>());
    let phases = unwrap_phases(&hest.iter().map(|h| h.arg()).collect::<Vec<_>>());
    (amps, phases)
}

/// Fits a residual pilot phase line over bin index.
///
/// Rationale:
/// The dominant residual we keep seeing is timing-like phase slope across bins.
/// A least-squares line,
///
/// $$\Delta \phi(k) \approx a + b k$$
///
/// is therefore a better constrained correction than independent per-bin pilot
/// updates.
fn fit_phase_line(pts: &[(f32, f32, f32)]) -> Option<(f32, f32)> {
    if pts.len() < 2 {
        return None;
    }
    let sw = pts.iter().map(|(_, _, w)| *w).sum::<f32>();
    let sx = pts.iter().map(|(x, _, w)| x * w).sum::<f32>();
    let sy = pts.iter().map(|(_, y, w)| y * w).sum::<f32>();
    let sxx = pts.iter().map(|(x, _, w)| x * x * w).sum::<f32>();
    let sxy = pts.iter().map(|(x, y, w)| x * y * w).sum::<f32>();
    let denom = sw * sxx - sx * sx;
    if denom.abs() <= 1.0e-9 {
        return None;
    }
    let slope = (sw * sxy - sx * sy) / denom;
    let intercept = (sy - slope * sx) / sw;
    Some((intercept, slope))
}

/// Fits a smooth real-valued residual over bin index.
///
/// Rationale:
/// Residual pilot amplitude tends to vary slowly across the active band, so a
/// first-order fit is a safer correction than per-bin inversion:
///
/// $$\Delta A(k) \approx c + d k$$
///
/// The fit is intentionally low-order and later clamped to avoid noisy gain
/// excursions on weak pilots.
fn fit_real_line(pts: &[(f32, f32, f32)]) -> Option<(f32, f32)> {
    if pts.is_empty() {
        return None;
    }
    if pts.len() == 1 {
        return Some((pts[0].1, 0.0));
    }
    let sw = pts.iter().map(|(_, _, w)| *w).sum::<f32>();
    let sx = pts.iter().map(|(x, _, w)| x * w).sum::<f32>();
    let sy = pts.iter().map(|(_, y, w)| y * w).sum::<f32>();
    let sxx = pts.iter().map(|(x, _, w)| x * x * w).sum::<f32>();
    let sxy = pts.iter().map(|(x, y, w)| x * y * w).sum::<f32>();
    let denom = sw * sxx - sx * sx;
    if denom.abs() <= 1.0e-9 {
        return Some((sy / sw.max(1.0e-9), 0.0));
    }
    let slope = (sw * sxy - sx * sy) / denom;
    let intercept = (sy - slope * sx) / sw;
    Some((intercept, slope))
}

/// Computes RMS EVM against a known reference constellation.
///
/// Rationale:
/// This reports
///
/// $$\mathrm{EVM}_{\mathrm{RMS}} = \sqrt{\frac{1}{N}\sum_k |\hat X[k]-X[k]|^2}$$
///
/// which is the compact scalar diagnostic used throughout the modem logs.
pub(crate) fn rms_evm(got: &[Complex32], want: &[Complex32]) -> f32 {
    if got.is_empty() || got.len() != want.len() {
        return 0.0;
    }
    let mse = got
        .iter()
        .zip(want.iter())
        .map(|(g, w)| (*g - *w).norm_sqr())
        .sum::<f32>()
        / (got.len() as f32);
    mse.sqrt()
}

/// Hard-slices noisy symbols onto the ideal constellation.
///
/// Rationale:
/// This gives the nearest nominal point used by decision-directed metrics and
/// diagnostics without altering the actual decoder decisions elsewhere.
pub(crate) fn hard_slice_symbols(syms: &[Complex32], modulation: Modulation) -> Vec<Complex32> {
    match modulation {
        Modulation::Bpsk => syms
            .iter()
            .map(|s| {
                if s.re < 0.0 {
                    Complex32::new(-1.0, 0.0)
                } else {
                    Complex32::new(1.0, 0.0)
                }
            })
            .collect(),
        Modulation::Qpsk => syms
            .iter()
            .map(|s| {
                let re = if s.re < 0.0 { -1.0 } else { 1.0 };
                let im = if s.im < 0.0 { -1.0 } else { 1.0 };
                Complex32::new(re, im) * (1.0 / 2.0f32.sqrt())
            })
            .collect(),
    }
}

/// Computes decision-directed EVM.
///
/// Rationale:
/// When the true transmitted data symbols are unknown, we compare against the
/// sliced constellation,
///
/// $$X_{\text{ref}}[k] = \mathcal{Q}(\hat X[k])$$
///
/// to get a pragmatic quality metric for post-equalized data bins.
pub(crate) fn decision_directed_evm(syms: &[Complex32], modulation: Modulation) -> f32 {
    if syms.is_empty() {
        return 0.0;
    }
    let refs = hard_slice_symbols(syms, modulation);
    rms_evm(syms, &refs)
}

fn unwrap_phase_points(pts: &mut [(f32, f32, f32)]) {
    for i in 1..pts.len() {
        let mut phi = pts[i].1;
        let prev = pts[i - 1].1;
        while phi - prev > std::f32::consts::PI {
            phi -= 2.0 * std::f32::consts::PI;
        }
        while phi - prev < -std::f32::consts::PI {
            phi += 2.0 * std::f32::consts::PI;
        }
        pts[i].1 = phi;
    }
}

fn pilot_weight(z: Complex32, pref: Complex32, cfg: &OfdmConfig) -> f32 {
    if !cfg
        .equalizer
        .features
        .contains(EqualizerFeatures::WEIGHTED_PILOTS)
    {
        return 1.0;
    }
    let err = (z - pref).norm_sqr();
    (1.0 / (0.05 + err)).clamp(0.1, 10.0)
}

/// Unwraps a phase sequence over bin index.
///
/// Rationale:
/// Wrapped phase jumps at $\pm \pi$ are artificial. Unwrapping restores the
/// continuous phase trend needed for fitting or smoothing across bins.
fn unwrap_phases(phases: &[f32]) -> Vec<f32> {
    if phases.is_empty() {
        return Vec::new();
    }
    let mut out = phases.to_vec();
    for i in 1..out.len() {
        let mut phi = out[i];
        let prev = out[i - 1];
        while phi - prev > std::f32::consts::PI {
            phi -= 2.0 * std::f32::consts::PI;
        }
        while phi - prev < -std::f32::consts::PI {
            phi += 2.0 * std::f32::consts::PI;
        }
        out[i] = phi;
    }
    out
}

/// Smooths a real-valued line with a short symmetric kernel.
///
/// Rationale:
/// Training magnitude is useful as a baseline, but noisy per-bin inverses are
/// fragile. A local average keeps the gross envelope while avoiding aggressive
/// gain excursions on isolated weak bins.
fn smooth_real_line(x: &[f32]) -> Vec<f32> {
    if x.is_empty() {
        return Vec::new();
    }
    let mut out = vec![0.0f32; x.len()];
    for i in 0..x.len() {
        let i0 = i.saturating_sub(1);
        let i1 = (i + 1).min(x.len() - 1);
        let mut sum = 0.0f32;
        let mut wsum = 0.0f32;
        for j in i0..=i1 {
            let w = if j == i { 0.5 } else { 0.25 };
            sum += w * x[j];
            wsum += w;
        }
        out[i] = (sum / wsum).max(1.0e-3);
    }
    out
}

fn pilot_residual_ifft_curve_with_cfg(
    cfg: &OfdmConfig,
    xeq_used: &[Complex32],
    used_bins: &[usize],
    pilot_bins: &[usize],
    pref: &[Complex32],
) -> Option<(Vec<Complex32>, f32)> {
    let pts = pilot_residual_points(xeq_used, used_bins, pilot_bins, pref);
    if pts.len() < 2 {
        return None;
    }

    let noise_var = pilot_residual_noise_var(&pts);
    let selected_order = match cfg.equalizer.residual_tap_order_mode {
        ResidualTapOrderMode::All => None,
        ResidualTapOrderMode::Fixed => Some(
            cfg.equalizer
                .residual_tap_order
                .max(1)
                .min(used_bins.len())
                .min(pts.len()),
        ),
        ResidualTapOrderMode::Mdl => Some(select_delay_tap_order_mdl(
            &pts,
            used_bins.len(),
            cfg.equalizer.residual_tap_order_max,
        )?),
    };

    let mut taps = match selected_order {
        None => {
            let curve = interpolate_residual_curve(&pts, used_bins.len());
            let mut taps = curve.clone();
            unitary_ifft_in_place(&mut taps);
            taps
        }
        Some(order) => fit_delay_taps_least_squares(&pts, used_bins.len(), order)?,
    };

    if let Some(keep) = selected_order {
        let keep = keep.min(taps.len()).min(pts.len());
        for tap in taps.iter_mut().skip(keep) {
            *tap = Complex32::new(0.0, 0.0);
        }
    }

    for tap in taps.iter_mut() {
        let power = tap.norm_sqr();
        let shrink = ((power - noise_var) / (power + 1.0e-9)).clamp(0.0, 1.0);
        *tap *= shrink;
    }

    let mut curve = taps.clone();
    unitary_fft_in_place(&mut curve);
    for v in &mut curve {
        let mag = v.norm().clamp(0.75, 1.5);
        let phase = v.arg();
        *v = Complex32::from_polar(mag, phase);
    }
    Some((curve, noise_var))
}

fn pilot_residual_points(
    xeq_used: &[Complex32],
    used_bins: &[usize],
    pilot_bins: &[usize],
    pref: &[Complex32],
) -> Vec<(usize, Complex32)> {
    let mut pts = Vec::<(usize, Complex32)>::new();
    for (k, pbin) in pilot_bins.iter().enumerate() {
        if let Some(pos) = used_bins.iter().position(|b| b == pbin) {
            let ref_sym = pref[k];
            if ref_sym.norm_sqr() > 1.0e-9 {
                pts.push((pos, xeq_used[pos] / ref_sym));
            }
        }
    }
    pts
}

fn interpolate_residual_curve(pts: &[(usize, Complex32)], curve_len: usize) -> Vec<Complex32> {
    let mut curve = vec![Complex32::new(1.0, 0.0); curve_len];
    for i in 0..curve.len() {
        if i <= pts[0].0 {
            curve[i] = pts[0].1;
            continue;
        }
        if i >= pts[pts.len() - 1].0 {
            curve[i] = pts[pts.len() - 1].1;
            continue;
        }
        let mut seg = 0usize;
        while seg + 1 < pts.len() && i > pts[seg + 1].0 {
            seg += 1;
        }
        let (i0, z0) = pts[seg];
        let (i1, z1) = pts[seg + 1];
        let t = ((i - i0) as f32 / (i1 - i0).max(1) as f32).clamp(0.0, 1.0);
        curve[i] = z0 * (1.0 - t) + z1 * t;
    }
    curve
}

/// Fits a leading-tap residual channel model directly in pilot space.
///
/// Rationale:
/// The pilot residual denoiser lives in a delay-domain basis, so model-order
/// selection should evaluate candidate tap counts in that same basis rather
/// than after a heuristic interpolation step.
fn fit_delay_taps_least_squares(
    pts: &[(usize, Complex32)],
    curve_len: usize,
    order: usize,
) -> Option<Vec<Complex32>> {
    if pts.is_empty() || curve_len == 0 {
        return None;
    }

    let order = order.max(1).min(curve_len).min(pts.len());
    let mut gram = vec![vec![Complex32::new(0.0, 0.0); order]; order];
    let mut rhs = vec![Complex32::new(0.0, 0.0); order];
    for &(pos, obs) in pts {
        let basis = delay_basis_row(curve_len, pos, order);
        for i in 0..order {
            rhs[i] += basis[i].conj() * obs;
            for j in 0..order {
                gram[i][j] += basis[i].conj() * basis[j];
            }
        }
    }
    for (i, row) in gram.iter_mut().enumerate() {
        row[i] += Complex32::new(1.0e-6, 0.0);
    }

    let sol = solve_complex_linear_system(gram, rhs)?;
    let mut taps = vec![Complex32::new(0.0, 0.0); curve_len];
    for (tap, value) in taps.iter_mut().zip(sol.iter()) {
        *tap = *value;
    }
    Some(taps)
}

fn select_delay_tap_order_mdl(
    pts: &[(usize, Complex32)],
    curve_len: usize,
    max_order: usize,
) -> Option<usize> {
    if pts.is_empty() || curve_len == 0 {
        return None;
    }

    let sample_count = (2 * pts.len()).max(1) as f32;
    let max_order = max_order.max(1).min(curve_len).min(pts.len());
    let mut best_order = 1usize;
    let mut best_score = f32::INFINITY;
    for order in 1..=max_order {
        let taps = fit_delay_taps_least_squares(pts, curve_len, order)?;
        let mse = delay_tap_fit_mse(pts, curve_len, &taps).clamp(1.0e-6, 1.0e3);
        let n_params = 2.0 * order as f32;
        let score = sample_count * mse.ln() + n_params * sample_count.ln();
        if score < best_score {
            best_score = score;
            best_order = order;
        }
    }
    Some(best_order)
}

fn delay_tap_fit_mse(pts: &[(usize, Complex32)], curve_len: usize, taps: &[Complex32]) -> f32 {
    if pts.is_empty() {
        return 1.0;
    }

    let mut err = 0.0f32;
    for &(pos, obs) in pts {
        let est = synthesize_delay_response_at(curve_len, taps, pos);
        err += (obs - est).norm_sqr();
    }
    err / pts.len() as f32
}

fn delay_basis_row(curve_len: usize, pos: usize, order: usize) -> Vec<Complex32> {
    let scale = 1.0 / (curve_len as f32).sqrt();
    (0..order)
        .map(|tap| {
            let phase = -2.0 * PI * (pos as f32) * (tap as f32) / curve_len as f32;
            Complex32::from_polar(scale, phase)
        })
        .collect()
}

fn synthesize_delay_response_at(curve_len: usize, taps: &[Complex32], pos: usize) -> Complex32 {
    let scale = 1.0 / (curve_len as f32).sqrt();
    taps.iter()
        .enumerate()
        .fold(Complex32::new(0.0, 0.0), |acc, (tap, value)| {
            let phase = -2.0 * PI * (pos as f32) * (tap as f32) / curve_len as f32;
            acc + *value * Complex32::from_polar(scale, phase)
        })
}

fn solve_complex_linear_system(
    mut a: Vec<Vec<Complex32>>,
    mut b: Vec<Complex32>,
) -> Option<Vec<Complex32>> {
    let n = b.len();
    for i in 0..n {
        let mut pivot = i;
        let mut pivot_norm = a[i][i].norm_sqr();
        for (row_idx, row) in a.iter().enumerate().skip(i + 1) {
            let cand = row[i].norm_sqr();
            if cand > pivot_norm {
                pivot = row_idx;
                pivot_norm = cand;
            }
        }
        if pivot_norm <= 1.0e-12 {
            return None;
        }
        if pivot != i {
            a.swap(i, pivot);
            b.swap(i, pivot);
        }

        let diag = a[i][i];
        for col in i..n {
            a[i][col] /= diag;
        }
        b[i] /= diag;

        let pivot_row = a[i].clone();
        let pivot_rhs = b[i];
        for row in 0..n {
            if row == i {
                continue;
            }
            let factor = a[row][i];
            if factor.norm_sqr() <= 1.0e-18 {
                continue;
            }
            for col in i..n {
                a[row][col] -= factor * pivot_row[col];
            }
            b[row] -= factor * pivot_rhs;
        }
    }
    Some(b)
}

fn unitary_fft_in_place(x: &mut [Complex32]) {
    if x.is_empty() {
        return;
    }
    let mut planner = FftPlanner::<f32>::new();
    let fft = planner.plan_fft_forward(x.len());
    fft.process(x);
    let scale = 1.0 / (x.len() as f32).sqrt();
    for v in x {
        *v *= scale;
    }
}

fn unitary_ifft_in_place(x: &mut [Complex32]) {
    if x.is_empty() {
        return;
    }
    let mut planner = FftPlanner::<f32>::new();
    let ifft = planner.plan_fft_inverse(x.len());
    ifft.process(x);
    let scale = 1.0 / (x.len() as f32).sqrt();
    for v in x {
        *v *= scale;
    }
}

fn pilot_residual_noise_var(pts: &[(usize, Complex32)]) -> f32 {
    if pts.len() < 3 {
        return 5.0e-3;
    }

    let mut err_sum = 0.0f32;
    let mut count = 0usize;
    for win in pts.windows(3) {
        let (i0, z0) = win[0];
        let (i1, z1) = win[1];
        let (i2, z2) = win[2];
        let denom = (i2 - i0).max(1) as f32;
        let t = ((i1 - i0) as f32 / denom).clamp(0.0, 1.0);
        let pred = z0 * (1.0 - t) + z2 * t;
        err_sum += (z1 - pred).norm_sqr();
        count += 1;
    }

    if count == 0 {
        5.0e-3
    } else {
        (err_sum / count as f32).clamp(1.0e-4, 2.5e-1)
    }
}

fn temporal_ema_residual_curve(
    cfg: &OfdmConfig,
    state: &mut EqualizerTrackingState,
    current: Vec<Complex32>,
    _current_var: f32,
) -> Vec<Complex32> {
    if !cfg
        .equalizer
        .features
        .contains(EqualizerFeatures::TEMPORAL_RESIDUAL_EMA)
    {
        return current;
    }

    let keep = cfg.equalizer.temporal_window.max(1);
    if keep <= 1 {
        state.residual_curve_est = Some(current.clone());
        return current;
    }

    let alpha = (1.0 / keep as f32).clamp(0.05, 1.0);
    let fused = match &state.residual_curve_est {
        Some(prev) if prev.len() == current.len() => {
            let phase_align = common_phase_delta(prev, &current);
            let rot = Complex32::from_polar(1.0, -phase_align);
            let mut out = current.clone();
            for ((dst, p), c) in out.iter_mut().zip(prev.iter()).zip(current.iter()) {
                let c_aligned = *c * rot;
                *dst = c_aligned * alpha + *p * (1.0 - alpha);
            }
            out
        }
        _ => current.clone(),
    };
    state.residual_curve_est = Some(fused.clone());
    fused
}

fn common_phase_delta(reference: &[Complex32], current: &[Complex32]) -> f32 {
    let mut acc = Complex32::new(0.0, 0.0);
    for (r, c) in reference.iter().zip(current.iter()) {
        let rn = r.norm();
        let cn = c.norm();
        if rn > 1.0e-6 && cn > 1.0e-6 {
            let ru = *r / rn;
            let cu = *c / cn;
            acc += cu * ru.conj();
        }
    }
    if acc.norm() > 1.0e-9 {
        acc.arg()
    } else {
        0.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn synthesize_points(
        curve_len: usize,
        taps: &[Complex32],
        pilot_pos: &[usize],
    ) -> Vec<(usize, Complex32)> {
        pilot_pos
            .iter()
            .map(|&pos| (pos, synthesize_delay_response_at(curve_len, taps, pos)))
            .collect()
    }

    #[test]
    fn mdl_prefers_two_tap_model_when_data_is_two_tap() {
        let curve_len = 8usize;
        let pilot_pos = [0usize, 1, 2, 4, 5, 6, 7];
        let mut taps = vec![Complex32::new(0.0, 0.0); curve_len];
        taps[0] = Complex32::new(1.0, 0.0);
        taps[1] = Complex32::new(0.2, -0.15);
        let pts = synthesize_points(curve_len, &taps, &pilot_pos);

        let order = select_delay_tap_order_mdl(&pts, curve_len, 6).expect("mdl order");
        assert_eq!(order, 2);
    }

    #[test]
    fn least_squares_delay_fit_matches_pilot_observations() {
        let curve_len = 8usize;
        let pilot_pos = [0usize, 2, 3, 4, 6, 7];
        let mut taps = vec![Complex32::new(0.0, 0.0); curve_len];
        taps[0] = Complex32::new(1.0, 0.0);
        taps[1] = Complex32::new(-0.1, 0.25);
        taps[2] = Complex32::new(0.05, -0.08);
        let pts = synthesize_points(curve_len, &taps, &pilot_pos);

        let fit = fit_delay_taps_least_squares(&pts, curve_len, 3).expect("ls fit");
        let mse = delay_tap_fit_mse(&pts, curve_len, &fit);
        assert!(mse < 1.0e-8, "mse={mse}");
    }
}

// vim: set ts=4 sw=4 et:
