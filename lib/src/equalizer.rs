// Copyright (c) 2026 Elias S. G. Carotti

use rustfft::num_complex::Complex32;

use crate::config::{EqualizerFeatures, Modulation, OfdmConfig};

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

/// Equalizes one OFDM symbol against a baseline channel model.
///
/// Rationale:
/// We explicitly separate the training estimate into amplitude and phase,
///
/// $$\hat H_0[k] = \hat A_0[k] e^{j \hat \phi_0[k]}$$
///
/// and invert that baseline before applying any pilot-derived residual.
fn equalize_symbol_with_baseline(
    y: &[Complex32],
    used_bins: &[usize],
    hest: &[Complex32],
) -> Vec<Complex32> {
    let (amp0, phase0) = baseline_channel_model(hest);
    let mut xeq_used = Vec::with_capacity(used_bins.len());
    for (k, &bin) in used_bins.iter().enumerate() {
        let h0 = Complex32::from_polar(amp0[k], phase0[k]);
        xeq_used.push(regularized_equalize(y[bin], h0));
    }
    xeq_used
}

/// Equalizes one OFDM symbol and applies pilot-derived residual correction.
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
pub(crate) fn equalize_symbol_with_pilots(
    cfg: &OfdmConfig,
    y: &[Complex32],
    used_bins: &[usize],
    pilot_bins: &[usize],
    pref: &[Complex32],
    hest: &[Complex32],
) -> Vec<Complex32> {
    if used_bins.is_empty() || hest.len() != used_bins.len() {
        return Vec::new();
    }

    let mut xeq_used = equalize_symbol_with_baseline(y, used_bins, hest);

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

// vim: set ts=4 sw=4 et:
