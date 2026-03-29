// Copyright (c) 2026 Elias S. G. Carotti

use rustfft::num_complex::Complex32;

use crate::config::{EqualizerMode, Modulation, OfdmConfig};

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
/// `TrainingPilot` uses the training estimate as the static baseline,
/// while `PilotOnly` starts from a flat unit channel,
///
/// $$\hat H_0[k] = 1$$
///
/// and lets later pilot corrections carry the residual alignment.
pub(crate) fn equalizer_initial_channel(
    cfg: &OfdmConfig,
    ytrain: &[Complex32],
    used_bins: &[usize],
    train_known: &[Complex32],
) -> Vec<Complex32> {
    match cfg.equalizer_mode {
        EqualizerMode::TrainingPilot => {
            estimate_channel_from_training(ytrain, used_bins, train_known)
        }
        EqualizerMode::PilotOnly => vec![Complex32::new(1.0, 0.0); used_bins.len()],
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
    if matches!(cfg.equalizer_mode, EqualizerMode::TrainingPilot) {
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

/// Equalizes one OFDM symbol and applies a pilot-derived phase model.
///
/// Rationale:
/// The training symbol gives the static baseline, while pilots capture the
/// symbol-local residual phase trend. The current model fits a linear phase
/// residual over bin index,
///
/// $$\Delta \phi(k) = a + b k$$
///
/// and derotates all used bins by that trend before a final common pilot
/// normalization.
pub(crate) fn equalize_symbol_with_pilots(
    y: &[Complex32],
    used_bins: &[usize],
    pilot_bins: &[usize],
    pref: &[Complex32],
    hest: &[Complex32],
) -> Vec<Complex32> {
    if used_bins.is_empty() || hest.len() != used_bins.len() {
        return Vec::new();
    }

    let mut xeq_used = Vec::with_capacity(used_bins.len());
    for (k, &bin) in used_bins.iter().enumerate() {
        xeq_used.push(regularized_equalize(y[bin], hest[k]));
    }

    if !pilot_bins.is_empty() && !pref.is_empty() {
        let mut pilot_phase_pts = Vec::<(f32, f32)>::new();
        for (k, pbin) in pilot_bins.iter().enumerate() {
            if let Some(pos) = used_bins.iter().position(|b| b == pbin) {
                let ref_sym = pref[k];
                if ref_sym.norm_sqr() > 1.0e-9 {
                    pilot_phase_pts.push((*pbin as f32, (xeq_used[pos] * ref_sym.conj()).arg()));
                }
            }
        }
        if pilot_phase_pts.len() >= 2 {
            unwrap_phase_points(&mut pilot_phase_pts);

            let n = pilot_phase_pts.len() as f32;
            let sx = pilot_phase_pts.iter().map(|(x, _)| *x).sum::<f32>();
            let sy = pilot_phase_pts.iter().map(|(_, y)| *y).sum::<f32>();
            let sxx = pilot_phase_pts.iter().map(|(x, _)| x * x).sum::<f32>();
            let sxy = pilot_phase_pts.iter().map(|(x, y)| x * y).sum::<f32>();
            let denom = n * sxx - sx * sx;
            if denom.abs() > 1.0e-9 {
                let slope = (n * sxy - sx * sy) / denom;
                let intercept = (sy - slope * sx) / n;
                for (k, &bin) in used_bins.iter().enumerate() {
                    let ph = intercept + slope * (bin as f32);
                    xeq_used[k] *= Complex32::from_polar(1.0, -ph);
                }
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

fn unwrap_phase_points(pts: &mut [(f32, f32)]) {
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

// vim: set ts=4 sw=4 et:
