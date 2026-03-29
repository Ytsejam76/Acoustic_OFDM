// Copyright (c) 2026 Elias S. G. Carotti

use rustfft::num_complex::Complex32;

use crate::baseband::{fft, known_pilot_symbols, known_training_symbols, ofdm_bin_plan};
use crate::config::OfdmConfig;
use crate::equalizer::{equalize_symbol_with_pilots, regularized_equalize, rms_evm};

/// Returns the active DSP sample rate seen by the modem.
///
/// Rationale:
/// `legacy` and `iq` paths operate at different internal rates. Any timing or
/// CFO estimator must use the rate of the signal it actually sees.
pub(crate) fn active_baseband_fs(cfg: &OfdmConfig) -> f32 {
    match cfg.passband_mode {
        crate::config::PassbandMode::Legacy => cfg.fs,
        crate::config::PassbandMode::Iq => cfg.fs_baseband,
    }
}

/// Samples a complex sequence with linear interpolation.
///
/// Rationale:
/// Fractional timing offsets are inevitable. We model them by evaluating the
/// sequence at non-integer positions,
///
/// $$x(n+\delta) \approx (1-\alpha)x[i_0] + \alpha x[i_1]$$
///
/// with $\alpha = n+\delta - i_0$.
pub(crate) fn sample_complex_linear(x: &[Complex32], pos: f32) -> Complex32 {
    if x.is_empty() || pos < 0.0 {
        return Complex32::new(0.0, 0.0);
    }
    let i0 = pos.floor() as usize;
    if i0 >= x.len() {
        return Complex32::new(0.0, 0.0);
    }
    let i1 = (i0 + 1).min(x.len() - 1);
    let a = pos - (i0 as f32);
    x[i0] * (1.0 - a) + x[i1] * a
}

/// Resamples a sequence from a fractional start and step.
///
/// Rationale:
/// This is the simplest way to express both a fixed sync offset and a small
/// sample-clock scale error,
///
/// $$x_k = x(t_0 + k\,\Delta t)$$
///
/// where `start = t_0` and `step = \Delta t`.
pub(crate) fn resample_from_offset_rate(x: &[Complex32], start: f32, step: f32) -> Vec<Complex32> {
    if x.is_empty() {
        return Vec::new();
    }
    if !step.is_finite() || step <= 0.0 {
        return Vec::new();
    }
    let start0 = start.max(0.0);
    let remain = (x.len() as f32 - start0).max(0.0);
    let n = (remain / step).floor().max(0.0) as usize;
    let mut out = Vec::with_capacity(n);
    for k in 0..n {
        out.push(sample_complex_linear(x, start0 + (k as f32) * step));
    }
    out
}

/// Convenience wrapper for unit-rate resampling.
pub(crate) fn resample_from_offset(x: &[Complex32], start: f32) -> Vec<Complex32> {
    resample_from_offset_rate(x, start, 1.0)
}

/// Locates the repeated-half preamble using a Schmidl-Cox metric.
///
/// Rationale:
/// For a repeated preamble of length $L$, the timing metric is
///
/// $$M(d)=\frac{|\sum_{n=0}^{L-1}r^*[d+n]r[d+n+L]|^2}{(\sum_{n=0}^{L-1}|r[d+n+L]|^2)^2}$$
///
/// and the first strong plateau onset is used as the coarse packet start.
pub(crate) fn find_repeated_half_sync_offset(rbb: &[Complex32], cfg: &OfdmConfig) -> usize {
    let metrics = repeated_half_sync_metrics(rbb, cfg);
    if metrics.is_empty() {
        return 0;
    }
    let best_m = metrics.iter().copied().fold(-1.0f32, f32::max);
    let thresh = 0.97 * best_m.max(0.0);
    metrics.iter().position(|&m| m >= thresh).unwrap_or(0)
}

/// Computes the Schmidl-Cox repeated-half metric over the capture window.
pub(crate) fn repeated_half_sync_metrics(rbb: &[Complex32], cfg: &OfdmConfig) -> Vec<f32> {
    let l = cfg.sync_half_len;
    if l == 0 || rbb.len() < 2 * l + 2 {
        return Vec::new();
    }
    let max_search = ((0.12 * active_baseband_fs(cfg)).round() as usize)
        .min(rbb.len().saturating_sub(2 * l + 1));
    if max_search == 0 {
        return Vec::new();
    }

    let n_terms = max_search + l;
    let mut pref_v = vec![Complex32::new(0.0, 0.0); n_terms + 1];
    let mut pref_e = vec![0.0f32; n_terms + 1];
    for n in 0..n_terms {
        let v = rbb[n].conj() * rbb[n + l];
        let e = rbb[n + l].norm_sqr();
        pref_v[n + 1] = pref_v[n] + v;
        pref_e[n + 1] = pref_e[n] + e;
    }

    let mut metrics = vec![0.0f32; max_search + 1];
    for d in 0..=max_search {
        let p = pref_v[d + l] - pref_v[d];
        let r = (pref_e[d + l] - pref_e[d]).max(1e-9);
        metrics[d] = p.norm_sqr() / (r * r);
    }
    metrics
}

/// Refines the coarse sync offset around the repeated-half estimate.
///
/// Rationale:
/// The Schmidl-Cox plateau is broad. We therefore rescore a small fractional
/// neighborhood with a training-aware quality objective and keep the best
/// offset.
pub(crate) fn refine_sync_offset(rbb: &[Complex32], cfg: &OfdmConfig, coarse_off: usize) -> f32 {
    let search = (cfg.ncp / 8).clamp(2, 8) as i32;
    let mut best_off = coarse_off as f32;
    let mut best_score = sync_quality_score_fractional(rbb, cfg, best_off);
    for di in -search..=search {
        let base = (coarse_off as i32 + di).max(0) as f32;
        for q in 0..4 {
            let off = base + 0.25 * (q as f32);
            let score = sync_quality_score_fractional(rbb, cfg, off);
            if score > best_score + 1e-4 {
                best_score = score;
                best_off = off;
            }
        }
    }
    best_off
}

/// Scores the cyclic-prefix alignment of the training symbol.
///
/// Rationale:
/// A correct FFT window should align CP and symbol tail,
///
/// $$S_{\mathrm{CP}} = \frac{|\sum c^*[n]t[n]|}{\sqrt{\sum|c[n]|^2\sum|t[n]|^2}}$$
///
/// giving a timing cue that is narrower than the repeated-half plateau.
pub(crate) fn training_cp_score_fractional(rbb: &[Complex32], cfg: &OfdmConfig, off: f32) -> f32 {
    let xsync_len = 2 * cfg.sync_half_len;
    let train_end = xsync_len + cfg.ncp + cfg.nfft;
    let need = off.ceil().max(0.0) as usize + train_end;
    if rbb.len() < need || cfg.ncp == 0 || cfg.nfft == 0 {
        return -1.0;
    }
    let mut num = Complex32::new(0.0, 0.0);
    let mut e1 = 0.0f32;
    let mut e2 = 0.0f32;
    for k in 0..cfg.ncp {
        let cp = sample_complex_linear(rbb, off + xsync_len as f32 + k as f32);
        let tail = sample_complex_linear(rbb, off + xsync_len as f32 + cfg.nfft as f32 + k as f32);
        num += cp.conj() * tail;
        e1 += cp.norm_sqr();
        e2 += tail.norm_sqr();
    }
    let den = (e1 * e2).sqrt().max(1e-9);
    num.norm() / den
}

/// Scores a fractional timing hypothesis.
///
/// Rationale:
/// We combine CP consistency, training quality, and first-symbol pilot quality
/// into a practical score that prefers windows with cleaner downstream
/// equalization.
pub(crate) fn sync_quality_score_fractional(rbb: &[Complex32], cfg: &OfdmConfig, off: f32) -> f32 {
    let rbb_sync = resample_from_offset(rbb, off);
    let rbb_cfo = coarse_cfo_correct(&rbb_sync, cfg);
    let (used_bins, pilot_bins, _data_bins) = ofdm_bin_plan(cfg);
    let xsync_len = 2 * cfg.sync_half_len;
    let train_len = cfg.nfft + cfg.ncp;
    let need = xsync_len + train_len;
    if rbb_cfo.len() < need || used_bins.is_empty() {
        return -1.0;
    }

    let cp_score = training_cp_score_fractional(rbb, cfg, off);
    let train_start = xsync_len;
    let train_no_cp = &rbb_cfo[train_start + cfg.ncp..train_start + cfg.ncp + cfg.nfft];
    let ytrain = fft(train_no_cp);
    let train_known = known_training_symbols(used_bins.len(), cfg.modulation);
    let mut hest = vec![Complex32::new(1.0, 0.0); used_bins.len()];
    let mut ytrain_eq = Vec::with_capacity(used_bins.len());
    let mut hmag_sum = 0.0f32;
    for (k, &bin) in used_bins.iter().enumerate() {
        let h = ytrain[bin] / train_known[k];
        hest[k] = h;
        hmag_sum += h.norm();
        ytrain_eq.push(regularized_equalize(ytrain[bin], h));
    }
    let train_evm = rms_evm(&ytrain_eq, &train_known);
    let hmag_mean = hmag_sum / (used_bins.len() as f32);

    let mut pilot_evm = 0.0f32;
    if !pilot_bins.is_empty() {
        let sym_len = cfg.nfft + cfg.ncp;
        let s0 = xsync_len + train_len;
        let s1 = s0 + sym_len;
        if s1 <= rbb_cfo.len() {
            let y = fft(&rbb_cfo[s0 + cfg.ncp..s0 + cfg.ncp + cfg.nfft]);
            let pref = known_pilot_symbols(pilot_bins.len(), 1);
            let xeq_used = equalize_symbol_with_pilots(&y, &used_bins, &pilot_bins, &pref, &hest);
            let mut pilot_eq = Vec::new();
            let mut pilot_ref = Vec::new();
            for (k, pbin) in pilot_bins.iter().enumerate() {
                if let Some(pos) = used_bins.iter().position(|b| b == pbin) {
                    pilot_eq.push(xeq_used[pos]);
                    pilot_ref.push(pref[k]);
                }
            }
            pilot_evm = rms_evm(&pilot_eq, &pilot_ref);
        }
    }

    2.0 * cp_score + 0.3 * hmag_mean - 2.0 * train_evm - 1.5 * pilot_evm
}

/// Applies the repeated-half CFO estimate.
///
/// Rationale:
/// The repeated preamble yields a phase increment per sample,
///
/// $$\Delta\phi = \frac{\arg\sum r^*[n]r[n+L]}{L}$$
///
/// which is removed as a complex derotation.
pub(crate) fn coarse_cfo_correct(rbb: &[Complex32], cfg: &OfdmConfig) -> Vec<Complex32> {
    let l = cfg.sync_half_len;
    if l == 0 || rbb.len() < 2 * l {
        return rbb.to_vec();
    }
    let mut p = Complex32::new(0.0, 0.0);
    for n in 0..l {
        p += rbb[n].conj() * rbb[n + l];
    }
    let ph_inc = p.arg() / (l as f32);
    rbb.iter()
        .enumerate()
        .map(|(n, &x)| x * Complex32::from_polar(1.0, -(n as f32) * ph_inc))
        .collect()
}

/// Applies an explicit CFO hypothesis in hertz.
///
/// Rationale:
/// This is the direct derotation
///
/// $$r'[n] = r[n] e^{-j 2\pi f_\Delta n / f_s}$$
///
/// used by the short CFO search around the coarse estimate.
pub(crate) fn apply_cfo_hz(rbb: &[Complex32], fs: f32, cfo_hz: f32) -> Vec<Complex32> {
    rbb.iter()
        .enumerate()
        .map(|(n, &x)| {
            let phase = -2.0 * std::f32::consts::PI * cfo_hz * (n as f32) / fs;
            x * Complex32::from_polar(1.0, phase)
        })
        .collect()
}

/// Estimates CFO from the repeated-half preamble in hertz.
///
/// Rationale:
/// The repeated-half phase increment is converted back to hertz as
///
/// $$\hat f_\Delta = \frac{\Delta\phi}{2\pi} f_s$$
///
/// using the active baseband rate.
pub(crate) fn estimate_coarse_cfo_hz(rbb: &[Complex32], cfg: &OfdmConfig) -> f32 {
    let l = cfg.sync_half_len;
    if l == 0 || rbb.len() < 2 * l {
        return 0.0;
    }
    let mut p = Complex32::new(0.0, 0.0);
    for n in 0..l {
        p += rbb[n].conj() * rbb[n + l];
    }
    let ph_inc = p.arg() / (l as f32);
    ph_inc * active_baseband_fs(cfg) / (2.0 * std::f32::consts::PI)
}

// vim: set ts=4 sw=4 et:
