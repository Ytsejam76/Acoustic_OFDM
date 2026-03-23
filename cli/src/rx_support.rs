// Copyright (c) 2026 Elias S. G. Carotti

use acoustic_ofdm::{
    decode_single_packet_passband, diagnose_passband_window, OfdmConfig, PassbandDiagnostics,
    WakePreamble,
};

use crate::audio::{fir_bandpass, fir_filter};
use crate::debug_line;

pub(crate) fn signal_diag(x: &[f32]) -> (f32, f32, usize) {
    if x.is_empty() {
        return (0.0, 0.0, 0);
    }
    let mut peak = 0.0f32;
    let mut pwr = 0.0f32;
    let mut first = x.len();
    for (i, &s) in x.iter().enumerate() {
        let a = s.abs();
        if a > peak {
            peak = a;
        }
        if first == x.len() && a > 0.02 {
            first = i;
        }
        pwr += s * s;
    }
    let rms = (pwr / (x.len() as f32)).sqrt();
    (rms, peak, first)
}

pub(crate) fn clipping_diag(x: &[f32]) -> (usize, f32) {
    if x.is_empty() {
        return (0, 0.0);
    }
    let clipped = x.iter().filter(|&&s| s.abs() >= 0.995).count();
    (clipped, (clipped as f32) / (x.len() as f32))
}

pub(crate) fn filter_for_sync_detection(
    x: &[f32],
    fs: f32,
    hp_hz: f32,
    lp_hz: f32,
    enabled: bool,
) -> Vec<f32> {
    if !enabled || x.is_empty() {
        return x.to_vec();
    }
    let hp = hp_hz.clamp(10.0, 0.45 * fs);
    let lp = lp_hz.clamp((hp + 10.0).min(0.49 * fs), 0.49 * fs);
    let h = fir_bandpass(129, hp, lp, fs);
    fir_filter(x, &h)
}

fn pn_sequence(n: usize) -> Vec<f32> {
    let mut state: u16 = 0x01FF;
    let mut out = Vec::with_capacity(n);
    for _ in 0..n {
        let bit = (state & 1) as u8;
        out.push(if bit == 0 { -1.0 } else { 1.0 });
        let fb = ((state >> 8) ^ (state >> 4)) & 1;
        state = (state >> 1) | (fb << 8);
    }
    out
}

fn gold_sequence(n: usize) -> Vec<f32> {
    let mut s1: u16 = 0x01FF;
    let mut s2: u16 = 0x0155;
    let mut out = Vec::with_capacity(n);
    for _ in 0..n {
        let b1 = (s1 & 1) as u8;
        let b2 = (s2 & 1) as u8;
        out.push(if (b1 ^ b2) == 0 { -1.0 } else { 1.0 });

        let fb1 = ((s1 >> 8) ^ (s1 >> 4)) & 1;
        let fb2 = ((s2 >> 8) ^ (s2 >> 7) ^ (s2 >> 4) ^ (s2 >> 1)) & 1;
        s1 = (s1 >> 1) | (fb1 << 8);
        s2 = (s2 >> 1) | (fb2 << 8);
    }
    out
}

pub(crate) fn make_wake_ref(cfg: &OfdmConfig) -> Vec<f32> {
    let n = (cfg.wake_ms * 1e-3 * cfg.fs) as usize;
    let ramp = ((0.001 * cfg.fs) as usize).min(n / 4);
    let pn = pn_sequence(n);
    let gold = gold_sequence(n);
    let mut out = vec![0.0f32; n];
    for i in 0..n {
        let t = i as f32 / cfg.fs;
        let w = match cfg.wake_preamble {
            WakePreamble::Tone => (2.0 * std::f32::consts::PI * cfg.wake_freq * t).sin(),
            WakePreamble::Chirp => {
                let tmax = ((n - 1) as f32 / cfg.fs).max(1.0 / cfg.fs);
                let k = (cfg.sync_chirp_f1 - cfg.sync_chirp_f0) / tmax;
                (2.0 * std::f32::consts::PI * (cfg.sync_chirp_f0 * t + 0.5 * k * t * t)).sin()
            }
            WakePreamble::Pn => {
                let chip = pn[i];
                chip * (2.0 * std::f32::consts::PI * cfg.wake_freq * t).sin()
            }
            WakePreamble::Gold => {
                let chip = gold[i];
                chip * (2.0 * std::f32::consts::PI * cfg.wake_freq * t).sin()
            }
        };
        let env = if ramp > 1 && i < ramp {
            i as f32 / ramp as f32
        } else if ramp > 1 && i >= n - ramp {
            (n - 1 - i) as f32 / ramp as f32
        } else {
            1.0
        };
        out[i] = 0.7 * w * env;
    }
    out
}

pub(crate) fn wake_candidates(
    rx: &[f32],
    wake: &[f32],
    step: usize,
    top_k: usize,
    packet_len: usize,
) -> Vec<(usize, f32)> {
    if rx.len() < wake.len() || wake.is_empty() || step == 0 || top_k == 0 {
        return Vec::new();
    }
    let n = rx.len();
    let m = wake.len();
    let mut pref = vec![0.0f32; n + 1];
    for (i, &x) in rx.iter().enumerate() {
        pref[i + 1] = pref[i] + x * x;
    }
    let w_energy = wake.iter().map(|v| v * v).sum::<f32>().max(1e-12);
    let mut scored: Vec<(usize, f32)> = Vec::new();
    let min_sep = (m / 2).max(1);
    let last = n - m;
    for i in (0..=last).step_by(step) {
        let e = (pref[i + m] - pref[i]).max(1e-12);
        let mut dot = 0.0f32;
        for k in 0..m {
            dot += rx[i + k] * wake[k];
        }
        let corr = dot.abs() / (e.sqrt() * w_energy.sqrt());

        let pkt_end = i.saturating_add(packet_len).min(n);
        let post_e = (pref[pkt_end] - pref[i + m]).max(1e-12);
        let post_len = pkt_end.saturating_sub(i + m).max(1);
        let post_rms = (post_e / (post_len as f32)).sqrt();
        let wake_rms = (e / (m as f32)).sqrt();
        let energy_ratio = (post_rms / wake_rms.max(1e-6)).clamp(0.0, 4.0);

        let time_bias = 1.0 - 0.15 * ((i as f32) / (n as f32));
        let score = corr * (0.35 + 0.65 * energy_ratio) * time_bias.max(0.5);
        scored.push((i, score));
    }
    scored.sort_by(|a, b| b.1.total_cmp(&a.1));
    let mut filtered: Vec<(usize, f32)> = Vec::new();
    for (idx, sc) in scored {
        if filtered.iter().any(|(j, _)| idx.abs_diff(*j) < min_sep) {
            continue;
        }
        filtered.push((idx, sc));
        if filtered.len() == top_k {
            break;
        }
    }
    filtered
}

fn sample_linear(x: &[f32], pos: f32) -> f32 {
    if x.is_empty() || pos < 0.0 {
        return 0.0;
    }
    let i0 = pos.floor() as usize;
    if i0 >= x.len() {
        return 0.0;
    }
    let i1 = (i0 + 1).min(x.len() - 1);
    let a = pos - (i0 as f32);
    x[i0] * (1.0 - a) + x[i1] * a
}

fn fractional_wake_score(rx: &[f32], wake: &[f32], start: f32) -> f32 {
    if wake.is_empty() {
        return 0.0;
    }
    let w_energy = wake.iter().map(|v| v * v).sum::<f32>().max(1e-12);
    let mut dot = 0.0f32;
    let mut e = 0.0f32;
    for (k, &wk) in wake.iter().enumerate() {
        let s = sample_linear(rx, start + (k as f32));
        dot += s * wk;
        e += s * s;
    }
    dot.abs() / (e.sqrt().max(1e-12) * w_energy.sqrt())
}

pub(crate) fn refine_wake_candidates_fractional(
    rx: &[f32],
    wake: &[f32],
    cands: &[(usize, f32)],
) -> Vec<(usize, f32)> {
    let mut refined = Vec::with_capacity(cands.len());
    for (idx, base_score) in cands {
        let mut best_pos = *idx as f32;
        let mut best_score = *base_score;
        for di in -2..=2 {
            for frac_q in 0..4 {
                let pos = (*idx as f32) + (di as f32) + 0.25 * (frac_q as f32);
                if pos < 0.0 {
                    continue;
                }
                let score = fractional_wake_score(rx, wake, pos);
                if score > best_score {
                    best_score = score;
                    best_pos = pos;
                }
            }
        }
        refined.push((best_pos.round().max(0.0) as usize, best_score));
    }
    refined.sort_by(|a, b| b.1.total_cmp(&a.1));
    let mut dedup = Vec::with_capacity(refined.len());
    for (idx, sc) in refined {
        if dedup.iter().any(|(j, _)| idx.abs_diff(*j) < 4) {
            continue;
        }
        dedup.push((idx, sc));
    }
    dedup
}

#[derive(Clone, Copy, Debug)]
pub(crate) struct ActiveRegion {
    pub(crate) start: usize,
    pub(crate) end: usize,
    pub(crate) mean_rms: f32,
    pub(crate) peak_rms: f32,
}

pub(crate) fn burst_active_regions(x: &[f32], fs: f32) -> Vec<ActiveRegion> {
    if x.is_empty() {
        return Vec::new();
    }
    let win = ((0.010 * fs).round() as usize).max(1);
    let hop = ((0.002 * fs).round() as usize).max(1);
    if x.len() < win {
        let rms = (x.iter().map(|v| v * v).sum::<f32>() / (x.len() as f32)).sqrt();
        return vec![ActiveRegion {
            start: 0,
            end: x.len(),
            mean_rms: rms,
            peak_rms: rms,
        }];
    }
    let mut sum = x[..win].iter().map(|v| v * v).sum::<f32>();
    let mut env = Vec::<(usize, f32)>::new();
    let mut start = 0usize;
    loop {
        env.push((start, (sum / (win as f32)).sqrt()));
        if start + hop + win > x.len() {
            break;
        }
        for k in 0..hop {
            sum += x[start + win + k] * x[start + win + k] - x[start + k] * x[start + k];
        }
        start += hop;
    }
    let mut vals = env.iter().map(|(_, e)| *e).collect::<Vec<_>>();
    vals.sort_by(|a, b| a.total_cmp(b));
    let noise = vals[vals.len() / 5].max(1e-4);
    let th = (2.5 * noise).max(noise + 0.015);
    let min_run = ((0.025 * fs).round() as usize).max(hop);
    let pre = ((0.020 * fs).round() as usize).max(1);
    let post = ((0.180 * fs).round() as usize).max(1);
    let merge_gap = ((0.250 * fs).round() as usize).max(1);
    let mut runs = Vec::<ActiveRegion>::new();
    let mut cur: Option<(usize, usize, f32, f32, usize)> = None;
    for (s, e) in env {
        let active = e >= th;
        match (cur, active) {
            (None, true) => cur = Some((s, s + win, e, e, 1)),
            (Some((a, _b, sum_rms, peak_rms, nframes)), true) => {
                cur = Some((a, s + win, sum_rms + e, peak_rms.max(e), nframes + 1));
            }
            (Some((a, b, sum_rms, peak_rms, nframes)), false) => {
                if b.saturating_sub(a) >= min_run {
                    runs.push(ActiveRegion {
                        start: a.saturating_sub(pre),
                        end: (b + post).min(x.len()),
                        mean_rms: sum_rms / (nframes as f32),
                        peak_rms,
                    });
                }
                cur = None;
            }
            (None, false) => {}
        }
    }
    if let Some((a, b, sum_rms, peak_rms, nframes)) = cur {
        if b.saturating_sub(a) >= min_run {
            runs.push(ActiveRegion {
                start: a.saturating_sub(pre),
                end: (b + post).min(x.len()),
                mean_rms: sum_rms / (nframes as f32),
                peak_rms,
            });
        }
    }
    let mut merged = Vec::<ActiveRegion>::new();
    for r in runs {
        if let Some(last) = merged.last_mut() {
            if r.start <= last.end.saturating_add(merge_gap) {
                let last_len = last.end.saturating_sub(last.start).max(1) as f32;
                let r_len = r.end.saturating_sub(r.start).max(1) as f32;
                last.end = last.end.max(r.end);
                last.mean_rms =
                    (last.mean_rms * last_len + r.mean_rms * r_len) / (last_len + r_len);
                last.peak_rms = last.peak_rms.max(r.peak_rms);
                continue;
            }
        }
        merged.push(r);
    }
    merged.sort_by(|a, b| {
        let sa = a.peak_rms * (0.5 + a.mean_rms) * ((a.end - a.start) as f32).sqrt();
        let sb = b.peak_rms * (0.5 + b.mean_rms) * ((b.end - b.start) as f32).sqrt();
        sb.total_cmp(&sa).then_with(|| a.start.cmp(&b.start))
    });
    merged
}

pub(crate) fn filter_candidates_by_regions(
    cands: &[(usize, f32)],
    regions: &[ActiveRegion],
) -> Vec<(usize, f32)> {
    if regions.is_empty() {
        return Vec::new();
    }
    let keep_regions = regions.len().min(3);
    cands
        .iter()
        .copied()
        .filter_map(|(idx, score)| {
            regions
                .iter()
                .take(keep_regions)
                .find(|r| idx >= r.start && idx < r.end)
                .map(|r| {
                    let region_boost = (1.0 + 2.0 * r.mean_rms + 1.5 * r.peak_rms).max(1.0);
                    (idx, score * region_boost)
                })
        })
        .collect::<Vec<_>>()
        .tap_mut(|v| v.sort_by(|a, b| b.1.total_cmp(&a.1)))
}

pub(crate) fn candidate_region_index(
    idx: usize,
    regions: &[ActiveRegion],
    keep_regions: usize,
) -> Option<usize> {
    regions
        .iter()
        .take(keep_regions)
        .position(|r| idx >= r.start && idx < r.end)
}

trait TapMut: Sized {
    fn tap_mut<F: FnOnce(&mut Self)>(mut self, f: F) -> Self {
        f(&mut self);
        self
    }
}

impl<T> TapMut for T {}

pub(crate) fn estimated_packet_len_samples(cfg: &OfdmConfig, packet_bytes: Option<usize>) -> usize {
    let bps = cfg.modulation.bits_per_symbol();
    let used_bins = cfg.used_bins.len();
    let pilots_on = cfg
        .use_pilots
        .unwrap_or(matches!(cfg.modulation, acoustic_ofdm::Modulation::Qpsk));
    let pilot_count = if pilots_on {
        cfg.num_pilots
            .unwrap_or(cfg.pilot_bins.len())
            .min(cfg.pilot_bins.len())
            .min(used_bins)
    } else {
        0
    };
    let n_data_carriers = used_bins.saturating_sub(pilot_count).max(1);
    let serialized_packet_bytes = packet_bytes.unwrap_or(cfg.packet_payload_bytes + 11);
    let max_bits = serialized_packet_bytes * 8;
    let bits_per_ofdm = n_data_carriers * bps;
    let n_data_ofdm = max_bits.div_ceil(bits_per_ofdm) + 2;
    let retrain_count = cfg
        .retrain_interval_data_symbols
        .filter(|&interval| interval > 0)
        .map(|interval| n_data_ofdm / interval)
        .unwrap_or(0);
    let terminal_training = usize::from(cfg.terminal_training_symbol);
    let n_training_symbols = 1 + retrain_count + terminal_training;
    let baseband_len =
        2 * cfg.sync_half_len + (n_training_symbols + n_data_ofdm) * (cfg.nfft + cfg.ncp);
    let wake_len = (cfg.wake_ms * 1e-3 * cfg.fs) as usize;
    let guard_len = (cfg.wake_guard_ms * 1e-3 * cfg.fs) as usize;
    wake_len + guard_len + baseband_len
}

pub(crate) fn ranked_offset_hypotheses(
    rx: &[f32],
    seeds: &[(usize, f32, f32, bool)],
    est_pkt: usize,
    pad: usize,
    cfg: &OfdmConfig,
) -> Vec<(usize, PassbandDiagnostics)> {
    let back = ((cfg.fs * 0.008).round() as isize).max(1);
    let fwd = ((cfg.fs * 0.020).round() as isize).max(1);
    let step = ((cfg.fs * 0.0005).round() as isize).max(1);
    let mut scored = Vec::new();
    for (idx, _, _, _) in seeds.iter().take(3) {
        for dj in (-back..=fwd).step_by(step as usize) {
            let off_i = *idx as isize + dj;
            if off_i < 0 {
                continue;
            }
            let off = off_i as usize;
            if off >= rx.len() {
                continue;
            }
            let end = off.saturating_add(est_pkt + pad).min(rx.len());
            if end <= off + cfg.nfft + cfg.ncp {
                continue;
            }
            let diag = diagnose_passband_window(&rx[off..end], cfg);
            let score = diagnostic_candidate_score(&diag);
            scored.push((off, score, diag));
        }
    }
    scored.sort_by(|a, b| {
        let pa = plausible_candidate(&a.2);
        let pb = plausible_candidate(&b.2);
        pb.cmp(&pa)
            .then_with(|| b.1.total_cmp(&a.1))
            .then_with(|| a.0.cmp(&b.0))
    });
    let mut dedup = Vec::new();
    let min_sep = ((cfg.fs * 0.0015).round() as usize).max(1);
    for (off, _score, diag) in scored {
        if dedup.iter().any(|(j, _)| off.abs_diff(*j) < min_sep) {
            continue;
        }
        dedup.push((off, diag));
        if dedup.len() >= 12 {
            break;
        }
    }
    dedup
}

pub(crate) fn quick_realtime_decode(rx_raw: &[f32], rx_sync: &[f32], cfg: &OfdmConfig) -> Option<Vec<u8>> {
    let est_pkt = estimated_packet_len_samples(cfg, None);
    let pad = ((0.050 * cfg.fs).round() as usize).max(1);
    if rx_raw.len() < est_pkt + pad || rx_sync.len() < est_pkt + pad {
        return None;
    }
    let rt_span = ((cfg.fs * 2.5).round() as usize).max(est_pkt + pad);
    let rt_start = rx_raw.len().saturating_sub(rt_span);
    let rx_raw = &rx_raw[rt_start..];
    let rx_sync = &rx_sync[rt_start..];
    let wake = make_wake_ref(cfg);
    let step = ((cfg.fs * 0.002).round() as usize).max(1);
    let cands0 = wake_candidates(rx_sync, &wake, step, 2, est_pkt);
    let cands = refine_wake_candidates_fractional(rx_sync, &wake, &cands0);
    if let Some((seed, _)) = cands.first() {
        let local_step = ((cfg.fs * 0.0015).round() as isize).max(1);
        let back = ((cfg.fs * 0.006).round() as isize).max(1);
        let fwd = ((cfg.fs * 0.012).round() as isize).max(1);
        for dj in (-back..=fwd).step_by(local_step as usize) {
            let off_i = (*seed as isize) + dj;
            if off_i < 0 {
                continue;
            }
            let off = off_i as usize;
            let end = off.saturating_add(est_pkt + pad).min(rx_raw.len());
            if end <= off + cfg.nfft + cfg.ncp {
                continue;
            }
            if let Some(bytes) = decode_single_packet_passband(&rx_raw[off..end], cfg) {
                return Some(bytes);
            }
        }
    }
    None
}

pub(crate) fn print_passband_diagnostics(label: &str, pkt_audio: &[f32], cfg: &OfdmConfig) {
    let d = diagnose_passband_window(pkt_audio, cfg);
    debug_line!(
        "{label}: enough={} sync_off={} cfo={:.1}Hz train_rms={:.4} hest[min/mean/max]=[{:.3}/{:.3}/{:.3}] evm[train/pilot_pre/pilot_post/data]=[{:.3}/{:.3}/{:.3}/{:.3}] decoded={}",
        d.enough_samples,
        d.sync_off,
        d.cfo_hz,
        d.train_rms,
        d.hest_mag_min,
        d.hest_mag_mean,
        d.hest_mag_max,
        d.train_recon_evm,
        d.pilot_residual_evm,
        d.pilot_post_evm,
        d.post_eq_evm,
        d.decoded
    );
}

pub(crate) fn diagnostic_candidate_score(diag: &PassbandDiagnostics) -> f32 {
    if !diag.enough_samples {
        return -1e9;
    }
    let sync_penalty = 0.0001 * (diag.sync_off as f32);
    let cfo_penalty = 0.02 * diag.cfo_hz.abs();
    let train_bonus = 120.0 * diag.train_rms;
    let hest_bonus = 5.0 * diag.hest_mag_mean - 0.6 * (diag.hest_mag_max - diag.hest_mag_min);
    let evm_penalty =
        4.0 * diag.train_recon_evm + 6.0 * diag.pilot_residual_evm + 3.0 * diag.post_eq_evm;
    let decoded_bonus = if diag.decoded { 1000.0 } else { 0.0 };
    decoded_bonus + train_bonus + hest_bonus - sync_penalty - cfo_penalty - evm_penalty
}

pub(crate) fn plausible_candidate(diag: &PassbandDiagnostics) -> bool {
    diag.enough_samples
        && diag.train_rms >= 0.01
        && diag.hest_mag_mean >= 0.15
        && diag.hest_mag_max >= 0.4
        && (diag.train_recon_evm == 0.0 || diag.train_recon_evm <= 0.75)
        && (diag.pilot_residual_evm == 0.0 || diag.pilot_residual_evm <= 1.00)
        && (diag.post_eq_evm == 0.0 || diag.post_eq_evm <= 1.00)
}
