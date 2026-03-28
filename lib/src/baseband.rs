// Copyright (c) 2026 Elias S. G. Carotti

use rustfft::{num_complex::Complex32, FftPlanner};

use crate::config::{EqualizerMode, Modulation, OfdmConfig, PassbandMode};
use crate::packet::{
    bits_to_bytes, build_packet_bytes, bytes_to_bits, parse_packet_bytes, PacketInfo,
};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum PacketSymbolKind {
    Training,
    Data,
}

pub fn encode_single_packet_baseband(payload: &[u8], cfg: &OfdmConfig) -> Vec<Complex32> {
    let pkt_bytes = build_packet_bytes(payload, 0, 1, cfg);
    tx_one_packet_baseband(&pkt_bytes, cfg)
}

pub fn decode_packet_baseband(rbb: &[Complex32], cfg: &OfdmConfig) -> Option<Vec<u8>> {
    decode_packet_info_baseband(rbb, cfg).map(|p| p.payload)
}

pub fn expected_single_packet_data_symbols(payload: &[u8], cfg: &OfdmConfig) -> Vec<Complex32> {
    let pkt_bytes = build_packet_bytes(payload, 0, 1, cfg);
    expected_packet_data_symbols(&pkt_bytes, cfg)
}

pub fn recover_single_packet_data_symbols(
    rbb: &[Complex32],
    cfg: &OfdmConfig,
) -> Option<Vec<Complex32>> {
    equalized_data_symbols_baseband(rbb, cfg)
}

pub(crate) fn tx_one_packet_baseband(pkt_bytes: &[u8], cfg: &OfdmConfig) -> Vec<Complex32> {
    let (used_bins, pilot_bins, data_bins) = ofdm_bin_plan(cfg);
    let bits = bytes_to_bits(pkt_bytes);
    let mut payload_syms = map_bits(&bits, cfg.modulation);
    let syms_per_ofdm = data_bins.len();
    if syms_per_ofdm == 0 {
        return Vec::new();
    }
    let n_data = payload_syms.len().div_ceil(syms_per_ofdm);
    payload_syms.resize(n_data * syms_per_ofdm, Complex32::new(0.0, 0.0));

    let mut xbb = Vec::<Complex32>::new();
    let sync_half = known_sync_half(cfg);
    xbb.extend_from_slice(&sync_half);
    xbb.extend_from_slice(&sync_half);

    let train_time = training_symbol_time_domain(&used_bins, cfg);
    append_cp_symbol(&mut xbb, &train_time, cfg.ncp);

    let symbol_plan = packet_symbol_plan(n_data, cfg);
    let mut data_idx = 0usize;
    let mut data_symbol_idx = 0usize;
    for kind in symbol_plan {
        if kind == PacketSymbolKind::Training {
            append_cp_symbol(&mut xbb, &train_time, cfg.ncp);
            continue;
        }
        let mut x = vec![Complex32::new(0.0, 0.0); cfg.nfft];
        for (k, &bin) in data_bins.iter().enumerate() {
            x[bin] = payload_syms[data_idx * syms_per_ofdm + k];
        }
        if !pilot_bins.is_empty() {
            let pref = known_pilot_symbols(pilot_bins.len(), data_symbol_idx + 1);
            for (k, &bin) in pilot_bins.iter().enumerate() {
                x[bin] = pref[k];
            }
        }
        let xt = ifft(&x);
        append_cp_symbol(&mut xbb, &xt, cfg.ncp);
        data_idx += 1;
        data_symbol_idx += 1;
    }

    xbb
}

pub(crate) fn decode_packet_info_baseband(rbb: &[Complex32], cfg: &OfdmConfig) -> Option<PacketInfo> {
    let rx_syms = equalized_data_symbols_baseband(rbb, cfg)?;
    recover_packet_from_symbols(&rx_syms, cfg)
}

fn equalized_data_symbols_baseband(rbb: &[Complex32], cfg: &OfdmConfig) -> Option<Vec<Complex32>> {
    let (used_bins, pilot_bins, data_bins) = ofdm_bin_plan(cfg);
    let n_data_carriers = data_bins.len();
    if n_data_carriers == 0 {
        return None;
    }
    let xsync_len = 2 * cfg.sync_half_len;
    let train_len = cfg.nfft + cfg.ncp;
    if rbb.len() < xsync_len + train_len {
        return None;
    }

    let train_start = xsync_len;
    let train_no_cp = &rbb[train_start + cfg.ncp..train_start + cfg.ncp + cfg.nfft];
    let ytrain = fft(train_no_cp);

    let train_known = known_training_symbols(used_bins.len(), cfg.modulation);
    let mut hest = equalizer_initial_channel(cfg, &ytrain, &used_bins, &train_known);

    let data_start = xsync_len + train_len;
    let max_payload_bytes = cfg.packet_payload_bytes + 16;
    let max_bits = max_payload_bytes * 8;
    let max_data_ofdm = max_bits.div_ceil(n_data_carriers * cfg.modulation.bits_per_symbol()) + 2;
    let symbol_plan = packet_symbol_plan(max_data_ofdm, cfg);
    let mut rx_syms = Vec::<Complex32>::new();
    let mut data_symbol_idx = 0usize;

    for (sym_idx, kind) in symbol_plan.into_iter().enumerate() {
        let s0 = data_start + sym_idx * (cfg.nfft + cfg.ncp);
        let s1 = s0 + cfg.nfft + cfg.ncp;
        if s1 > rbb.len() {
            break;
        }
        let y = fft(&rbb[s0 + cfg.ncp..s0 + cfg.ncp + cfg.nfft]);
        if kind == PacketSymbolKind::Training {
            equalizer_refresh_channel(cfg, &mut hest, &y, &used_bins, &train_known);
            continue;
        }
        let pref = known_pilot_symbols(pilot_bins.len(), data_symbol_idx + 1);
        let xeq_used = equalize_symbol_with_pilots(&y, &used_bins, &pilot_bins, &pref, &hest);
        for dbin in &data_bins {
            if let Some(pos) = used_bins.iter().position(|b| b == dbin) {
                rx_syms.push(xeq_used[pos]);
            }
        }
        data_symbol_idx += 1;
    }

    Some(rx_syms)
}

fn expected_packet_data_symbols(pkt_bytes: &[u8], cfg: &OfdmConfig) -> Vec<Complex32> {
    let (_, _, data_bins) = ofdm_bin_plan(cfg);
    let syms_per_ofdm = data_bins.len();
    if syms_per_ofdm == 0 {
        return Vec::new();
    }
    let bits = bytes_to_bits(pkt_bytes);
    let mut payload_syms = map_bits(&bits, cfg.modulation);
    let n_data = payload_syms.len().div_ceil(syms_per_ofdm);
    payload_syms.resize(n_data * syms_per_ofdm, Complex32::new(0.0, 0.0));
    payload_syms
}

pub(crate) fn packet_symbol_plan(n_data_symbols: usize, cfg: &OfdmConfig) -> Vec<PacketSymbolKind> {
    let mut plan = Vec::with_capacity(
        n_data_symbols
            + cfg
                .retrain_interval_data_symbols
                .map(|intv| if intv > 0 { n_data_symbols / intv } else { 0 })
                .unwrap_or(0)
            + usize::from(cfg.terminal_training_symbol),
    );
    for data_idx in 0..n_data_symbols {
        if let Some(interval) = cfg.retrain_interval_data_symbols {
            if interval > 0 && data_idx > 0 && data_idx % interval == 0 {
                plan.push(PacketSymbolKind::Training);
            }
        }
        plan.push(PacketSymbolKind::Data);
    }
    if cfg.terminal_training_symbol {
        plan.push(PacketSymbolKind::Training);
    }
    plan
}

pub(crate) fn training_symbol_time_domain(used_bins: &[usize], cfg: &OfdmConfig) -> Vec<Complex32> {
    let train_known = known_training_symbols(used_bins.len(), cfg.modulation);
    let mut xtrain = vec![Complex32::new(0.0, 0.0); cfg.nfft];
    for (k, &bin) in used_bins.iter().enumerate() {
        xtrain[bin] = train_known[k];
    }
    ifft(&xtrain)
}

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

pub(crate) fn equalizer_initial_channel(
    cfg: &OfdmConfig,
    ytrain: &[Complex32],
    used_bins: &[usize],
    train_known: &[Complex32],
) -> Vec<Complex32> {
    match cfg.equalizer_mode {
        EqualizerMode::TrainingPilot => estimate_channel_from_training(ytrain, used_bins, train_known),
        EqualizerMode::PilotOnly => vec![Complex32::new(1.0, 0.0); used_bins.len()],
    }
}

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
            for i in 1..pilot_phase_pts.len() {
                let mut phi = pilot_phase_pts[i].1;
                let prev = pilot_phase_pts[i - 1].1;
                while phi - prev > std::f32::consts::PI {
                    phi -= 2.0 * std::f32::consts::PI;
                }
                while phi - prev < -std::f32::consts::PI {
                    phi += 2.0 * std::f32::consts::PI;
                }
                pilot_phase_pts[i].1 = phi;
            }

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

pub(crate) fn map_bits(bits: &[u8], modulation: Modulation) -> Vec<Complex32> {
    match modulation {
        Modulation::Bpsk => bits
            .iter()
            .map(|&b| {
                if b == 0 {
                    Complex32::new(1.0, 0.0)
                } else {
                    Complex32::new(-1.0, 0.0)
                }
            })
            .collect(),
        Modulation::Qpsk => {
            let mut out = Vec::with_capacity(bits.len().div_ceil(2));
            let mut i = 0;
            while i < bits.len() {
                let b0 = bits[i];
                let b1 = if i + 1 < bits.len() { bits[i + 1] } else { 0 };
                let re = if b0 == 0 { 1.0 } else { -1.0 };
                let im = if b1 == 0 { 1.0 } else { -1.0 };
                out.push(Complex32::new(re, im) * (1.0 / 2.0f32.sqrt()));
                i += 2;
            }
            out
        }
    }
}

pub(crate) fn demap_bits(syms: &[Complex32], modulation: Modulation) -> Vec<u8> {
    match modulation {
        Modulation::Bpsk => syms.iter().map(|s| (s.re < 0.0) as u8).collect(),
        Modulation::Qpsk => {
            let mut bits = Vec::with_capacity(syms.len() * 2);
            for s in syms {
                bits.push((s.re < 0.0) as u8);
                bits.push((s.im < 0.0) as u8);
            }
            bits
        }
    }
}

fn recover_packet_from_symbols(syms: &[Complex32], cfg: &OfdmConfig) -> Option<PacketInfo> {
    let bits = demap_bits(syms, cfg.modulation);
    let bytes = bits_to_bytes(&bits);
    parse_packet_bytes(&bytes).map(|(pkt, _)| pkt)
}

pub(crate) fn known_training_symbols(n: usize, modulation: Modulation) -> Vec<Complex32> {
    match modulation {
        Modulation::Bpsk => (0..n)
            .map(|i| {
                if i % 2 == 0 {
                    Complex32::new(1.0, 0.0)
                } else {
                    Complex32::new(-1.0, 0.0)
                }
            })
            .collect(),
        Modulation::Qpsk => {
            let base = [
                Complex32::new(1.0, 1.0),
                Complex32::new(1.0, -1.0),
                Complex32::new(-1.0, 1.0),
                Complex32::new(-1.0, -1.0),
            ];
            (0..n)
                .map(|i| base[i % 4] * (1.0 / 2.0f32.sqrt()))
                .collect()
        }
    }
}

pub(crate) fn known_pilot_symbols(n: usize, sym_idx: usize) -> Vec<Complex32> {
    (0..n)
        .map(|k| {
            let m = ((sym_idx - 1) + k) % 4;
            Complex32::from_polar(1.0, std::f32::consts::FRAC_PI_2 * (m as f32))
        })
        .collect()
}

pub(crate) fn ofdm_bin_plan(cfg: &OfdmConfig) -> (Vec<usize>, Vec<usize>, Vec<usize>) {
    let (used_bins, pilot_candidates) = resolve_bins_with_base_freq(cfg);
    let pilots_on = cfg
        .use_pilots
        .unwrap_or(matches!(cfg.modulation, Modulation::Qpsk));
    let pilot_bins = if pilots_on {
        let mut pilots = pilot_candidates
            .iter()
            .copied()
            .filter(|b| used_bins.contains(b))
            .collect::<Vec<_>>();
        if let Some(n) = cfg.num_pilots {
            pilots.truncate(n.min(pilots.len()));
        }
        pilots
    } else {
        Vec::new()
    };
    let data_bins = used_bins
        .iter()
        .copied()
        .filter(|b| !pilot_bins.contains(b))
        .collect::<Vec<_>>();
    (used_bins, pilot_bins, data_bins)
}

fn dedup_stable(v: Vec<usize>) -> Vec<usize> {
    let mut out = Vec::with_capacity(v.len());
    for x in v {
        if !out.contains(&x) {
            out.push(x);
        }
    }
    out
}

fn resolve_bins_with_base_freq(cfg: &OfdmConfig) -> (Vec<usize>, Vec<usize>) {
    let mut used_bins = cfg.used_bins.clone();
    let mut pilot_candidates = if cfg.pilot_bins.is_empty() {
        used_bins.clone()
    } else {
        cfg.pilot_bins.clone()
    };

    if let Some(base_hz) = cfg.base_freq_hz {
        let numerology_fs = match cfg.passband_mode {
            PassbandMode::Legacy => cfg.fs,
            PassbandMode::Iq => cfg.fs_baseband,
        };
        let df = numerology_fs / (cfg.nfft as f32);
        if df > 0.0 && !used_bins.is_empty() {
            let target_bin = (base_hz / df).round() as isize;
            let current_min = *used_bins.iter().min().unwrap_or(&1) as isize;
            let shift = target_bin.max(1) - current_min;
            let kmax = (cfg.nfft / 2).saturating_sub(1) as isize;
            used_bins = used_bins
                .iter()
                .map(|&b| ((b as isize + shift).clamp(1, kmax)) as usize)
                .collect();
            pilot_candidates = pilot_candidates
                .iter()
                .map(|&b| ((b as isize + shift).clamp(1, kmax)) as usize)
                .collect();
        }
    }

    (dedup_stable(used_bins), dedup_stable(pilot_candidates))
}

pub(crate) fn known_sync_half(cfg: &OfdmConfig) -> Vec<Complex32> {
    let l = cfg.sync_half_len;
    if l == 0 {
        return Vec::new();
    }

    let (used_bins, _pilot_bins, _data_bins) = ofdm_bin_plan(cfg);
    if used_bins.is_empty() {
        return vec![Complex32::new(0.0, 0.0); l];
    }

    let mut out = Vec::with_capacity(l);
    for n in 0..l {
        let mut acc = Complex32::new(0.0, 0.0);
        for (i, &bin) in used_bins.iter().enumerate() {
            let phase0 = std::f32::consts::FRAC_PI_2 * (((3 * i + 1) % 4) as f32);
            let phase =
                phase0 + 2.0 * std::f32::consts::PI * (bin as f32) * (n as f32) / (cfg.nfft as f32);
            acc += Complex32::from_polar(1.0, phase);
        }
        out.push(acc);
    }

    let scale = 1.0 / (used_bins.len() as f32).sqrt();
    for v in &mut out {
        *v *= scale;
    }
    out
}

pub(crate) fn append_cp_symbol(out: &mut Vec<Complex32>, x: &[Complex32], ncp: usize) {
    out.extend_from_slice(&x[x.len() - ncp..]);
    out.extend_from_slice(x);
}

pub(crate) fn ifft(x: &[Complex32]) -> Vec<Complex32> {
    let mut planner = FftPlanner::<f32>::new();
    let fft = planner.plan_fft_inverse(x.len());
    let mut buf = x.to_vec();
    fft.process(&mut buf);
    let scale = 1.0 / (x.len() as f32).sqrt();
    for v in &mut buf {
        *v *= scale;
    }
    buf
}

pub(crate) fn fft(x: &[Complex32]) -> Vec<Complex32> {
    let mut planner = FftPlanner::<f32>::new();
    let fft = planner.plan_fft_forward(x.len());
    let mut buf = x.to_vec();
    fft.process(&mut buf);
    let scale = 1.0 / (x.len() as f32).sqrt();
    for v in &mut buf {
        *v *= scale;
    }
    buf
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn decode_single_packet_baseband() {
        let cfg = OfdmConfig::default();
        let payload: Vec<u8> = (100..120).collect();
        let xbb = encode_single_packet_baseband(&payload, &cfg);
        let out = decode_packet_baseband(&xbb, &cfg).expect("baseband decode failed");
        assert_eq!(out, payload);
    }

    #[test]
    fn recover_single_packet_bpsk_symbols() {
        let cfg = OfdmConfig::default();
        let payload: Vec<u8> = (0..24).collect();
        let xbb = encode_single_packet_baseband(&payload, &cfg);
        let got = recover_single_packet_data_symbols(&xbb, &cfg).expect("symbol recovery failed");
        let expect = expected_single_packet_data_symbols(&payload, &cfg);
        assert_eq!(got.len(), expect.len());
        for (g, e) in got.iter().zip(expect.iter()) {
            assert!((*g - *e).norm() < 0.05, "got={g:?} expect={e:?}");
        }
    }

    #[test]
    fn recover_single_packet_qpsk_symbols() {
        let mut cfg = OfdmConfig::default();
        cfg.modulation = Modulation::Qpsk;
        cfg.use_pilots = Some(true);
        let payload: Vec<u8> = (0..24).map(|x| x ^ 0x5a).collect();
        let xbb = encode_single_packet_baseband(&payload, &cfg);
        let got = recover_single_packet_data_symbols(&xbb, &cfg).expect("symbol recovery failed");
        let expect = expected_single_packet_data_symbols(&payload, &cfg);
        assert_eq!(got.len(), expect.len());
        for (g, e) in got.iter().zip(expect.iter()) {
            assert!((*g - *e).norm() < 0.08, "got={g:?} expect={e:?}");
        }
    }

    #[test]
    fn pilot_count_capped_by_num_pilots() {
        let mut cfg = OfdmConfig::default();
        cfg.modulation = Modulation::Qpsk;
        cfg.use_pilots = Some(true);
        cfg.used_bins = vec![2, 3, 4, 5];
        cfg.pilot_bins = vec![2, 4, 5];
        cfg.num_pilots = Some(2);
        let (_used, pilots, data) = ofdm_bin_plan(&cfg);
        assert_eq!(pilots, vec![2, 4]);
        assert_eq!(data, vec![3, 5]);
    }

    #[test]
    fn base_frequency_shifts_bins() {
        let mut cfg = OfdmConfig::default();
        cfg.base_freq_hz = Some(2_000.0);
        cfg.use_pilots = Some(false);
        let (used, pilots, data) = ofdm_bin_plan(&cfg);
        assert_eq!(used, vec![93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104]);
        assert!(pilots.is_empty());
        assert_eq!(data, vec![93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104]);
    }
}
