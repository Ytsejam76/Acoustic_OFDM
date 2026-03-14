// Copyright (c) 2026 Elias S. G. Carotti

use rustfft::{num_complex::Complex32, FftPlanner};

use crate::config::{Modulation, OfdmConfig};
use crate::packet::{bits_to_bytes, build_packet_bytes, bytes_to_bits, modulation_from_id, parse_packet_bytes, split_payload, PacketInfo};

#[derive(Clone, Debug)]
pub struct EncodedPacketMeta {
    pub frag_index: usize,
    pub frag_count: usize,
    pub packet_start: usize,
    pub packet_len: usize,
    pub xbb: Vec<Complex32>,
}

#[derive(Clone, Debug)]
pub struct EncodedBurst {
    pub audio: Vec<f32>,
    pub packet_meta: Vec<EncodedPacketMeta>,
}

/// Encodes a full payload into one multi-packet OFDM burst.
///
/// Parameters:
/// - `payload`: full application payload bytes.
/// - `cfg`: modem configuration.
/// Returns:
/// - `EncodedBurst`: audio samples and per-packet metadata.
pub fn encode_payload(payload: &[u8], cfg: &OfdmConfig) -> EncodedBurst {
    let chunks = split_payload(payload, cfg.packet_payload_bytes);
    let num_packets = chunks.len();
    let mut audio = Vec::<f32>::new();
    let mut meta = Vec::<EncodedPacketMeta>::with_capacity(num_packets);

    let pre_sil = (0.015 * cfg.fs) as usize;
    let post_sil = (0.020 * cfg.fs) as usize;

    for (i, chunk) in chunks.iter().enumerate() {
        let pkt_bytes = build_packet_bytes(chunk, i as u8, num_packets as u8, cfg);
        let xbb = tx_one_packet_baseband(&pkt_bytes, cfg);
        let pkt_audio = tx_one_packet(&pkt_bytes, cfg);

        audio.extend(std::iter::repeat(0.0).take(pre_sil));
        let packet_start = audio.len();
        audio.extend_from_slice(&pkt_audio);
        let packet_len = pkt_audio.len();
        audio.extend(std::iter::repeat(0.0).take(post_sil));

        meta.push(EncodedPacketMeta {
            frag_index: i,
            frag_count: num_packets,
            packet_start,
            packet_len,
            xbb,
        });
    }

    normalize_in_place(&mut audio, 0.85);
    EncodedBurst { audio, packet_meta: meta }
}

/// Decodes an encoded burst using oracle packet boundaries and baseband symbols.
///
/// Parameters:
/// - `burst`: encoded burst with per-packet metadata.
/// - `cfg`: modem configuration.
/// Returns:
/// - `Option<Vec<u8>>`: reconstructed payload on success, otherwise `None`.
pub fn decode_encoded_burst_oracle(burst: &EncodedBurst, cfg: &OfdmConfig) -> Option<Vec<u8>> {
    let mut frags: Vec<Option<Vec<u8>>> = vec![None; burst.packet_meta.len()];
    for m in &burst.packet_meta {
        let info = decode_packet_info_baseband(&m.xbb, cfg)?;
        if modulation_from_id(info.mod_id)? != cfg.modulation {
            return None;
        }
        let idx = info.frag_index as usize;
        if idx >= frags.len() || info.frag_count as usize != frags.len() {
            return None;
        }
        frags[idx] = Some(info.payload);
    }

    let mut out = Vec::new();
    for f in frags {
        out.extend(f?);
    }
    Some(out)
}

/// Encodes one payload fragment (single packet) into passband samples.
///
/// Parameters:
/// - `payload`: single-packet payload bytes (must fit packet payload budget).
/// - `cfg`: modem configuration.
/// Returns:
/// - `Vec<f32>`: passband waveform for one packet.
pub fn encode_single_packet_passband(payload: &[u8], cfg: &OfdmConfig) -> Vec<f32> {
    let pkt_bytes = build_packet_bytes(payload, 0, 1, cfg);
    tx_one_packet(&pkt_bytes, cfg)
}

/// Decodes one passband packet waveform into payload bytes.
///
/// Parameters:
/// - `pkt_audio`: passband packet waveform (wake + guard + OFDM body).
/// - `cfg`: modem configuration.
/// Returns:
/// - `Option<Vec<u8>>`: decoded payload bytes, or `None` on failure.
pub fn decode_single_packet_passband(pkt_audio: &[f32], cfg: &OfdmConfig) -> Option<Vec<u8>> {
    decode_packet_from_passband(pkt_audio, cfg).map(|p| p.payload)
}

/// Encodes one packet into passband samples (wake + guard + OFDM body).
///
/// Parameters:
/// - `pkt_bytes`: serialized packet bytes.
/// - `cfg`: modem configuration.
/// Returns:
/// - `Vec<f32>`: real passband waveform for one packet.
fn tx_one_packet(pkt_bytes: &[u8], cfg: &OfdmConfig) -> Vec<f32> {
    let xbb = tx_one_packet_baseband(pkt_bytes, cfg);
    let passband = iq_upconvert(&xbb, cfg.fs, cfg.fc);
    let wake = make_wake_tone(cfg);
    let guard_len = (cfg.wake_guard_ms * 1e-3 * cfg.fs) as usize;

    let mut out = Vec::with_capacity(wake.len() + guard_len + passband.len());
    out.extend_from_slice(&wake);
    out.extend(std::iter::repeat(0.0).take(guard_len));
    out.extend(passband);
    normalize_in_place(&mut out, 0.85);
    out
}

/// Encodes one packet into complex baseband OFDM symbols.
///
/// Parameters:
/// - `pkt_bytes`: serialized packet bytes.
/// - `cfg`: modem configuration.
/// Returns:
/// - `Vec<Complex32>`: baseband complex samples.
fn tx_one_packet_baseband(pkt_bytes: &[u8], cfg: &OfdmConfig) -> Vec<Complex32> {
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
    let sync_half = known_sync_half(cfg.sync_half_len);
    xbb.extend_from_slice(&sync_half);
    xbb.extend_from_slice(&sync_half);

    let train_known = known_training_symbols(used_bins.len(), cfg.modulation);
    let mut xtrain = vec![Complex32::new(0.0, 0.0); cfg.nfft];
    for (k, &bin) in used_bins.iter().enumerate() {
        xtrain[bin] = train_known[k];
    }
    let train_time = ifft(&xtrain);
    append_cp_symbol(&mut xbb, &train_time, cfg.ncp);

    for i in 0..n_data {
        let mut x = vec![Complex32::new(0.0, 0.0); cfg.nfft];
        for (k, &bin) in data_bins.iter().enumerate() {
            x[bin] = payload_syms[i * syms_per_ofdm + k];
        }
        if !pilot_bins.is_empty() {
            let pref = known_pilot_symbols(pilot_bins.len(), i + 1);
            for (k, &bin) in pilot_bins.iter().enumerate() {
                x[bin] = pref[k];
            }
        }
        let xt = ifft(&x);
        append_cp_symbol(&mut xbb, &xt, cfg.ncp);
    }

    xbb
}

/// Decodes one baseband packet and returns payload bytes.
///
/// Parameters:
/// - `rbb`: baseband complex packet samples.
/// - `cfg`: modem configuration.
/// Returns:
/// - `Option<Vec<u8>>`: decoded payload bytes, or `None` on failure.
pub fn decode_packet_baseband(rbb: &[Complex32], cfg: &OfdmConfig) -> Option<Vec<u8>> {
    decode_packet_info_baseband(rbb, cfg).map(|p| p.payload)
}

/// Decodes one baseband packet and returns full packet metadata.
///
/// Parameters:
/// - `rbb`: baseband complex packet samples.
/// - `cfg`: modem configuration.
/// Returns:
/// - `Option<PacketInfo>`: parsed packet metadata, or `None` on failure.
fn decode_packet_info_baseband(rbb: &[Complex32], cfg: &OfdmConfig) -> Option<PacketInfo> {
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
    let mut hest = vec![Complex32::new(1.0, 0.0); used_bins.len()];
    for (k, &bin) in used_bins.iter().enumerate() {
        hest[k] = ytrain[bin] / train_known[k];
    }

    let data_start = xsync_len + train_len;
    let sym_len = cfg.nfft + cfg.ncp;
    let max_payload_bytes = cfg.packet_payload_bytes + 16;
    let max_bits = max_payload_bytes * 8;
    let max_data_ofdm = max_bits.div_ceil(n_data_carriers * cfg.modulation.bits_per_symbol()) + 2;
    let mut rx_syms = Vec::<Complex32>::new();

    for i in 0..max_data_ofdm {
        let s0 = data_start + i * sym_len;
        let s1 = s0 + sym_len;
        if s1 > rbb.len() {
            break;
        }
        let y = fft(&rbb[s0 + cfg.ncp..s0 + cfg.ncp + cfg.nfft]);
        let mut xeq_used = Vec::with_capacity(used_bins.len());
        for (k, &bin) in used_bins.iter().enumerate() {
            xeq_used.push(y[bin] / hest[k]);
        }
        if !pilot_bins.is_empty() {
            let pref = known_pilot_symbols(pilot_bins.len(), i + 1);
            let mut num = Complex32::new(0.0, 0.0);
            let mut den = 0.0f32;
            for (k, pbin) in pilot_bins.iter().enumerate() {
                if let Some(pos) = used_bins.iter().position(|b| b == pbin) {
                    num += xeq_used[pos] * pref[k].conj();
                    den += pref[k].norm_sqr();
                }
            }
            if den > 0.0 {
                let g = num / den;
                if g.norm() > 1e-9 {
                    for v in &mut xeq_used {
                        *v /= g;
                    }
                }
            }
        }
        for dbin in &data_bins {
            if let Some(pos) = used_bins.iter().position(|b| b == dbin) {
                rx_syms.push(xeq_used[pos]);
            }
        }

        let bits = demap_bits(&rx_syms, cfg.modulation);
        let bytes = bits_to_bytes(&bits);
        if let Some((pkt, _)) = parse_packet_bytes(&bytes) {
            return Some(pkt);
        }
    }

    None
}

/// Decodes one passband packet (wake/guard prefixed) into packet metadata.
///
/// Parameters:
/// - `pkt_audio`: real passband packet waveform (wake + guard + body).
/// - `cfg`: modem configuration.
/// Returns:
/// - `Option<PacketInfo>`: parsed packet metadata, or `None` on failure.
fn decode_packet_from_passband(pkt_audio: &[f32], cfg: &OfdmConfig) -> Option<PacketInfo> {
    let wake_len = (cfg.wake_ms * 1e-3 * cfg.fs) as usize;
    let guard_len = (cfg.wake_guard_ms * 1e-3 * cfg.fs) as usize;
    if pkt_audio.len() <= wake_len + guard_len + cfg.nfft + cfg.ncp {
        return None;
    }

    let passband = &pkt_audio[wake_len + guard_len..];
    let pre = 128usize.min(passband.len());
    let mut chunk = vec![0.0f32; pre];
    chunk.extend_from_slice(passband);
    let rbb_full = iq_downconvert(&chunk, cfg.fs, cfg.fc);
    let rbb = &rbb_full[pre..];
    decode_packet_info_baseband(rbb, cfg)
}

/// Maps bits to complex symbols for the selected modulation.
///
/// Parameters:
/// - `bits`: input bits (`0/1` values).
/// - `modulation`: modulation scheme.
/// Returns:
/// - `Vec<Complex32>`: constellation symbols.
fn map_bits(bits: &[u8], modulation: Modulation) -> Vec<Complex32> {
    match modulation {
        Modulation::Bpsk => bits
            .iter()
            .map(|&b| if b == 0 { Complex32::new(1.0, 0.0) } else { Complex32::new(-1.0, 0.0) })
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

/// Demaps complex symbols back to hard bits for the selected modulation.
///
/// Parameters:
/// - `syms`: input constellation symbols.
/// - `modulation`: modulation scheme.
/// Returns:
/// - `Vec<u8>`: hard-decoded bits (`0/1` values).
fn demap_bits(syms: &[Complex32], modulation: Modulation) -> Vec<u8> {
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

/// Hard-slices symbols to nearest ideal constellation point.
///
/// Parameters:
/// - `syms`: noisy constellation symbols.
/// - `modulation`: modulation scheme.
/// Returns:
/// - `Vec<Complex32>`: nearest ideal symbols.
fn hard_slice_symbols(syms: &[Complex32], modulation: Modulation) -> Vec<Complex32> {
    match modulation {
        Modulation::Bpsk => syms
            .iter()
            .map(|s| if s.re < 0.0 { Complex32::new(-1.0, 0.0) } else { Complex32::new(1.0, 0.0) })
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

/// Returns the known training sequence used on active carriers.
///
/// Parameters:
/// - `n`: number of active carriers.
/// - `modulation`: modulation scheme.
/// Returns:
/// - `Vec<Complex32>`: known training symbols.
fn known_training_symbols(n: usize, modulation: Modulation) -> Vec<Complex32> {
    match modulation {
        Modulation::Bpsk => (0..n)
            .map(|i| if i % 2 == 0 { Complex32::new(1.0, 0.0) } else { Complex32::new(-1.0, 0.0) })
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

/// Returns known pilot symbols for one OFDM data symbol index.
///
/// Parameters:
/// - `n`: number of pilot carriers.
/// - `sym_idx`: OFDM symbol index, starting at 1.
/// Returns:
/// - `Vec<Complex32>`: deterministic pilot symbols.
fn known_pilot_symbols(n: usize, sym_idx: usize) -> Vec<Complex32> {
    (0..n)
        .map(|k| {
            let m = ((sym_idx - 1) + k) % 4;
            Complex32::from_polar(1.0, std::f32::consts::FRAC_PI_2 * (m as f32))
        })
        .collect()
}

/// Splits active bins into used, pilot and data bins.
///
/// Parameters:
/// - `cfg`: modem configuration.
/// Returns:
/// - `(Vec<usize>, Vec<usize>, Vec<usize>)`: `(used_bins, pilot_bins, data_bins)`.
fn ofdm_bin_plan(cfg: &OfdmConfig) -> (Vec<usize>, Vec<usize>, Vec<usize>) {
    let used_bins = cfg.used_bins.clone();
    let pilots_on = cfg.use_pilots.unwrap_or(matches!(cfg.modulation, Modulation::Qpsk));
    let pilot_bins = if pilots_on {
        let mut pilots = if cfg.pilot_bins.is_empty() {
            used_bins.clone()
        } else {
            cfg.pilot_bins
                .iter()
                .copied()
                .filter(|b| used_bins.contains(b))
                .collect::<Vec<_>>()
        };
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

/// Returns one half of the repeated sync preamble.
///
/// Parameters:
/// - `l`: half-preamble length in samples.
/// Returns:
/// - `Vec<Complex32>`: sync half sequence.
fn known_sync_half(l: usize) -> Vec<Complex32> {
    (0..l)
        .map(|idx| {
            let m = (3 * idx + 1) % 4;
            Complex32::from_polar(1.0, std::f32::consts::FRAC_PI_2 * (m as f32))
        })
        .collect()
}

/// Appends cyclic prefix and symbol samples to output.
///
/// Parameters:
/// - `out`: destination sample vector.
/// - `x`: one OFDM symbol (no CP).
/// - `ncp`: cyclic prefix length.
/// Returns:
/// - none.
fn append_cp_symbol(out: &mut Vec<Complex32>, x: &[Complex32], ncp: usize) {
    out.extend_from_slice(&x[x.len() - ncp..]);
    out.extend_from_slice(x);
}

/// Computes inverse FFT with unitary time-domain scaling.
///
/// Parameters:
/// - `x`: frequency-domain bins.
/// Returns:
/// - `Vec<Complex32>`: time-domain samples.
fn ifft(x: &[Complex32]) -> Vec<Complex32> {
    let mut planner = FftPlanner::<f32>::new();
    let fft = planner.plan_fft_inverse(x.len());
    let mut buf = x.to_vec();
    fft.process(&mut buf);
    let scale = 1.0 / (x.len() as f32);
    for v in &mut buf {
        *v *= scale;
    }
    buf
}

/// Computes forward FFT.
///
/// Parameters:
/// - `x`: time-domain samples.
/// Returns:
/// - `Vec<Complex32>`: frequency-domain bins.
fn fft(x: &[Complex32]) -> Vec<Complex32> {
    let mut planner = FftPlanner::<f32>::new();
    let fft = planner.plan_fft_forward(x.len());
    let mut buf = x.to_vec();
    fft.process(&mut buf);
    buf
}

/// IQ-upconverts complex baseband to real passband.
///
/// Parameters:
/// - `xbb`: baseband complex samples.
/// - `fs`: sample rate (Hz).
/// - `fc`: carrier frequency (Hz).
/// Returns:
/// - `Vec<f32>`: real passband samples.
fn iq_upconvert(xbb: &[Complex32], fs: f32, fc: f32) -> Vec<f32> {
    xbb.iter()
        .enumerate()
        .map(|(n, x)| {
            let ph = 2.0 * std::f32::consts::PI * fc * (n as f32) / fs;
            (x * Complex32::from_polar(1.0, ph)).re
        })
        .collect()
}

/// IQ-downconverts real passband to complex baseband and low-pass filters it.
///
/// Parameters:
/// - `y`: real passband samples.
/// - `fs`: sample rate (Hz).
/// - `fc`: carrier frequency (Hz).
/// Returns:
/// - `Vec<Complex32>`: filtered complex baseband samples.
fn iq_downconvert(y: &[f32], fs: f32, fc: f32) -> Vec<Complex32> {
    let mixed: Vec<Complex32> = y
        .iter()
        .enumerate()
        .map(|(n, &s)| {
            let ph = -2.0 * std::f32::consts::PI * fc * (n as f32) / fs;
            Complex32::new(s, 0.0) * Complex32::from_polar(1.0, ph)
        })
        .collect();
    lowpass_fir(&mixed, fs, 5_000.0, 65)
}

/// Builds wake preamble waveform (chirp or tone with ramp envelope).
///
/// Parameters:
/// - `cfg`: modem configuration.
/// Returns:
/// - `Vec<f32>`: wake waveform samples.
fn make_wake_tone(cfg: &OfdmConfig) -> Vec<f32> {
    let n = (cfg.wake_ms * 1e-3 * cfg.fs) as usize;
    let ramp = ((0.001 * cfg.fs) as usize).min(n / 4);
    let mut out = vec![0.0f32; n];
    for i in 0..n {
        let t = i as f32 / cfg.fs;
        let w = if cfg.use_chirp_sync {
            let tmax = ((n - 1) as f32 / cfg.fs).max(1.0 / cfg.fs);
            let k = (cfg.sync_chirp_f1 - cfg.sync_chirp_f0) / tmax;
            (2.0 * std::f32::consts::PI * (cfg.sync_chirp_f0 * t + 0.5 * k * t * t)).sin()
        } else {
            (2.0 * std::f32::consts::PI * cfg.wake_freq * t).sin()
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

/// Scales samples so max absolute value equals `target_peak`.
///
/// Parameters:
/// - `x`: samples modified in place.
/// - `target_peak`: desired max absolute amplitude.
/// Returns:
/// - none.
fn normalize_in_place(x: &mut [f32], target_peak: f32) {
    let peak = x
        .iter()
        .fold(0.0f32, |a, &b| if b.abs() > a { b.abs() } else { a });
    if peak > 0.0 {
        let g = target_peak / peak;
        for v in x {
            *v *= g;
        }
    }
}

/// Applies FIR low-pass filtering to complex samples.
///
/// Parameters:
/// - `x`: input complex samples.
/// - `fs`: sample rate (Hz).
/// - `cutoff_hz`: cutoff frequency (Hz).
/// - `len`: FIR length (taps).
/// Returns:
/// - `Vec<Complex32>`: filtered samples.
fn lowpass_fir(x: &[Complex32], fs: f32, cutoff_hz: f32, len: usize) -> Vec<Complex32> {
    let h = fir_lowpass(len, cutoff_hz, fs);
    let mut y = vec![Complex32::new(0.0, 0.0); x.len()];
    for n in 0..x.len() {
        let mut acc = Complex32::new(0.0, 0.0);
        let kmax = (n + 1).min(h.len());
        for k in 0..kmax {
            acc += x[n - k] * h[k];
        }
        y[n] = acc;
    }
    // Group delay compensation (match Octave behavior).
    let gd = (len - 1) / 2;
    let mut out = Vec::with_capacity(y.len());
    out.extend_from_slice(&y[gd..]);
    out.extend(std::iter::repeat(Complex32::new(0.0, 0.0)).take(gd));
    out
}

/// Designs a windowed-sinc low-pass FIR.
///
/// Parameters:
/// - `len`: FIR length (taps).
/// - `cutoff_hz`: cutoff frequency (Hz).
/// - `fs`: sample rate (Hz).
/// Returns:
/// - `Vec<f32>`: FIR coefficients.
fn fir_lowpass(len: usize, cutoff_hz: f32, fs: f32) -> Vec<f32> {
    let m = (len - 1) as f32 / 2.0;
    let mut h = vec![0.0f32; len];
    for (i, hi) in h.iter_mut().enumerate() {
        let n = i as f32 - m;
        let sinc = if n.abs() < 1e-6 {
            2.0 * cutoff_hz / fs
        } else {
            let x = 2.0 * std::f32::consts::PI * cutoff_hz * n / fs;
            (2.0 * cutoff_hz / fs) * (x.sin() / x)
        };
        let w = 0.54 - 0.46 * (2.0 * std::f32::consts::PI * i as f32 / (len as f32 - 1.0)).cos();
        *hi = sinc * w;
    }
    let sum: f32 = h.iter().sum();
    if sum != 0.0 {
        for v in &mut h {
            *v /= sum;
        }
    }
    h
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    /// Verifies full BPSK oracle burst round-trip.
    fn round_trip_oracle_burst() {
        let cfg = OfdmConfig::default();
        let payload = (0..80u8).collect::<Vec<_>>();
        let tx = encode_payload(&payload, &cfg);
        let rx = decode_encoded_burst_oracle(&tx, &cfg).expect("decode failed");
        assert_eq!(rx, payload);
    }

    #[test]
    /// Verifies full QPSK oracle burst round-trip.
    fn round_trip_oracle_burst_qpsk() {
        let mut cfg = OfdmConfig::default();
        cfg.modulation = Modulation::Qpsk;
        let payload = (0..96u8).collect::<Vec<_>>();
        let tx = encode_payload(&payload, &cfg);
        let rx = decode_encoded_burst_oracle(&tx, &cfg).expect("decode failed");
        assert_eq!(rx, payload);
    }

    #[test]
    /// Verifies QPSK oracle burst round-trip with pilots explicitly disabled.
    fn round_trip_oracle_burst_qpsk_no_pilots() {
        let mut cfg = OfdmConfig::default();
        cfg.modulation = Modulation::Qpsk;
        cfg.use_pilots = Some(false);
        let payload = (0..80u8).collect::<Vec<_>>();
        let tx = encode_payload(&payload, &cfg);
        let rx = decode_encoded_burst_oracle(&tx, &cfg).expect("decode failed");
        assert_eq!(rx, payload);
    }

    #[test]
    /// Verifies baseband packet encode/decode path.
    fn decode_single_packet_baseband() {
        let cfg = OfdmConfig::default();
        let payload: Vec<u8> = (100..120).collect();
        let pkt = build_packet_bytes(&payload, 0, 1, &cfg);
        let xbb = tx_one_packet_baseband(&pkt, &cfg);
        let out = decode_packet_baseband(&xbb, &cfg).expect("baseband decode failed");
        assert_eq!(out, payload);
    }

    #[test]
    /// Verifies passband packet encode/decode path.
    fn decode_single_packet_passband_test() {
        let cfg = OfdmConfig::default();
        let payload: Vec<u8> = (40..64).collect();
        let y = encode_single_packet_passband(&payload, &cfg);
        let out = decode_single_packet_passband(&y, &cfg).expect("passband decode failed");
        assert_eq!(out, payload);
    }

    #[test]
    /// Verifies BPSK/QPSK hard slicer output points.
    fn hard_slice_bpsk_and_qpsk() {
        let bpsk = vec![Complex32::new(0.3, 0.0), Complex32::new(-0.2, 1.0)];
        let bpsk_s = hard_slice_symbols(&bpsk, Modulation::Bpsk);
        assert_eq!(bpsk_s[0], Complex32::new(1.0, 0.0));
        assert_eq!(bpsk_s[1], Complex32::new(-1.0, 0.0));

        let qpsk = vec![
            Complex32::new(0.1, 0.2),
            Complex32::new(0.1, -0.2),
            Complex32::new(-0.1, 0.2),
            Complex32::new(-0.1, -0.2),
        ];
        let qpsk_s = hard_slice_symbols(&qpsk, Modulation::Qpsk);
        let k = 1.0 / 2.0f32.sqrt();
        assert_eq!(qpsk_s[0], Complex32::new(1.0 * k, 1.0 * k));
        assert_eq!(qpsk_s[1], Complex32::new(1.0 * k, -1.0 * k));
        assert_eq!(qpsk_s[2], Complex32::new(-1.0 * k, 1.0 * k));
        assert_eq!(qpsk_s[3], Complex32::new(-1.0 * k, -1.0 * k));
    }

    #[test]
    /// Verifies pilot count is capped by `num_pilots`.
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
}

// vim: set ts=4 sw=4 et:
