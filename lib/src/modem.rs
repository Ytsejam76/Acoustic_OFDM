// Copyright (c) 2026 Elias S. G. Carotti

use rustfft::num_complex::Complex32;

use crate::baseband::{
    decode_packet_info_baseband, fft, known_pilot_symbols, known_training_symbols, ofdm_bin_plan,
    packet_symbol_plan, recover_decided_packet_bytes_baseband, tx_one_packet_baseband,
    PacketSymbolKind,
};
use crate::config::{OfdmConfig, PassbandMode};
use crate::debug::{
    EncodedBurst, EncodedPacketMeta, PassbandBinDump, PassbandBinDumpRow,
    PassbandChannelCompareDump, PassbandChannelCompareRow, PassbandConstellationDump,
    PassbandDiagnostics, PassbandIqChainDump, PassbandPilotTrackDump, PassbandSyncDump,
};
use crate::equalizer::{
    decision_directed_evm, equalize_symbol_with_pilots, equalizer_initial_channel,
    equalizer_refresh_channel, equalizer_reset_tracking, regularized_equalize, rms_evm,
    EqualizerTrackingState,
};
use crate::packet::{
    build_packet_bytes, fec_encoded_bits_len, modulation_from_id, split_payload, PacketInfo,
};
use crate::sync::{
    active_baseband_fs, apply_cfo_hz, coarse_cfo_correct, estimate_coarse_cfo_hz,
    find_repeated_half_sync_offset, refine_sync_offset, repeated_half_sync_metrics,
    resample_from_offset, resample_from_offset_rate, sample_complex_linear,
};
use crate::wake::make_wake_tone;

fn pilot_phase_error(
    xeq_used: &[Complex32],
    used_bins: &[usize],
    pilot_bins: &[usize],
    pref: &[Complex32],
) -> Option<f32> {
    if pilot_bins.is_empty() || pref.is_empty() {
        return None;
    }
    let mut acc = Complex32::new(0.0, 0.0);
    for (k, pbin) in pilot_bins.iter().enumerate() {
        if let Some(pos) = used_bins.iter().position(|b| b == pbin) {
            acc += xeq_used[pos] * pref[k].conj();
        }
    }
    if acc.norm() <= 1.0e-9 {
        None
    } else {
        Some(acc.arg())
    }
}

fn resample_complex_linear_rate(x: &[Complex32], fs_in: f32, fs_out: f32) -> Vec<Complex32> {
    if x.is_empty() || fs_in <= 0.0 || fs_out <= 0.0 {
        return Vec::new();
    }
    if (fs_in - fs_out).abs() <= 1.0e-6 {
        return x.to_vec();
    }
    let out_len = ((x.len() as f32) * fs_out / fs_in).round().max(1.0) as usize;
    let step = fs_in / fs_out;
    let mut out = Vec::with_capacity(out_len);
    for n in 0..out_len {
        out.push(sample_complex_linear(x, n as f32 * step));
    }
    out
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
    EncodedBurst {
        audio,
        packet_meta: meta,
    }
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

/// Encodes one payload fragment into passband OFDM-body samples only.
///
/// Parameters:
/// - `payload`: single-packet payload bytes (must fit packet payload budget).
/// - `cfg`: modem configuration.
/// Returns:
/// - `Vec<f32>`: passband waveform for the OFDM body only, without wake/guard.
pub fn encode_single_packet_passband_body(payload: &[u8], cfg: &OfdmConfig) -> Vec<f32> {
    let pkt_bytes = build_packet_bytes(payload, 0, 1, cfg);
    tx_one_packet_body(&pkt_bytes, cfg)
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

/// Decodes one passband packet waveform into payload bytes using a known sync offset.
///
/// Parameters:
/// - `pkt_audio`: passband packet waveform (wake + guard + OFDM body).
/// - `cfg`: modem configuration.
/// - `sync_off`: known repeated-half/training alignment offset in baseband samples.
/// Returns:
/// - `Option<Vec<u8>>`: decoded payload bytes, or `None` on failure.
pub fn decode_single_packet_passband_with_sync(
    pkt_audio: &[f32],
    cfg: &OfdmConfig,
    sync_off: f32,
) -> Option<Vec<u8>> {
    decode_single_packet_passband_with_sync_rate(pkt_audio, cfg, sync_off, 1.0)
}

pub fn decode_single_packet_passband_with_sync_rate(
    pkt_audio: &[f32],
    cfg: &OfdmConfig,
    sync_off: f32,
    time_scale: f32,
) -> Option<Vec<u8>> {
    decode_packet_from_passband_with_sync_rate(pkt_audio, cfg, sync_off, time_scale)
        .map(|p| p.payload)
}

pub fn recover_decided_packet_bytes_passband_with_sync(
    pkt_audio: &[f32],
    cfg: &OfdmConfig,
    sync_off: f32,
) -> Option<Vec<u8>> {
    recover_decided_packet_bytes_passband_with_sync_rate(pkt_audio, cfg, sync_off, 1.0)
}

pub fn recover_decided_packet_bytes_passband_with_sync_rate(
    pkt_audio: &[f32],
    cfg: &OfdmConfig,
    sync_off: f32,
    time_scale: f32,
) -> Option<Vec<u8>> {
    let wake_len = (cfg.wake_ms * 1e-3 * cfg.fs) as usize;
    let guard_len = (cfg.wake_guard_ms * 1e-3 * cfg.fs) as usize;
    if pkt_audio.len() <= wake_len + guard_len + cfg.nfft + cfg.ncp {
        return None;
    }

    let passband = &pkt_audio[wake_len + guard_len..];
    let pre = 128usize.min(passband.len());
    let mut chunk = vec![0.0f32; pre];
    chunk.extend_from_slice(passband);
    let rbb_full = downconvert_passband(&chunk, cfg);
    let rbb = &rbb_full[pre..];
    let rbb_sync = resample_from_offset_rate(rbb, sync_off, time_scale);
    let rbb_cfo = coarse_cfo_correct(&rbb_sync, cfg);
    recover_decided_packet_bytes_baseband(&rbb_cfo, cfg)
}

/// Produces diagnostics for one passband packet window.
///
/// Parameters:
/// - `pkt_audio`: passband packet waveform window.
/// - `cfg`: modem configuration.
/// Returns:
/// - `PassbandDiagnostics`: sync/CFO/equalization diagnostics.
pub fn diagnose_passband_window(pkt_audio: &[f32], cfg: &OfdmConfig) -> PassbandDiagnostics {
    diagnose_passband_window_with_sync_opt(pkt_audio, cfg, None)
}

/// Produces diagnostics for one passband packet window using a known sync offset.
///
/// Parameters:
/// - `pkt_audio`: passband packet waveform window.
/// - `cfg`: modem configuration.
/// - `sync_off`: known repeated-half/training alignment offset in baseband samples.
/// Returns:
/// - `PassbandDiagnostics`: sync/CFO/equalization diagnostics.
pub fn diagnose_passband_window_with_sync(
    pkt_audio: &[f32],
    cfg: &OfdmConfig,
    sync_off: f32,
) -> PassbandDiagnostics {
    diagnose_passband_window_with_sync_rate(pkt_audio, cfg, sync_off, 1.0)
}

pub fn diagnose_passband_window_with_sync_rate(
    pkt_audio: &[f32],
    cfg: &OfdmConfig,
    sync_off: f32,
    time_scale: f32,
) -> PassbandDiagnostics {
    diagnose_passband_window_with_sync_opt(pkt_audio, cfg, Some((sync_off, time_scale)))
}

fn diagnose_passband_window_with_sync_opt(
    pkt_audio: &[f32],
    cfg: &OfdmConfig,
    forced_sync: Option<(f32, f32)>,
) -> PassbandDiagnostics {
    let wake_len = (cfg.wake_ms * 1e-3 * cfg.fs) as usize;
    let guard_len = (cfg.wake_guard_ms * 1e-3 * cfg.fs) as usize;
    if pkt_audio.len() <= wake_len + guard_len + cfg.nfft + cfg.ncp {
        return PassbandDiagnostics {
            enough_samples: false,
            sync_off: 0,
            cfo_hz: 0.0,
            sync_rms: 0.0,
            sync_peak: 0.0,
            post_rms: 0.0,
            post_peak: 0.0,
            train_rms: 0.0,
            hest_mag_min: 0.0,
            hest_mag_mean: 0.0,
            hest_mag_max: 0.0,
            train_recon_evm: 0.0,
            pilot_residual_evm: 0.0,
            pilot_post_evm: 0.0,
            post_eq_evm: 0.0,
            decoded: false,
            decoded_payload_len: None,
        };
    }

    let passband = &pkt_audio[wake_len + guard_len..];
    let pre = 128usize.min(passband.len());
    let mut chunk = vec![0.0f32; pre];
    chunk.extend_from_slice(passband);
    let rbb_full = downconvert_passband(&chunk, cfg);
    let rbb = &rbb_full[pre..];
    let coarse_sync_off = find_repeated_half_sync_offset(rbb, cfg);
    let (sync_off, time_scale) =
        forced_sync.unwrap_or_else(|| (refine_sync_offset(rbb, cfg, coarse_sync_off), 1.0));
    let rbb_sync = resample_from_offset_rate(rbb, sync_off, time_scale);
    let cfo_hz = estimate_coarse_cfo_hz(&rbb_sync, cfg);
    let rbb_cfo = coarse_cfo_correct(&rbb_sync, cfg);

    let (used_bins, pilot_bins, data_bins) = ofdm_bin_plan(cfg);
    let xsync_len = 2 * cfg.sync_half_len;
    let train_len = cfg.nfft + cfg.ncp;
    if rbb_cfo.len() < xsync_len + train_len || used_bins.is_empty() {
        return PassbandDiagnostics {
            enough_samples: false,
            sync_off: sync_off.round() as usize,
            cfo_hz,
            sync_rms: 0.0,
            sync_peak: 0.0,
            post_rms: 0.0,
            post_peak: 0.0,
            train_rms: 0.0,
            hest_mag_min: 0.0,
            hest_mag_mean: 0.0,
            hest_mag_max: 0.0,
            train_recon_evm: 0.0,
            pilot_residual_evm: 0.0,
            pilot_post_evm: 0.0,
            post_eq_evm: 0.0,
            decoded: false,
            decoded_payload_len: None,
        };
    }

    let train_start = xsync_len;
    let sync_slice = &rbb_cfo[..xsync_len.min(rbb_cfo.len())];
    let data_post_start = (xsync_len + train_len).min(rbb_cfo.len());
    let post_slice = &rbb_cfo[data_post_start..];
    let sync_rms = if sync_slice.is_empty() {
        0.0
    } else {
        (sync_slice.iter().map(|v| v.norm_sqr()).sum::<f32>() / sync_slice.len() as f32).sqrt()
    };
    let sync_peak = sync_slice.iter().map(|v| v.norm()).fold(0.0f32, f32::max);
    let post_rms = if post_slice.is_empty() {
        0.0
    } else {
        (post_slice.iter().map(|v| v.norm_sqr()).sum::<f32>() / post_slice.len() as f32).sqrt()
    };
    let post_peak = post_slice.iter().map(|v| v.norm()).fold(0.0f32, f32::max);
    let train_no_cp = &rbb_cfo[train_start + cfg.ncp..train_start + cfg.ncp + cfg.nfft];
    let train_rms =
        (train_no_cp.iter().map(|v| v.norm_sqr()).sum::<f32>() / (train_no_cp.len() as f32)).sqrt();
    let ytrain = fft(train_no_cp);
    let train_known = known_training_symbols(used_bins.len(), cfg.modulation);
    let mut hest = equalizer_initial_channel(cfg, &ytrain, &used_bins, &train_known);
    let mut eq_state = EqualizerTrackingState::default();
    equalizer_reset_tracking(&mut eq_state);
    let mut mags = hest.iter().map(|h| h.norm()).collect::<Vec<_>>();
    let mut ytrain_eq = Vec::with_capacity(used_bins.len());
    for (k, &bin) in used_bins.iter().enumerate() {
        ytrain_eq.push(regularized_equalize(ytrain[bin], hest[k]));
    }
    let train_recon_evm = rms_evm(&ytrain_eq, &train_known);
    let data_start = xsync_len + train_len;
    let sym_len = cfg.nfft + cfg.ncp;
    let max_payload_bytes = cfg.packet_payload_bytes + 16;
    let max_bits = fec_encoded_bits_len(max_payload_bytes * 8, cfg.fec_mode);
    let max_data_ofdm =
        max_bits.div_ceil(data_bins.len().max(1) * cfg.modulation.bits_per_symbol()) + 2;
    let symbol_plan = packet_symbol_plan(max_data_ofdm, cfg);
    let mut pilot_eq = Vec::new();
    let mut pilot_ref = Vec::new();
    let mut pilot_eq_post = Vec::new();
    let mut post_eq_data = Vec::new();
    let mut data_symbol_idx = 0usize;
    for (sym_idx, kind) in symbol_plan.into_iter().enumerate() {
        let s0 = data_start + sym_idx * sym_len;
        let s1 = s0 + sym_len;
        if s1 > rbb_cfo.len() || data_symbol_idx >= 3 {
            break;
        }
        let y = fft(&rbb_cfo[s0 + cfg.ncp..s0 + cfg.ncp + cfg.nfft]);
        if kind == PacketSymbolKind::Training {
            equalizer_refresh_channel(cfg, &mut hest, &y, &used_bins, &train_known);
            equalizer_reset_tracking(&mut eq_state);
            mags.extend(hest.iter().map(|h| h.norm()));
            continue;
        }
        let pref = known_pilot_symbols(pilot_bins.len(), data_symbol_idx + 1);
        let xeq_used = equalize_symbol_with_pilots(
            cfg,
            &mut eq_state,
            &y,
            &used_bins,
            &pilot_bins,
            &pref,
            &hest,
        );
        if !pilot_bins.is_empty() {
            for (k, pbin) in pilot_bins.iter().enumerate() {
                if let Some(pos) = used_bins.iter().position(|b| b == pbin) {
                    pilot_eq.push(xeq_used[pos]);
                    pilot_ref.push(pref[k]);
                    pilot_eq_post.push(xeq_used[pos]);
                }
            }
        }
        for dbin in &data_bins {
            if let Some(pos) = used_bins.iter().position(|b| b == dbin) {
                post_eq_data.push(xeq_used[pos]);
            }
        }
        data_symbol_idx += 1;
    }
    let hest_mag_min = mags.iter().copied().fold(f32::INFINITY, f32::min);
    let hest_mag_max = mags.iter().copied().fold(0.0f32, f32::max);
    let hest_mag_mean = mags.iter().sum::<f32>() / (mags.len() as f32);
    let pilot_residual_evm = rms_evm(&pilot_eq, &pilot_ref);
    let pilot_post_evm = rms_evm(&pilot_eq_post, &pilot_ref);
    let post_eq_evm = decision_directed_evm(&post_eq_data, cfg.modulation);
    let decoded = decode_packet_info_baseband(&rbb_cfo, cfg);

    PassbandDiagnostics {
        enough_samples: true,
        sync_off: sync_off.round() as usize,
        cfo_hz,
        sync_rms,
        sync_peak,
        post_rms,
        post_peak,
        train_rms,
        hest_mag_min,
        hest_mag_mean,
        hest_mag_max,
        train_recon_evm,
        pilot_residual_evm,
        pilot_post_evm,
        post_eq_evm,
        decoded: decoded.is_some(),
        decoded_payload_len: decoded.map(|p| p.payload.len()),
    }
}

/// Extracts equalizer input/output constellation samples for one passband window.
///
/// Parameters:
/// - `pkt_audio`: passband packet waveform window.
/// - `cfg`: modem configuration.
/// Returns:
/// - `Option<PassbandConstellationDump>`: constellation samples when extraction succeeds.
pub(crate) fn dump_passband_constellation_impl(
    pkt_audio: &[f32],
    cfg: &OfdmConfig,
    sync_off: f32,
) -> Option<PassbandConstellationDump> {
    let wake_len = (cfg.wake_ms * 1e-3 * cfg.fs) as usize;
    let guard_len = (cfg.wake_guard_ms * 1e-3 * cfg.fs) as usize;
    if pkt_audio.len() <= wake_len + guard_len + cfg.nfft + cfg.ncp {
        return None;
    }

    let passband = &pkt_audio[wake_len + guard_len..];
    let pre = 128usize.min(passband.len());
    let mut chunk = vec![0.0f32; pre];
    chunk.extend_from_slice(passband);
    let rbb_full = downconvert_passband(&chunk, cfg);
    let rbb = &rbb_full[pre..];
    let rbb_sync = resample_from_offset(rbb, sync_off);
    let rbb_cfo = coarse_cfo_correct(&rbb_sync, cfg);

    let (used_bins, pilot_bins, data_bins) = ofdm_bin_plan(cfg);
    if used_bins.is_empty() || data_bins.is_empty() {
        return None;
    }
    let xsync_len = 2 * cfg.sync_half_len;
    let train_len = cfg.nfft + cfg.ncp;
    if rbb_cfo.len() < xsync_len + train_len + cfg.nfft + cfg.ncp {
        return None;
    }

    let train_start = xsync_len;
    let train_no_cp = &rbb_cfo[train_start + cfg.ncp..train_start + cfg.ncp + cfg.nfft];
    let ytrain = fft(train_no_cp);
    let train_known = known_training_symbols(used_bins.len(), cfg.modulation);
    let mut hest = equalizer_initial_channel(cfg, &ytrain, &used_bins, &train_known);
    let mut eq_state = EqualizerTrackingState::default();
    equalizer_reset_tracking(&mut eq_state);

    let data_start = xsync_len + train_len;
    let sym_len = cfg.nfft + cfg.ncp;
    let mut pre_eq = Vec::new();
    let mut post_eq = Vec::new();
    let max_payload_bytes = cfg.packet_payload_bytes + 16;
    let max_bits = fec_encoded_bits_len(max_payload_bytes * 8, cfg.fec_mode);
    let max_data_ofdm =
        max_bits.div_ceil(data_bins.len().max(1) * cfg.modulation.bits_per_symbol()) + 2;
    let symbol_plan = packet_symbol_plan(max_data_ofdm, cfg);
    let mut data_symbol_idx = 0usize;
    for (sym_idx, kind) in symbol_plan.into_iter().enumerate() {
        let s0 = data_start + sym_idx * sym_len;
        let s1 = s0 + sym_len;
        if s1 > rbb_cfo.len() || data_symbol_idx >= 8 {
            break;
        }
        let y = fft(&rbb_cfo[s0 + cfg.ncp..s0 + cfg.ncp + cfg.nfft]);
        if kind == PacketSymbolKind::Training {
            equalizer_refresh_channel(cfg, &mut hest, &y, &used_bins, &train_known);
            equalizer_reset_tracking(&mut eq_state);
            continue;
        }
        let pref = known_pilot_symbols(pilot_bins.len(), data_symbol_idx + 1);
        let xeq_used = equalize_symbol_with_pilots(
            cfg,
            &mut eq_state,
            &y,
            &used_bins,
            &pilot_bins,
            &pref,
            &hest,
        );
        for dbin in &data_bins {
            if let Some(pos) = used_bins.iter().position(|b| b == dbin) {
                pre_eq.push(y[*dbin]);
                post_eq.push(xeq_used[pos]);
            }
        }
        data_symbol_idx += 1;
    }

    Some(PassbandConstellationDump { pre_eq, post_eq })
}

pub(crate) fn dump_passband_pilot_tracking_impl(
    pkt_audio: &[f32],
    cfg: &OfdmConfig,
) -> Option<PassbandPilotTrackDump> {
    let wake_len = (cfg.wake_ms * 1e-3 * cfg.fs) as usize;
    let guard_len = (cfg.wake_guard_ms * 1e-3 * cfg.fs) as usize;
    if pkt_audio.len() <= wake_len + guard_len + cfg.nfft + cfg.ncp {
        return None;
    }

    let passband = &pkt_audio[wake_len + guard_len..];
    let pre = 128usize.min(passband.len());
    let mut chunk = vec![0.0f32; pre];
    chunk.extend_from_slice(passband);
    let rbb_full = downconvert_passband(&chunk, cfg);
    let rbb = &rbb_full[pre..];
    let coarse_sync_off = find_repeated_half_sync_offset(rbb, cfg);
    let sync_off = refine_sync_offset(rbb, cfg, coarse_sync_off);
    let rbb_sync = resample_from_offset(rbb, sync_off);
    let rbb_cfo = coarse_cfo_correct(&rbb_sync, cfg);

    let (used_bins, pilot_bins, data_bins) = ofdm_bin_plan(cfg);
    if used_bins.is_empty() || pilot_bins.is_empty() || data_bins.is_empty() {
        return Some(PassbandPilotTrackDump {
            pilot_phase_rad: Vec::new(),
            pilot_evm_pre: Vec::new(),
            pilot_evm_post: Vec::new(),
            hest_mag_mean: Vec::new(),
            hest_mag_max: Vec::new(),
        });
    }
    let xsync_len = 2 * cfg.sync_half_len;
    let train_len = cfg.nfft + cfg.ncp;
    if rbb_cfo.len() < xsync_len + train_len + cfg.nfft + cfg.ncp {
        return None;
    }

    let train_start = xsync_len;
    let train_no_cp = &rbb_cfo[train_start + cfg.ncp..train_start + cfg.ncp + cfg.nfft];
    let ytrain = fft(train_no_cp);
    let train_known = known_training_symbols(used_bins.len(), cfg.modulation);
    let mut hest = equalizer_initial_channel(cfg, &ytrain, &used_bins, &train_known);
    let mut eq_state = EqualizerTrackingState::default();
    equalizer_reset_tracking(&mut eq_state);

    let data_start = xsync_len + train_len;
    let sym_len = cfg.nfft + cfg.ncp;
    let max_payload_bytes = cfg.packet_payload_bytes + 16;
    let max_bits = fec_encoded_bits_len(max_payload_bytes * 8, cfg.fec_mode);
    let max_data_ofdm =
        max_bits.div_ceil(data_bins.len().max(1) * cfg.modulation.bits_per_symbol()) + 2;
    let symbol_plan = packet_symbol_plan(max_data_ofdm, cfg);
    let mut data_symbol_idx = 0usize;
    let mut pilot_phase_rad = Vec::new();
    let mut pilot_evm_pre = Vec::new();
    let mut pilot_evm_post = Vec::new();
    let mut hest_mag_mean = Vec::new();
    let mut hest_mag_max = Vec::new();
    for (sym_idx, kind) in symbol_plan.into_iter().enumerate() {
        let s0 = data_start + sym_idx * sym_len;
        let s1 = s0 + sym_len;
        if s1 > rbb_cfo.len() {
            break;
        }
        let y = fft(&rbb_cfo[s0 + cfg.ncp..s0 + cfg.ncp + cfg.nfft]);
        if kind == PacketSymbolKind::Training {
            equalizer_refresh_channel(cfg, &mut hest, &y, &used_bins, &train_known);
            equalizer_reset_tracking(&mut eq_state);
            let mags = hest.iter().map(|h| h.norm()).collect::<Vec<_>>();
            hest_mag_mean.push(mags.iter().sum::<f32>() / (mags.len() as f32));
            hest_mag_max.push(mags.iter().copied().fold(0.0f32, f32::max));
            continue;
        }
        let pref = known_pilot_symbols(pilot_bins.len(), data_symbol_idx + 1);
        let xeq_used = equalize_symbol_with_pilots(
            cfg,
            &mut eq_state,
            &y,
            &used_bins,
            &pilot_bins,
            &pref,
            &hest,
        );
        let phase = pilot_phase_error(&xeq_used, &used_bins, &pilot_bins, &pref).unwrap_or(0.0);
        let mut pilot_eq_pre = Vec::new();
        let mut pilot_ref = Vec::new();
        for (k, pbin) in pilot_bins.iter().enumerate() {
            if let Some(pos) = used_bins.iter().position(|b| b == pbin) {
                pilot_eq_pre.push(xeq_used[pos]);
                pilot_ref.push(pref[k]);
            }
        }
        let mut pilot_eq_post = Vec::new();
        for (k, pbin) in pilot_bins.iter().enumerate() {
            if let Some(pos) = used_bins.iter().position(|b| b == pbin) {
                let _ = k;
                pilot_eq_post.push(xeq_used[pos]);
            }
        }
        pilot_phase_rad.push(phase);
        pilot_evm_pre.push(rms_evm(&pilot_eq_pre, &pilot_ref));
        pilot_evm_post.push(rms_evm(&pilot_eq_post, &pilot_ref));
        data_symbol_idx += 1;
    }

    Some(PassbandPilotTrackDump {
        pilot_phase_rad,
        pilot_evm_pre,
        pilot_evm_post,
        hest_mag_mean,
        hest_mag_max,
    })
}

pub(crate) fn dump_passband_bins_impl(
    pkt_audio: &[f32],
    cfg: &OfdmConfig,
) -> Option<PassbandBinDump> {
    dump_passband_bins_with_sync_opt(pkt_audio, cfg, None)
}

pub(crate) fn dump_passband_bins_with_sync_impl(
    pkt_audio: &[f32],
    cfg: &OfdmConfig,
    sync_off: f32,
) -> Option<PassbandBinDump> {
    dump_passband_bins_with_sync_opt(pkt_audio, cfg, Some(sync_off))
}

fn dump_passband_bins_with_sync_opt(
    pkt_audio: &[f32],
    cfg: &OfdmConfig,
    forced_sync_off: Option<f32>,
) -> Option<PassbandBinDump> {
    let wake_len = (cfg.wake_ms * 1e-3 * cfg.fs) as usize;
    let guard_len = (cfg.wake_guard_ms * 1e-3 * cfg.fs) as usize;
    if pkt_audio.len() <= wake_len + guard_len + cfg.nfft + cfg.ncp {
        return None;
    }

    let passband = &pkt_audio[wake_len + guard_len..];
    let pre = 128usize.min(passband.len());
    let mut chunk = vec![0.0f32; pre];
    chunk.extend_from_slice(passband);
    let rbb_full = downconvert_passband(&chunk, cfg);
    let rbb = &rbb_full[pre..];
    let coarse_sync_off = find_repeated_half_sync_offset(rbb, cfg);
    let sync_off = forced_sync_off.unwrap_or_else(|| refine_sync_offset(rbb, cfg, coarse_sync_off));
    let rbb_sync = resample_from_offset(rbb, sync_off);
    let rbb_cfo = coarse_cfo_correct(&rbb_sync, cfg);

    let (used_bins, pilot_bins, data_bins) = ofdm_bin_plan(cfg);
    if used_bins.is_empty() || data_bins.is_empty() {
        return None;
    }
    let xsync_len = 2 * cfg.sync_half_len;
    let train_len = cfg.nfft + cfg.ncp;
    if rbb_cfo.len() < xsync_len + train_len + cfg.nfft + cfg.ncp {
        return None;
    }

    let train_start = xsync_len;
    let train_no_cp = &rbb_cfo[train_start + cfg.ncp..train_start + cfg.ncp + cfg.nfft];
    let ytrain = fft(train_no_cp);
    let train_known = known_training_symbols(used_bins.len(), cfg.modulation);
    let mut hest = equalizer_initial_channel(cfg, &ytrain, &used_bins, &train_known);
    let mut eq_state = EqualizerTrackingState::default();
    equalizer_reset_tracking(&mut eq_state);

    let data_start = xsync_len + train_len;
    let sym_len = cfg.nfft + cfg.ncp;
    let max_payload_bytes = cfg.packet_payload_bytes + 16;
    let max_bits = fec_encoded_bits_len(max_payload_bytes * 8, cfg.fec_mode);
    let max_data_ofdm =
        max_bits.div_ceil(data_bins.len().max(1) * cfg.modulation.bits_per_symbol()) + 2;
    let symbol_plan = packet_symbol_plan(max_data_ofdm, cfg);
    let mut data_symbol_idx = 0usize;
    let mut rows = Vec::new();
    for (sym_idx, kind) in symbol_plan.into_iter().enumerate() {
        let s0 = data_start + sym_idx * sym_len;
        let s1 = s0 + sym_len;
        if s1 > rbb_cfo.len() || data_symbol_idx >= 8 {
            break;
        }
        let y = fft(&rbb_cfo[s0 + cfg.ncp..s0 + cfg.ncp + cfg.nfft]);
        if kind == PacketSymbolKind::Training {
            equalizer_refresh_channel(cfg, &mut hest, &y, &used_bins, &train_known);
            equalizer_reset_tracking(&mut eq_state);
            continue;
        }
        let pref = known_pilot_symbols(pilot_bins.len(), data_symbol_idx + 1);
        let xeq_used = equalize_symbol_with_pilots(
            cfg,
            &mut eq_state,
            &y,
            &used_bins,
            &pilot_bins,
            &pref,
            &hest,
        );
        let mut pre_eq_used = Vec::with_capacity(used_bins.len());
        let mut post_eq_used = Vec::with_capacity(used_bins.len());
        for (k, &bin) in used_bins.iter().enumerate() {
            pre_eq_used.push(y[bin]);
            post_eq_used.push(xeq_used[k]);
        }
        for (k, &ubin) in used_bins.iter().enumerate() {
            let role = if pilot_bins.contains(&ubin) {
                "pilot"
            } else {
                "data"
            };
            let reference = if role == "pilot" {
                pilot_bins
                    .iter()
                    .position(|b| *b == ubin)
                    .map(|pos| pref[pos])
            } else {
                None
            };
            rows.push(PassbandBinDumpRow {
                data_symbol_idx: data_symbol_idx + 1,
                used_bin: ubin,
                role,
                pre_eq: pre_eq_used[k],
                post_eq: post_eq_used[k],
                reference,
            });
        }
        data_symbol_idx += 1;
    }

    Some(PassbandBinDump { rows })
}

pub(crate) fn dump_passband_channel_compare_with_sync_impl(
    payload: &[u8],
    pkt_audio: &[f32],
    cfg: &OfdmConfig,
    sync_off: f32,
) -> Option<PassbandChannelCompareDump> {
    let wake_len = (cfg.wake_ms * 1e-3 * cfg.fs) as usize;
    let guard_len = (cfg.wake_guard_ms * 1e-3 * cfg.fs) as usize;
    if pkt_audio.len() <= wake_len + guard_len + cfg.nfft + cfg.ncp {
        return None;
    }

    let pkt_bytes = build_packet_bytes(payload, 0, 1, cfg);
    let xbb = tx_one_packet_baseband(&pkt_bytes, cfg);
    let passband = &pkt_audio[wake_len + guard_len..];
    let pre = 128usize.min(passband.len());
    let mut chunk = vec![0.0f32; pre];
    chunk.extend_from_slice(passband);
    let rbb_full = downconvert_passband(&chunk, cfg);
    let rbb = &rbb_full[pre..];
    let rbb_sync = resample_from_offset(rbb, sync_off);
    let rbb_cfo = coarse_cfo_correct(&rbb_sync, cfg);

    let (used_bins, pilot_bins, data_bins) = ofdm_bin_plan(cfg);
    if used_bins.is_empty() || data_bins.is_empty() {
        return None;
    }
    let xsync_len = 2 * cfg.sync_half_len;
    let train_len = cfg.nfft + cfg.ncp;
    if rbb_cfo.len() < xsync_len + train_len + cfg.nfft + cfg.ncp
        || xbb.len() < xsync_len + train_len
    {
        return None;
    }

    let train_start = xsync_len;
    let train_no_cp = &rbb_cfo[train_start + cfg.ncp..train_start + cfg.ncp + cfg.nfft];
    let ytrain = fft(train_no_cp);
    let train_known = known_training_symbols(used_bins.len(), cfg.modulation);
    let mut hest = equalizer_initial_channel(cfg, &ytrain, &used_bins, &train_known);
    let mut eq_state = EqualizerTrackingState::default();
    equalizer_reset_tracking(&mut eq_state);

    let data_start = xsync_len + train_len;
    let sym_len = cfg.nfft + cfg.ncp;
    let max_payload_bytes = cfg.packet_payload_bytes + 16;
    let max_bits = fec_encoded_bits_len(max_payload_bytes * 8, cfg.fec_mode);
    let max_data_ofdm =
        max_bits.div_ceil(data_bins.len().max(1) * cfg.modulation.bits_per_symbol()) + 2;
    let symbol_plan = packet_symbol_plan(max_data_ofdm, cfg);
    let mut data_symbol_idx = 0usize;
    let mut rows = Vec::new();
    for (sym_idx, kind) in symbol_plan.into_iter().enumerate() {
        let s0 = data_start + sym_idx * sym_len;
        let s1 = s0 + sym_len;
        if s1 > rbb_cfo.len() || s1 > xbb.len() || data_symbol_idx >= 8 {
            break;
        }
        let y = fft(&rbb_cfo[s0 + cfg.ncp..s0 + cfg.ncp + cfg.nfft]);
        let x = fft(&xbb[s0 + cfg.ncp..s0 + cfg.ncp + cfg.nfft]);
        if kind == PacketSymbolKind::Training {
            equalizer_refresh_channel(cfg, &mut hest, &y, &used_bins, &train_known);
            equalizer_reset_tracking(&mut eq_state);
            continue;
        }
        let pref = known_pilot_symbols(pilot_bins.len(), data_symbol_idx + 1);
        let xeq_used = equalize_symbol_with_pilots(
            cfg,
            &mut eq_state,
            &y,
            &used_bins,
            &pilot_bins,
            &pref,
            &hest,
        );
        let mut phase_by_bin = vec![0.0f32; used_bins.len()];
        if !pilot_bins.is_empty() && !pref.is_empty() {
            let mut pilot_phase_pts = Vec::<(f32, f32)>::new();
            for (k, pbin) in pilot_bins.iter().enumerate() {
                if let Some(pos) = used_bins.iter().position(|b| b == pbin) {
                    let z0 = regularized_equalize(y[*pbin], hest[pos]);
                    let ref_sym = pref[k];
                    if ref_sym.norm_sqr() > 1.0e-9 {
                        pilot_phase_pts.push((*pbin as f32, (z0 * ref_sym.conj()).arg()));
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
                        phase_by_bin[k] = intercept + slope * (bin as f32);
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
        let g = if den > 1.0e-9 {
            num / den
        } else {
            Complex32::new(1.0, 0.0)
        };
        for (k, &ubin) in used_bins.iter().enumerate() {
            let xref = x[ubin];
            if xref.norm() <= 1.0e-9 {
                continue;
            }
            let role = if pilot_bins.contains(&ubin) {
                "pilot"
            } else {
                "data"
            };
            rows.push(PassbandChannelCompareRow {
                data_symbol_idx: data_symbol_idx + 1,
                used_bin: ubin,
                role,
                actual_h: y[ubin] / xref,
                estimated_h_train: hest[k],
                estimated_h_pilot: hest[k] * Complex32::from_polar(1.0, phase_by_bin[k]) * g,
            });
        }
        data_symbol_idx += 1;
    }

    Some(PassbandChannelCompareDump { rows })
}

/// Extracts Schmidl-Cox timing metrics for one passband window.
///
/// Parameters:
/// - `pkt_audio`: passband packet waveform window.
/// - `cfg`: modem configuration.
/// Returns:
/// - `Option<PassbandSyncDump>`: sync metric samples and chosen offsets.
pub(crate) fn dump_passband_sync_metric_impl(
    pkt_audio: &[f32],
    cfg: &OfdmConfig,
) -> Option<PassbandSyncDump> {
    let wake_len = (cfg.wake_ms * 1e-3 * cfg.fs) as usize;
    let guard_len = (cfg.wake_guard_ms * 1e-3 * cfg.fs) as usize;
    if pkt_audio.len() <= wake_len + guard_len + cfg.nfft + cfg.ncp {
        return None;
    }

    let passband = &pkt_audio[wake_len + guard_len..];
    let pre = 128usize.min(passband.len());
    let mut chunk = vec![0.0f32; pre];
    chunk.extend_from_slice(passband);
    let rbb_full = downconvert_passband(&chunk, cfg);
    let rbb = &rbb_full[pre..];
    let metrics = repeated_half_sync_metrics(rbb, cfg);
    let coarse_sync_off = find_repeated_half_sync_offset(rbb, cfg);
    let refined_sync_off = refine_sync_offset(rbb, cfg, coarse_sync_off);
    Some(PassbandSyncDump {
        coarse_sync_off,
        refined_sync_off: refined_sync_off.round() as usize,
        metrics,
    })
}

pub(crate) fn dump_passband_iq_chain_impl(
    pkt_audio: &[f32],
    cfg: &OfdmConfig,
) -> Option<PassbandIqChainDump> {
    let wake_len = (cfg.wake_ms * 1e-3 * cfg.fs) as usize;
    let guard_len = (cfg.wake_guard_ms * 1e-3 * cfg.fs) as usize;
    if pkt_audio.len() <= wake_len + guard_len + cfg.nfft + cfg.ncp {
        return None;
    }

    let passband = &pkt_audio[wake_len + guard_len..];
    let pre = 128usize.min(passband.len());
    let mut chunk = vec![0.0f32; pre];
    chunk.extend_from_slice(passband);
    let down_audio = iq_downconvert(&chunk, cfg.fs, cfg.fc, passband_lpf_cutoff_hz(cfg));
    let baseband = match cfg.passband_mode {
        PassbandMode::Legacy => down_audio.clone(),
        PassbandMode::Iq => resample_complex_linear_rate(&down_audio, cfg.fs, cfg.fs_baseband),
    };
    let pre_baseband = match cfg.passband_mode {
        PassbandMode::Legacy => pre,
        PassbandMode::Iq => ((pre as f32) * (cfg.fs_baseband / cfg.fs)).round() as usize,
    };
    Some(PassbandIqChainDump {
        downconverted_audio_rate: down_audio[pre..].to_vec(),
        baseband_rate: baseband[pre_baseband.min(baseband.len())..].to_vec(),
        fs_audio: cfg.fs,
        fs_baseband: active_baseband_fs(cfg),
    })
}

/// Encodes one packet into passband samples (wake + guard + OFDM body).
///
/// Parameters:
/// - `pkt_bytes`: serialized packet bytes.
/// - `cfg`: modem configuration.
/// Returns:
/// - `Vec<f32>`: real passband waveform for one packet.
fn tx_one_packet(pkt_bytes: &[u8], cfg: &OfdmConfig) -> Vec<f32> {
    let passband = tx_one_packet_body(pkt_bytes, cfg);
    let mut wake = make_wake_tone(cfg);
    normalize_in_place(&mut wake, 0.25);
    let guard_len = (cfg.wake_guard_ms * 1e-3 * cfg.fs) as usize;

    let mut out = Vec::with_capacity(wake.len() + guard_len + passband.len());
    out.extend_from_slice(&wake);
    out.extend(std::iter::repeat(0.0).take(guard_len));
    out.extend(passband);
    out
}

fn tx_one_packet_body(pkt_bytes: &[u8], cfg: &OfdmConfig) -> Vec<f32> {
    let xbb = tx_one_packet_baseband(pkt_bytes, cfg);
    let mut passband = upconvert_passband(&xbb, cfg);
    if (cfg.payload_gain - 1.0).abs() > f32::EPSILON {
        for s in &mut passband {
            *s *= cfg.payload_gain;
        }
    }
    let sync_len = (2 * cfg.sync_half_len).min(passband.len());
    if sync_len > 0 {
        normalize_in_place(&mut passband[..sync_len], 0.35);
    }
    if sync_len < passband.len() {
        normalize_in_place(&mut passband[sync_len..], 0.96);
    }
    passband
}

fn upconvert_passband(xbb: &[Complex32], cfg: &OfdmConfig) -> Vec<f32> {
    match cfg.passband_mode {
        PassbandMode::Legacy => iq_upconvert(xbb, cfg.fs, cfg.fc),
        PassbandMode::Iq => {
            let xbb_audio = resample_complex_linear_rate(xbb, cfg.fs_baseband, cfg.fs);
            iq_upconvert(&xbb_audio, cfg.fs, cfg.fc)
        }
    }
}

fn passband_lpf_cutoff_hz(cfg: &OfdmConfig) -> f32 {
    let (used_bins, _, _) = ofdm_bin_plan(cfg);
    let max_bin = used_bins.iter().copied().max().unwrap_or(1) as f32;
    let bw = ((max_bin + 2.0) * active_baseband_fs(cfg) / cfg.nfft as f32).max(500.0);
    bw.min(0.45 * cfg.fs).min(0.45 * active_baseband_fs(cfg))
}

fn downconvert_passband(pkt_audio: &[f32], cfg: &OfdmConfig) -> Vec<Complex32> {
    let mixed = iq_downconvert(pkt_audio, cfg.fs, cfg.fc, passband_lpf_cutoff_hz(cfg));
    match cfg.passband_mode {
        PassbandMode::Legacy => mixed,
        PassbandMode::Iq => resample_complex_linear_rate(&mixed, cfg.fs, cfg.fs_baseband),
    }
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
    let rbb_full = downconvert_passband(&chunk, cfg);
    let rbb = &rbb_full[pre..];
    let coarse_sync_off = find_repeated_half_sync_offset(rbb, cfg);
    let sync_off = refine_sync_offset(rbb, cfg, coarse_sync_off);
    for off in [sync_off, coarse_sync_off as f32] {
        let rbb_sync = resample_from_offset(rbb, off);
        if let Some(pkt) = decode_packet_from_synced_baseband(&rbb_sync, cfg) {
            return Some(pkt);
        }
    }
    None
}

fn decode_packet_from_passband_with_sync_rate(
    pkt_audio: &[f32],
    cfg: &OfdmConfig,
    sync_off: f32,
    time_scale: f32,
) -> Option<PacketInfo> {
    let wake_len = (cfg.wake_ms * 1e-3 * cfg.fs) as usize;
    let guard_len = (cfg.wake_guard_ms * 1e-3 * cfg.fs) as usize;
    if pkt_audio.len() <= wake_len + guard_len + cfg.nfft + cfg.ncp {
        return None;
    }

    let passband = &pkt_audio[wake_len + guard_len..];
    let pre = 128usize.min(passband.len());
    let mut chunk = vec![0.0f32; pre];
    chunk.extend_from_slice(passband);
    let rbb_full = downconvert_passband(&chunk, cfg);
    let rbb = &rbb_full[pre..];
    let rbb_sync = resample_from_offset_rate(rbb, sync_off, time_scale);
    decode_packet_from_synced_baseband(&rbb_sync, cfg)
}

fn decode_packet_from_synced_baseband(
    rbb_sync: &[Complex32],
    cfg: &OfdmConfig,
) -> Option<PacketInfo> {
    let coarse_cfo_hz = estimate_coarse_cfo_hz(rbb_sync, cfg);
    let mut tried = Vec::new();
    tried.push(coarse_cfo_hz);
    for hz in [-12.0f32, -8.0, -4.0, 4.0, 8.0, 12.0] {
        tried.push(coarse_cfo_hz + hz);
    }
    for cfo_hz in tried {
        let rbb_cfo = apply_cfo_hz(rbb_sync, active_baseband_fs(cfg), cfo_hz);
        if let Some(pkt) = decode_packet_info_baseband(&rbb_cfo, cfg) {
            return Some(pkt);
        }
    }
    None
}

/// Finds sync start using Schmidl-Cox metric on repeated-half preamble.
///
/// Parameters:
/// - `rbb`: baseband complex samples near packet start.
/// - `cfg`: modem configuration.
/// Returns:
/// - `usize`: estimated sync start offset in samples (relative to `rbb`).
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
/// - `cutoff_hz`: low-pass cutoff (Hz).
/// Returns:
/// - `Vec<Complex32>`: filtered complex baseband samples.
fn iq_downconvert(y: &[f32], fs: f32, fc: f32, cutoff_hz: f32) -> Vec<Complex32> {
    let mixed: Vec<Complex32> = y
        .iter()
        .enumerate()
        .map(|(n, &s)| {
            let ph = -2.0 * std::f32::consts::PI * fc * (n as f32) / fs;
            Complex32::new(s, 0.0) * Complex32::from_polar(1.0, ph)
        })
        .collect();
    let mut out = lowpass_fir(&mixed, fs, cutoff_hz.max(100.0), 97);
    for z in &mut out {
        *z *= 2.0;
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
    use crate::config::Modulation;
    use crate::equalizer::hard_slice_symbols;

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
    /// Verifies passband packet encode/decode path.
    fn decode_single_packet_passband_test() {
        let cfg = OfdmConfig::default();
        let payload: Vec<u8> = (40..64).collect();
        let y = encode_single_packet_passband(&payload, &cfg);
        let out =
            decode_single_packet_passband_with_sync(&y, &cfg, 0.0).expect("passband decode failed");
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
}

// vim: set ts=4 sw=4 et:
