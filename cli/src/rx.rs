// Copyright (c) 2026 Elias S. G. Carotti

use std::error::Error;
use std::io::Write;
use std::path::Path;
use std::sync::{Arc, Mutex};
use std::time::{Duration, Instant};

use acoustic_ofdm::{
    decode_single_packet_passband, diagnose_passband_window, dump_passband_bins,
    dump_passband_constellation, dump_passband_pilot_tracking, dump_passband_sync_metric,
    save_constellation_comparison_png, save_spectrogram_png, save_spectrogram_png_with_options,
    save_wav_mono_i16, OfdmConfig,
};
use cpal::traits::{DeviceTrait, HostTrait, StreamTrait};
use ringbuf::{traits::*, HeapRb};

use crate::audio::build_input_stream;
use crate::cli_args::{AudioOpts, ORACLE_PAYLOAD};
use crate::rx_support::{
    burst_active_regions, candidate_region_index, clipping_diag, diagnostic_candidate_score,
    estimated_packet_len_samples, filter_candidates_by_regions, filter_for_sync_detection,
    make_wake_ref, plausible_candidate, print_passband_diagnostics, quick_realtime_decode,
    ranked_offset_hypotheses, refine_wake_candidates_fractional, signal_diag, wake_candidates,
};
use crate::{debug_line, info_line, warn_line};

pub(crate) fn cmd_rx(
    cfg: &OfdmConfig,
    opts: &AudioOpts,
    stdout_raw: bool,
) -> Result<(), Box<dyn Error>> {
    let host = cpal::default_host();
    let in_dev = host.default_input_device().ok_or("no input device")?;
    let in_cfg = in_dev.default_input_config()?;
    let mut cfg_rt = cfg.clone();
    cfg_rt.fs = in_cfg.sample_rate().0 as f32;
    cfg_rt.sync_half_len = ((0.25 * cfg_rt.fs * 0.5).round() as usize).max(64);
    cfg_rt.use_pilots = Some(true);

    let in_cap = (opts.duration_sec * cfg_rt.fs).ceil() as usize + 4096;
    let in_rb = HeapRb::<f32>::new(in_cap.max(4096));
    let (in_prod, mut in_cons) = in_rb.split();
    let in_prod = Arc::new(Mutex::new(in_prod));
    let in_stream = build_input_stream(
        &in_dev,
        &in_cfg.clone().into(),
        in_cfg.sample_format(),
        opts.mic_gain,
        in_prod,
    )?;

    info_line!("Input device : {}", in_dev.name()?);
    info_line!(
        "Stream config: {} Hz, in {:?}",
        in_cfg.sample_rate().0,
        in_cfg.sample_format()
    );
    info_line!("RX detector: wake-correlation-v2");
    info_line!("Wake preamble: {}", cfg_rt.wake_preamble.as_str());
    if opts.oracle {
        info_line!(
            "Oracle mode: enabled (expect {} bytes)",
            ORACLE_PAYLOAD.len()
        );
    }
    info_line!(
        "RX sync filter: {} (hp={:.1}Hz, lp={:.1}Hz)",
        if opts.input_filter { "on" } else { "off" },
        opts.input_hp_hz,
        opts.input_lp_hz
    );
    info_line!("Listening for {:.2}s ...", opts.duration_sec);
    in_stream.play()?;
    let mut rx = Vec::<f32>::new();
    let t0 = Instant::now();
    let deadline = t0 + Duration::from_secs_f32(opts.duration_sec);
    let mut last_rt_try = Instant::now();
    while Instant::now() < deadline {
        while let Some(s) = in_cons.try_pop() {
            rx.push(s);
        }
        if Instant::now() >= deadline {
            break;
        }
        if last_rt_try.elapsed() >= Duration::from_millis(300) {
            let tail_span = ((cfg_rt.fs * 3.0).round() as usize).max(1);
            let st = rx.len().saturating_sub(tail_span);
            if let Some(bytes) = quick_realtime_decode(&rx[st..], &rx[st..], &cfg_rt) {
                let rt = t0.elapsed().as_secs_f32();
                info_line!("Realtime decode: OK at t={rt:.3}s ({} bytes)", bytes.len());
                if opts.oracle {
                    info_line!(
                        "Oracle verdict: {}",
                        if bytes.as_slice() == ORACLE_PAYLOAD {
                            "match"
                        } else {
                            "mismatch"
                        }
                    );
                }
                if stdout_raw {
                    let mut out = std::io::stdout().lock();
                    out.write_all(&bytes)?;
                    out.flush()?;
                } else {
                    let hex = bytes
                        .iter()
                        .map(|b| format!("{b:02X}"))
                        .collect::<Vec<_>>()
                        .join(" ");
                    info_line!("HEX: {hex}");
                    info_line!("UTF8(lossy): {}", String::from_utf8_lossy(&bytes));
                }
                drop(in_stream);
                return Ok(());
            }
            if opts.verbose {
                debug_line!(
                    "Realtime: checked t={:.3}s captured={}",
                    t0.elapsed().as_secs_f32().min(opts.duration_sec),
                    rx.len()
                );
            }
            last_rt_try = Instant::now();
        }
        if let Some(rem) = deadline.checked_duration_since(Instant::now()) {
            std::thread::sleep(rem.min(Duration::from_millis(10)));
        } else {
            break;
        }
    }
    drop(in_stream);

    let target_samples = (opts.duration_sec * cfg_rt.fs).round().max(0.0) as usize;
    while let Some(s) = in_cons.try_pop() {
        if rx.len() >= target_samples {
            break;
        }
        rx.push(s);
    }
    if rx.len() > target_samples {
        rx.truncate(target_samples);
    }
    info_line!("Captured samples: {}", rx.len());
    if let Some(path) = &opts.dump_wav {
        save_wav_mono_i16(Path::new(path), &rx, cfg_rt.fs.round() as u32)?;
        info_line!("Saved RX capture: {path}");
    }
    if opts.spectrogram {
        let spec_path = Path::new(&opts.spectrogram_path);
        let default_spec = acoustic_ofdm::SpectrogramOptions::default();
        if opts.spectrogram_opts.nfft == default_spec.nfft
            && opts.spectrogram_opts.hop == default_spec.hop
            && opts.spectrogram_opts.window == default_spec.window
        {
            save_spectrogram_png(spec_path, &rx, cfg_rt.fs)?;
        } else {
            save_spectrogram_png_with_options(spec_path, &rx, cfg_rt.fs, opts.spectrogram_opts)?;
        }
        info_line!("Saved spectrogram PNG: {}", spec_path.display());
    }
    let rx_sync = filter_for_sync_detection(
        &rx,
        cfg_rt.fs,
        opts.input_hp_hz,
        opts.input_lp_hz,
        opts.input_filter,
    );
    let (rms, peak, first_loud) = signal_diag(&rx);
    let (clipped, clipped_frac) = clipping_diag(&rx);
    if opts.verbose {
        debug_line!(
            "RX diagnostics: rms={rms:.5} peak={peak:.5} first_loud_sample={first_loud} ({:.3}s)",
            (first_loud as f32) / cfg_rt.fs
        );
        debug_line!(
            "RX clipping: {clipped} samples ({:.2}%) at |x| >= 0.995",
            100.0 * clipped_frac
        );
        if peak < 0.01 {
            warn_line!("RX warning: very low capture level; increase speaker volume or mic gain.");
        }
    }

    let packet_bytes = if opts.oracle {
        Some(11 + ORACLE_PAYLOAD.len())
    } else {
        None
    };
    let est_pkt = estimated_packet_len_samples(&cfg_rt, packet_bytes);
    let pad = ((0.050 * cfg_rt.fs).round() as usize).max(1);
    if opts.verbose {
        debug_line!(
            "Decode window: est_packet={est_pkt} samples ({:.3}s), pad={pad} samples",
            (est_pkt as f32) / cfg_rt.fs
        );
    }
    let first_end = (est_pkt + pad).min(rx.len());
    let mut dec = decode_single_packet_passband(&rx[..first_end], &cfg_rt);
    let mut attempts = 1usize;
    if dec.is_none() {
        let wake = make_wake_ref(&cfg_rt);
        let coarse_step = ((cfg_rt.fs * 0.0005).round() as usize).max(1);
        let gate = filter_for_sync_detection(&rx, cfg_rt.fs, 12_000.0, 19_000.0, true);
        let active_regions = burst_active_regions(&gate, cfg_rt.fs);
        let cands0 = filter_candidates_by_regions(
            &wake_candidates(&rx_sync, &wake, coarse_step, 12, est_pkt),
            &active_regions,
        );
        let cands = refine_wake_candidates_fractional(&rx_sync, &wake, &cands0);
        let keep_regions = active_regions.len().min(3);
        let mut ranked_cands: Vec<(usize, f32, f32, bool, usize)> = Vec::new();
        for (idx, wake_score) in cands.iter().copied().take(8) {
            let Some(region_idx) = candidate_region_index(idx, &active_regions, keep_regions)
            else {
                continue;
            };
            let end = idx.saturating_add(est_pkt + pad).min(rx.len());
            if end <= idx + cfg_rt.nfft + cfg_rt.ncp {
                continue;
            }
            let diag = diagnose_passband_window(&rx[idx..end], &cfg_rt);
            let diag_score = diagnostic_candidate_score(&diag);
            let plausible = plausible_candidate(&diag);
            let cand = (idx, wake_score, diag_score, plausible, region_idx);
            match ranked_cands
                .iter_mut()
                .find(|(_, _, _, _, ridx)| *ridx == region_idx)
            {
                Some(best) => {
                    if cand
                        .3
                        .cmp(&best.3)
                        .then_with(|| cand.2.total_cmp(&best.2))
                        .then_with(|| cand.1.total_cmp(&best.1))
                        .then_with(|| best.0.cmp(&cand.0))
                        .is_gt()
                    {
                        *best = cand;
                    }
                }
                None => ranked_cands.push(cand),
            }
        }
        let ranked_cands_simple: Vec<(usize, f32, f32, bool)> = ranked_cands
            .iter()
            .map(|(idx, wake_score, diag_score, plausible, _)| {
                (*idx, *wake_score, *diag_score, *plausible)
            })
            .collect();
        let local_ranked_simple: Vec<(usize, f32, f32, bool)> = cands
            .iter()
            .take(2)
            .filter_map(|(idx, wake_score)| {
                let end = idx.saturating_add(est_pkt + pad).min(rx.len());
                if end <= *idx + cfg_rt.nfft + cfg_rt.ncp {
                    return None;
                }
                let diag = diagnose_passband_window(&rx[*idx..end], &cfg_rt);
                Some((
                    *idx,
                    *wake_score,
                    diagnostic_candidate_score(&diag),
                    plausible_candidate(&diag),
                ))
            })
            .collect();
        let local_refined =
            ranked_offset_hypotheses(&rx, &local_ranked_simple, est_pkt, pad, &cfg_rt);
        let refined = if local_refined.is_empty() {
            ranked_offset_hypotheses(&rx, &ranked_cands_simple, est_pkt, pad, &cfg_rt)
        } else {
            local_refined
        };
        let primary_off = refined
            .first()
            .map(|(off, _)| *off)
            .or_else(|| ranked_cands.first().map(|(idx, _, _, _, _)| *idx));

        if opts.verbose {
            debug_line!(
                "Wake search: {} candidates (step={} samples, active_regions={})",
                cands.len(),
                coarse_step,
                active_regions.len()
            );
            for (i, r) in active_regions.iter().take(5).enumerate() {
                debug_line!(
                    "  region {:2}: [{:.3}s, {:.3}s] mean_rms={:.4} peak_rms={:.4}",
                    i + 1,
                    (r.start as f32) / cfg_rt.fs,
                    (r.end as f32) / cfg_rt.fs,
                    r.mean_rms,
                    r.peak_rms
                );
            }
            for (i, (idx, sc)) in cands.iter().enumerate() {
                debug_line!(
                    "  cand {:2}: idx={} t={:.3}s score={:.4}",
                    i + 1,
                    idx,
                    (*idx as f32) / cfg_rt.fs,
                    sc
                );
            }
            for (i, (idx, wake_score, diag_score, plausible, region_idx)) in
                ranked_cands.iter().take(3).enumerate()
            {
                debug_line!(
                    "  rank {:2}: idx={} t={:.3}s region={} wake_score={:.4} diag_score={:.4} plausible={}",
                    i + 1,
                    idx,
                    (*idx as f32) / cfg_rt.fs,
                    region_idx + 1,
                    wake_score,
                    diag_score,
                    plausible
                );
            }
            for (i, (idx, _, _, _, _)) in ranked_cands.iter().take(3).enumerate() {
                let off = *idx;
                let end = off.saturating_add(est_pkt + pad).min(rx.len());
                if end > off + cfg_rt.nfft + cfg_rt.ncp {
                    print_passband_diagnostics(
                        &format!("  cand {:2} diag", i + 1),
                        &rx[off..end],
                        &cfg_rt,
                    );
                }
            }
            if let Some(off) = primary_off {
                debug_line!(
                    "Primary hypothesis: off={} t={:.3}s",
                    off,
                    (off as f32) / cfg_rt.fs
                );
                let end = off.saturating_add(est_pkt + pad).min(rx.len());
                if end > off + cfg_rt.nfft + cfg_rt.ncp {
                    let (local_clipped, local_clipped_frac) = clipping_diag(&rx[off..end]);
                    debug_line!(
                        "Top-window clipping: {local_clipped} samples ({:.2}%) at |x| >= 0.995",
                        100.0 * local_clipped_frac
                    );
                    if local_clipped_frac >= 0.001 {
                        warn_line!("RX warning: selected packet window is clipping; reduce speaker volume, mic gain, or disable AGC.");
                    } else if local_clipped_frac >= 0.0001 {
                        warn_line!("RX warning: selected packet window is close to clipping.");
                    } else if clipped_frac >= 0.001 {
                        debug_line!("Global clipping exists outside the selected packet window; gain reduction may not be necessary for this decode.");
                    }
                    if let Some(dump) = dump_passband_constellation(&rx[off..end], &cfg_rt) {
                        let pre_path = Path::new("/tmp/ofdm_constellation_pre_eq.csv");
                        let post_path = Path::new("/tmp/ofdm_constellation_post_eq.csv");
                        let plot_png_path = Path::new("/tmp/ofdm_constellation.png");
                        {
                            let mut out = std::io::BufWriter::new(std::fs::File::create(pre_path)?);
                            writeln!(out, "re,im")?;
                            for z in &dump.pre_eq {
                                writeln!(out, "{},{}", z.re, z.im)?;
                            }
                            out.flush()?;
                        }
                        {
                            let mut out =
                                std::io::BufWriter::new(std::fs::File::create(post_path)?);
                            writeln!(out, "re,im")?;
                            for z in &dump.post_eq {
                                writeln!(out, "{},{}", z.re, z.im)?;
                            }
                            out.flush()?;
                        }
                        save_constellation_comparison_png(
                            plot_png_path,
                            &dump.pre_eq,
                            &dump.post_eq,
                        )?;
                        debug_line!("Saved constellation CSV: {}", pre_path.display());
                        debug_line!("Saved constellation CSV: {}", post_path.display());
                        debug_line!("Saved constellation PNG: {}", plot_png_path.display());
                    }
                    if let Some(bin_dump) = dump_passband_bins(&rx[off..end], &cfg_rt) {
                        let bins_path = Path::new("/tmp/ofdm_bins.csv");
                        let mut out = std::io::BufWriter::new(std::fs::File::create(bins_path)?);
                        writeln!(out, "data_symbol_idx,used_bin,role,pre_re,pre_im,post_re,post_im,ref_re,ref_im")?;
                        for row in &bin_dump.rows {
                            let (ref_re, ref_im) = row
                                .reference
                                .map(|z| (z.re, z.im))
                                .unwrap_or((f32::NAN, f32::NAN));
                            writeln!(
                                out,
                                "{},{},{},{},{},{},{},{},{}",
                                row.data_symbol_idx,
                                row.used_bin,
                                row.role,
                                row.pre_eq.re,
                                row.pre_eq.im,
                                row.post_eq.re,
                                row.post_eq.im,
                                ref_re,
                                ref_im
                            )?;
                        }
                        out.flush()?;
                        debug_line!("Saved bin dump CSV: {}", bins_path.display());
                    }
                    if let Some(track) = dump_passband_pilot_tracking(&rx[off..end], &cfg_rt) {
                        for (i, (hmean, hmax)) in track
                            .hest_mag_mean
                            .iter()
                            .zip(track.hest_mag_max.iter())
                            .take(8)
                            .enumerate()
                        {
                            debug_line!(
                                "Channel refresh {:2}: hest_mean={hmean:.3} hest_max={hmax:.3}",
                                i + 1
                            );
                        }
                        for (i, ((phase_rad, pilot_evm_pre), pilot_evm_post)) in track
                            .pilot_phase_rad
                            .iter()
                            .zip(track.pilot_evm_pre.iter())
                            .zip(track.pilot_evm_post.iter())
                            .take(8)
                            .enumerate()
                        {
                            debug_line!(
                                "Pilot track sym={:2}: phase={:.3}rad ({:.1}deg) pilot_evm_pre={:.3} pilot_evm_post={:.3}",
                                i + 1,
                                phase_rad,
                                phase_rad.to_degrees(),
                                pilot_evm_pre,
                                pilot_evm_post
                            );
                        }
                    }
                    if let Some(sync_dump) = dump_passband_sync_metric(&rx[off..end], &cfg_rt) {
                        let sync_path = Path::new("/tmp/ofdm_sync_metric.csv");
                        let mut out = std::io::BufWriter::new(std::fs::File::create(sync_path)?);
                        writeln!(out, "offset,metric")?;
                        for (i, m) in sync_dump.metrics.iter().enumerate() {
                            writeln!(out, "{i},{m}")?;
                        }
                        out.flush()?;
                        debug_line!(
                            "Saved sync metric CSV: {} (coarse_sync_off={}, refined_sync_off={})",
                            sync_path.display(),
                            sync_dump.coarse_sync_off,
                            sync_dump.refined_sync_off
                        );
                    }
                }
            }
            debug_line!(
                "Metric refinement ({}): {} shortlisted offsets from top {} seeds",
                if !local_ranked_simple.is_empty() {
                    "wake-anchored"
                } else {
                    "fallback"
                },
                refined.len(),
                if !local_ranked_simple.is_empty() {
                    local_ranked_simple.len()
                } else {
                    ranked_cands_simple.len().min(3)
                }
            );
            for (i, (off, diag)) in refined.iter().take(6).enumerate() {
                debug_line!(
                    "  refine {:2}: off={} t={:.3}s sync_off={} train_rms={:.4} evm[train/pilot_pre/pilot_post/data]=[{:.3}/{:.3}/{:.3}/{:.3}]",
                    i + 1,
                    off,
                    (*off as f32) / cfg_rt.fs,
                    diag.sync_off,
                    diag.train_rms,
                    diag.train_recon_evm,
                    diag.pilot_residual_evm,
                    diag.pilot_post_evm,
                    diag.post_eq_evm,
                );
            }
        }
        for (off, _diag) in refined.iter().take(12) {
            let end = off.saturating_add(est_pkt + pad).min(rx.len());
            if end <= *off + cfg_rt.nfft + cfg_rt.ncp {
                continue;
            }
            attempts += 1;
            if let Some(bytes) = decode_single_packet_passband(&rx[*off..end], &cfg_rt) {
                if opts.verbose {
                    debug_line!(
                        "Decode recovered at refined offset {} ({:.3}s)",
                        off,
                        (*off as f32) / cfg_rt.fs
                    );
                }
                dec = Some(bytes);
                break;
            }
        }
    }
    if opts.verbose {
        info_line!("Decode attempts: {attempts}");
    }

    match dec {
        Some(bytes) => {
            info_line!("Decode: OK ({} bytes)", bytes.len());
            if opts.oracle {
                info_line!(
                    "Oracle verdict: {}",
                    if bytes.as_slice() == ORACLE_PAYLOAD {
                        "match"
                    } else {
                        "mismatch"
                    }
                );
            }
            if stdout_raw {
                let mut out = std::io::stdout().lock();
                out.write_all(&bytes)?;
                out.flush()?;
            } else {
                let hex = bytes
                    .iter()
                    .map(|b| format!("{b:02X}"))
                    .collect::<Vec<_>>()
                    .join(" ");
                info_line!("HEX: {hex}");
                info_line!("UTF8(lossy): {}", String::from_utf8_lossy(&bytes));
            }
        }
        None => warn_line!("Decode: FAIL (sync/CFO/equalization/CRC path)"),
    }
    Ok(())
}
