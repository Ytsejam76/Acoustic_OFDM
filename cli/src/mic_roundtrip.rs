// Copyright (c) 2026 Elias S. G. Carotti

use std::error::Error;
use std::io::Write;
use std::path::Path;
use std::sync::{Arc, Mutex};
use std::time::{Duration, Instant};

use acoustic_ofdm::{
    decode_single_packet_passband_with_sync, diagnose_passband_window_with_sync,
    dump_passband_bins_with_sync, dump_passband_channel_compare_with_sync,
    dump_passband_constellation, dump_passband_iq_chain, inspect_packet_bytes,
    recover_decided_packet_bytes_passband_with_sync, save_channel_compare_png,
    save_constellation_comparison_png, save_spectrogram_png, save_spectrogram_png_with_options,
    save_wav_mono_i16, Complex32, OfdmConfig, PassbandBinDump, PassbandChannelCompareDump,
};
use cpal::traits::{DeviceTrait, HostTrait, StreamTrait};
use ringbuf::{traits::*, HeapRb};

use crate::audio::{build_input_stream, build_output_stream, clipping_diag, signal_diag};
use crate::cli_args::{AudioOpts, ORACLE_PAYLOAD};
use crate::tx::build_tx_plan;
use crate::{debug_line, info_line, warn_line};

fn save_bin_dump_csv(path: &Path, dump: &PassbandBinDump) -> Result<(), Box<dyn Error>> {
    let mut file = std::fs::File::create(path)?;
    writeln!(
        file,
        "data_symbol_idx,used_bin,role,pre_re,pre_im,post_re,post_im,ref_re,ref_im"
    )?;
    for row in &dump.rows {
        let (rr, ri) = row
            .reference
            .map(|z| (z.re, z.im))
            .unwrap_or((f32::NAN, f32::NAN));
        writeln!(
            file,
            "{},{},{},{},{},{},{},{},{}",
            row.data_symbol_idx,
            row.used_bin,
            row.role,
            row.pre_eq.re,
            row.pre_eq.im,
            row.post_eq.re,
            row.post_eq.im,
            rr,
            ri
        )?;
    }
    Ok(())
}

fn save_channel_compare_csv(
    path: &Path,
    dump: &PassbandChannelCompareDump,
) -> Result<(), Box<dyn Error>> {
    let mut file = std::fs::File::create(path)?;
    writeln!(
        file,
        "data_symbol_idx,used_bin,role,actual_re,actual_im,est_train_re,est_train_im,est_pilot_re,est_pilot_im"
    )?;
    for row in &dump.rows {
        writeln!(
            file,
            "{},{},{},{},{},{},{},{},{}",
            row.data_symbol_idx,
            row.used_bin,
            row.role,
            row.actual_h.re,
            row.actual_h.im,
            row.estimated_h_train.re,
            row.estimated_h_train.im,
            row.estimated_h_pilot.re,
            row.estimated_h_pilot.im
        )?;
    }
    Ok(())
}

fn save_complex_parts_wav(
    prefix: &Path,
    samples: &[Complex32],
    sample_rate: u32,
) -> Result<(), Box<dyn Error>> {
    let re = samples.iter().map(|z| z.re).collect::<Vec<_>>();
    let im = samples.iter().map(|z| z.im).collect::<Vec<_>>();
    let re_path = prefix.with_extension("re.wav");
    let im_path = prefix.with_extension("im.wav");
    save_wav_mono_i16(&re_path, &re, sample_rate)?;
    save_wav_mono_i16(&im_path, &im, sample_rate)?;
    Ok(())
}

pub(crate) fn cmd_mic_roundtrip(
    payload: &[u8],
    cfg: &OfdmConfig,
    opts: &AudioOpts,
    stdout_raw: bool,
) -> Result<(), Box<dyn Error>> {
    let host = cpal::default_host();
    let in_dev = host.default_input_device().ok_or("no input device")?;
    let in_cfg = in_dev.default_input_config()?;
    let out_dev = host.default_output_device().ok_or("no output device")?;
    let out_cfg = out_dev.default_output_config()?;

    let tx_plan = build_tx_plan(payload, cfg, opts, out_cfg.sample_rate().0 as f32)?;
    let mut cfg_rx = cfg.clone();
    cfg_rx.fs = in_cfg.sample_rate().0 as f32;

    let tx_total_sec = tx_plan.scheduled.len() as f32 / tx_plan.fs;
    let capture_sec = opts.duration_sec.max(tx_total_sec + 0.5);
    let in_cap = (capture_sec * cfg_rx.fs).ceil() as usize + 4096;
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
    let out_stream = build_output_stream(
        &out_dev,
        &out_cfg.clone().into(),
        out_cfg.sample_format(),
        tx_plan.scheduled.clone(),
    )?;

    info_line!("Input device : {}", in_dev.name()?);
    info_line!(
        "Input config : {} Hz, {} ch, in {:?}",
        in_cfg.sample_rate().0,
        in_cfg.channels(),
        in_cfg.sample_format()
    );
    info_line!("Output device: {}", out_dev.name()?);
    info_line!(
        "Output config: {} Hz, {} ch, out {:?}",
        out_cfg.sample_rate().0,
        out_cfg.channels(),
        out_cfg.sample_format()
    );
    info_line!("Wake preamble: {}", cfg_rx.wake_preamble.as_str());
    if opts.oracle {
        info_line!(
            "Oracle mode: enabled (expect {} bytes)",
            ORACLE_PAYLOAD.len()
        );
    }
    info_line!(
        "Mic roundtrip: capture={capture_sec:.2}s packet={:.3}s repeats={} pre_delay={:.2}s gap={:.2}s",
        tx_plan.packet.len() as f32 / tx_plan.fs,
        opts.repeats,
        opts.pre_delay_sec,
        opts.gap_sec
    );

    in_stream.play()?;
    std::thread::sleep(Duration::from_millis(50));
    let tx_start_in_capture = 0.05f32;
    out_stream.play()?;

    let deadline = Instant::now() + Duration::from_secs_f32(capture_sec);
    let mut rx = Vec::<f32>::new();
    let mut next_level_log = Instant::now() + Duration::from_millis(250);
    let mut showed_progress = false;
    while Instant::now() < deadline {
        while let Some(s) = in_cons.try_pop() {
            rx.push(s);
        }
        if opts.verbose && Instant::now() >= next_level_log {
            let tail_len = ((0.25 * cfg_rx.fs).round() as usize).max(1).min(rx.len());
            let tail = if tail_len > 0 {
                &rx[rx.len().saturating_sub(tail_len)..]
            } else {
                &rx[..0]
            };
            let (rms, peak, _) = signal_diag(tail);
            eprint!(
                "\rCapture level: samples={} tail={:.3}s rms={rms:.5} peak={peak:.5}",
                rx.len(),
                tail_len as f32 / cfg_rx.fs
            );
            let _ = std::io::stderr().flush();
            showed_progress = true;
            next_level_log += Duration::from_millis(250);
        }
        if let Some(rem) = deadline.checked_duration_since(Instant::now()) {
            std::thread::sleep(rem.min(Duration::from_millis(10)));
        } else {
            break;
        }
    }
    if showed_progress {
        eprint!("\r{: <96}\r", "");
        let _ = std::io::stderr().flush();
    }
    drop(out_stream);
    drop(in_stream);

    let target_samples = (capture_sec * cfg_rx.fs).round().max(0.0) as usize;
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
    if let Some(path) = &opts.dump_tx_wav {
        save_wav_mono_i16(
            Path::new(path),
            &tx_plan.scheduled,
            tx_plan.fs.round() as u32,
        )?;
        info_line!("Saved TX schedule: {path}");
    }
    if let Some(path) = &opts.dump_wav {
        save_wav_mono_i16(Path::new(path), &rx, cfg_rx.fs.round() as u32)?;
        info_line!("Saved RX capture: {path}");
    }
    if opts.spectrogram {
        let spec_path = Path::new(&opts.spectrogram_path);
        let default_spec = acoustic_ofdm::SpectrogramOptions::default();
        if opts.spectrogram_opts.nfft == default_spec.nfft
            && opts.spectrogram_opts.hop == default_spec.hop
            && opts.spectrogram_opts.window == default_spec.window
        {
            save_spectrogram_png(spec_path, &rx, cfg_rx.fs)?;
        } else {
            save_spectrogram_png_with_options(spec_path, &rx, cfg_rx.fs, opts.spectrogram_opts)?;
        }
        info_line!("Saved spectrogram PNG: {}", spec_path.display());
    }

    if opts.verbose {
        let (rms, peak, first_loud) = signal_diag(&rx);
        let (clipped, clipped_frac) = clipping_diag(&rx);
        debug_line!(
            "RX diagnostics: rms={rms:.5} peak={peak:.5} first_loud_sample={first_loud} ({:.3}s)",
            (first_loud as f32) / cfg_rx.fs
        );
        debug_line!(
            "RX clipping: {clipped} samples ({:.2}%) at |x| >= 0.995",
            100.0 * clipped_frac
        );
        if peak < 0.01 {
            warn_line!("RX warning: very low capture level; increase speaker volume or mic gain.");
        }
    }

    let packet_sec = tx_plan.packet.len() as f32 / tx_plan.fs;
    let window_sec = packet_sec + 0.10;
    let sync_min = 0usize;
    let sync_max = 48usize;
    let sync_step = 2usize;

    #[derive(Clone)]
    struct Candidate {
        burst_idx: usize,
        nominal_start_sec: f32,
        start_sec: f32,
        sync_off: usize,
        sync_rms: f32,
        sync_peak: f32,
        post_rms: f32,
        post_peak: f32,
        train_rms: f32,
        train_evm: f32,
        pilot_evm: f32,
        data_evm: f32,
        decoded: bool,
        payload: Option<Vec<u8>>,
    }

    let mut best_overall: Option<Candidate> = None;
    for (i, &burst_start_tx_sec) in tx_plan.burst_starts_sec.iter().enumerate() {
        let nominal_start_sec = tx_start_in_capture + burst_start_tx_sec;
        let mut best_local: Option<Candidate> = None;
        let start_sec = nominal_start_sec;
        let start = (start_sec * cfg_rx.fs).round() as usize;
        let win_len = (window_sec * cfg_rx.fs).round() as usize;
        if start >= rx.len() || start + win_len > rx.len() {
            continue;
        }
        let window = &rx[start..start + win_len];
        let sync_candidates: Vec<usize> = if opts.oracle {
            vec![0]
        } else {
            (sync_min..=sync_max).step_by(sync_step).collect()
        };
        for sync_off in sync_candidates {
            let payload = decode_single_packet_passband_with_sync(window, &cfg_rx, sync_off as f32);
            let diag = diagnose_passband_window_with_sync(window, &cfg_rx, sync_off as f32);
            let decoded = payload.is_some();
            let cand = Candidate {
                burst_idx: i + 1,
                nominal_start_sec,
                start_sec,
                sync_off,
                sync_rms: diag.sync_rms,
                sync_peak: diag.sync_peak,
                post_rms: diag.post_rms,
                post_peak: diag.post_peak,
                train_rms: diag.train_rms,
                train_evm: diag.train_recon_evm,
                pilot_evm: diag.pilot_residual_evm,
                data_evm: diag.post_eq_evm,
                decoded,
                payload,
            };
            let oracle_match = opts.oracle
                && cand
                    .payload
                    .as_deref()
                    .map(|p| p == ORACLE_PAYLOAD)
                    .unwrap_or(false);
            let better = match &best_local {
                Some(best) => {
                    cand.decoded.cmp(&best.decoded).is_gt()
                        || (cand.decoded == best.decoded
                            && cand
                                .data_evm
                                .total_cmp(&best.data_evm)
                                .then_with(|| best.train_rms.total_cmp(&cand.train_rms))
                                .is_lt())
                }
                None => true,
            };
            if better || oracle_match {
                best_local = Some(cand.clone());
            }
            if oracle_match {
                break;
            }
        }
        if let Some(best) = best_local {
            info_line!(
                "Burst {}: nominal={:.3}s best_start={:.3}s sync_off={} level[sync/post]=[peak {:.3}/{:.3}, rms {:.3}/{:.3}] evm[train/pilot/data]=[{:.3}/{:.3}/{:.3}] decoded={}",
                best.burst_idx,
                best.nominal_start_sec,
                best.start_sec,
                best.sync_off,
                best.sync_peak,
                best.post_peak,
                best.sync_rms,
                best.post_rms,
                best.train_evm,
                best.pilot_evm,
                best.data_evm,
                best.decoded
            );
            let better = match &best_overall {
                Some(global) => {
                    best.decoded.cmp(&global.decoded).is_gt()
                        || (best.decoded == global.decoded
                            && best
                                .data_evm
                                .total_cmp(&global.data_evm)
                                .then_with(|| global.train_rms.total_cmp(&best.train_rms))
                                .is_lt())
                }
                None => true,
            };
            if better {
                best_overall = Some(best.clone());
            }
            if opts.oracle
                && best
                    .payload
                    .as_deref()
                    .map(|p| p == ORACLE_PAYLOAD)
                    .unwrap_or(false)
            {
                best_overall = Some(best);
                break;
            }
        }
    }

    let Some(best) = best_overall else {
        return Err("no valid burst window found in roundtrip capture".into());
    };

    let start = (best.start_sec * cfg_rx.fs).round() as usize;
    let win_len = (window_sec * cfg_rx.fs).round() as usize;
    let window = &rx[start..start + win_len];
    if let Some(path) = &opts.dump_wav {
        if let Some(dir) = Path::new(path).parent() {
            if let Some(cd) = dump_passband_constellation(window, &cfg_rx, best.sync_off as f32) {
                let png = dir.join("ofdm_constellation.png");
                let pre_csv = dir.join("ofdm_constellation_pre_eq.csv");
                let post_csv = dir.join("ofdm_constellation_post_eq.csv");
                save_constellation_comparison_png(&png, &cd.pre_eq, &cd.post_eq)?;
                {
                    let mut file = std::fs::File::create(&pre_csv)?;
                    writeln!(file, "re,im")?;
                    for z in &cd.pre_eq {
                        writeln!(file, "{},{}", z.re, z.im)?;
                    }
                }
                {
                    let mut file = std::fs::File::create(&post_csv)?;
                    writeln!(file, "re,im")?;
                    for z in &cd.post_eq {
                        writeln!(file, "{},{}", z.re, z.im)?;
                    }
                }
                info_line!("Saved constellation PNG: {}", png.display());
                info_line!("Saved constellation CSV: {}", pre_csv.display());
                info_line!("Saved constellation CSV: {}", post_csv.display());
            }
            if opts.oracle {
                if let Some(ch) = dump_passband_channel_compare_with_sync(
                    payload,
                    window,
                    &cfg_rx,
                    best.sync_off as f32,
                ) {
                    let csv_path = dir.join("ofdm_channel_compare.csv");
                    let png_path = dir.join("ofdm_channel_compare.png");
                    save_channel_compare_csv(&csv_path, &ch)?;
                    save_channel_compare_png(&png_path, &ch)?;
                    info_line!("Saved channel compare CSV: {}", csv_path.display());
                    info_line!("Saved channel compare PNG: {}", png_path.display());
                }
            }
            if let Some(chain) = dump_passband_iq_chain(window, &cfg_rx) {
                let audio_prefix = dir.join("iq_down_audio_rate");
                let bb_prefix = dir.join("iq_down_baseband_rate");
                save_complex_parts_wav(
                    &audio_prefix,
                    &chain.downconverted_audio_rate,
                    chain.fs_audio.round() as u32,
                )?;
                save_complex_parts_wav(
                    &bb_prefix,
                    &chain.baseband_rate,
                    chain.fs_baseband.round() as u32,
                )?;
                info_line!(
                    "Saved IQ downconverted WAVs: {}.(re|im).wav",
                    audio_prefix.display()
                );
                info_line!(
                    "Saved IQ baseband WAVs: {}.(re|im).wav",
                    bb_prefix.display()
                );
            }
        }
    }
    if let Some(dump) = dump_passband_bins_with_sync(window, &cfg_rx, best.sync_off as f32) {
        let path = Path::new("/tmp/ofdm_decode_bins.csv");
        save_bin_dump_csv(path, &dump)?;
        info_line!("Saved bin dump CSV: {}", path.display());
    }

    if let Some(payload) = best.payload {
        if opts.oracle {
            info_line!(
                "Oracle verdict: {}",
                if payload.as_slice() == ORACLE_PAYLOAD {
                    "match"
                } else {
                    "mismatch"
                }
            );
        }
        if stdout_raw {
            let mut out = std::io::stdout().lock();
            out.write_all(&payload)?;
            out.flush()?;
        } else {
            let hex = payload
                .iter()
                .map(|b| format!("{b:02X}"))
                .collect::<Vec<_>>()
                .join(" ");
            info_line!(
                "Decoded {} bytes from burst={} start={:.3}s sync_off={}",
                payload.len(),
                best.burst_idx,
                best.start_sec,
                best.sync_off
            );
            info_line!("HEX: {hex}");
            info_line!("UTF8(lossy): {}", String::from_utf8_lossy(&payload));
        }
        Ok(())
    } else {
        if let Some(path) = &opts.dump_wav {
            if let Some(dir) = Path::new(path).parent() {
                if let Some(raw) = recover_decided_packet_bytes_passband_with_sync(
                    window,
                    &cfg_rx,
                    best.sync_off as f32,
                ) {
                    let raw_path = dir.join("ofdm_pre_crc_bytes.bin");
                    std::fs::write(&raw_path, &raw)?;
                    let inspect = inspect_packet_bytes(&raw);
                    info_line!("Saved pre-CRC bytes: {}", raw_path.display());
                    info_line!(
                        "Pre-CRC packet inspect: preamble_ok={} header_ok={} payload_len={:?} total_len={:?} enough_total={} crc_ok={}",
                        inspect.preamble_ok,
                        inspect.enough_for_header,
                        inspect.payload_len,
                        inspect.total_len,
                        inspect.enough_for_total,
                        inspect.crc_ok
                    );
                }
            }
        }
        Err(format!(
            "roundtrip decode failed: best burst={} start={:.3}s sync_off={} evm[train/pilot/data]=[{:.3}/{:.3}/{:.3}]",
            best.burst_idx,
            best.start_sec,
            best.sync_off,
            best.train_evm,
            best.pilot_evm,
            best.data_evm
        )
        .into())
    }
}
