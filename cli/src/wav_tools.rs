// Copyright (c) 2026 Elias S. G. Carotti

use std::error::Error;
use std::io::Write;
use std::path::Path;

use acoustic_ofdm::{
    decode_single_packet_passband, decode_single_packet_passband_with_sync,
    diagnose_passband_window, diagnose_passband_window_with_sync, dump_passband_bins,
    dump_passband_bins_with_sync, encode_single_packet_passband,
    encode_single_packet_passband_body, load_wav_mono_f32, save_spectrogram_png,
    save_spectrogram_png_with_options, save_wav_mono_i16, OfdmConfig, PassbandBinDump,
    SpectrogramOptions,
};

use crate::cli_args::{
    apply_common_cfg, apply_profile_cfg, CommonCfgArgs, DecodeCmd, ScanCmd, SpectrogramCmd,
};

fn save_bin_dump_csv(path: &str, dump: &PassbandBinDump) -> Result<(), Box<dyn Error>> {
    let mut file = std::fs::File::create(path)?;
    writeln!(file, "data_symbol_idx,used_bin,role,pre_re,pre_im,post_re,post_im,ref_re,ref_im")?;
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

fn print_decode_payload(
    payload: &[u8],
    stdout_raw: bool,
    window_info: Option<(f32, f32)>,
) -> Result<(), Box<dyn Error>> {
    if stdout_raw {
        let mut out = std::io::stdout().lock();
        out.write_all(payload)?;
        out.flush()?;
        return Ok(());
    }
    let hex = payload
        .iter()
        .map(|b| format!("{b:02X}"))
        .collect::<Vec<_>>()
        .join(" ");
    match window_info {
        Some((start_sec, window_sec)) => {
            println!(
                "Decoded {} bytes from start={:.3}s window={:.3}s",
                payload.len(),
                start_sec,
                window_sec
            );
        }
        None => println!("Decoded {} bytes", payload.len()),
    }
    println!("HEX: {hex}");
    println!("UTF8(lossy): {}", String::from_utf8_lossy(payload));
    Ok(())
}

fn cfg_for_decode_profile(
    common: &CommonCfgArgs,
    profile: crate::live_profile::LiveProfileArg,
) -> Result<OfdmConfig, Box<dyn Error>> {
    let mut cfg = OfdmConfig::default();
    apply_common_cfg(&mut cfg, common)?;
    apply_profile_cfg(&mut cfg, profile, common.wake_preamble.is_none());
    Ok(cfg)
}

pub(crate) fn cmd_encode(
    out_path: &Path,
    payload: &[u8],
    cfg: &OfdmConfig,
) -> Result<(), Box<dyn Error>> {
    if payload.len() > cfg.packet_payload_bytes {
        return Err(format!(
            "payload too long for single-packet app: {} > {}",
            payload.len(),
            cfg.packet_payload_bytes
        )
        .into());
    }
    let y = encode_single_packet_passband(payload, cfg);
    save_wav_mono_i16(out_path, &y, cfg.fs as u32)?;
    println!("Wrote WAV: {}", out_path.display());
    Ok(())
}

pub(crate) fn cmd_encode_body(
    out_path: &Path,
    payload: &[u8],
    cfg: &OfdmConfig,
) -> Result<(), Box<dyn Error>> {
    if payload.len() > cfg.packet_payload_bytes {
        return Err(format!(
            "payload too long for single-packet app: {} > {}",
            payload.len(),
            cfg.packet_payload_bytes
        )
        .into());
    }
    let y = encode_single_packet_passband_body(payload, cfg);
    save_wav_mono_i16(out_path, &y, cfg.fs as u32)?;
    println!("Wrote OFDM-body WAV: {}", out_path.display());
    Ok(())
}

pub(crate) fn cmd_decode(
    in_path: &Path,
    cfg: &OfdmConfig,
    stdout_raw: bool,
) -> Result<(), Box<dyn Error>> {
    let (samples, sr) = load_wav_mono_f32(in_path)?;
    let mut cfg = cfg.clone();
    cfg.fs = sr as f32;
    let payload = match decode_single_packet_passband(&samples, &cfg) {
        Some(payload) => payload,
        None => {
            let diag = diagnose_passband_window(&samples, &cfg);
            return Err(format!(
                "decode failed: sync_off={} train_rms={:.4} evm[train/pilot/data]=[{:.3}/{:.3}/{:.3}]",
                diag.sync_off,
                diag.train_rms,
                diag.train_recon_evm,
                diag.pilot_residual_evm,
                diag.post_eq_evm
            )
            .into());
        }
    };
    print_decode_payload(&payload, stdout_raw, None)
}

pub(crate) fn cmd_decode_window(
    in_path: &Path,
    cfg: &OfdmConfig,
    start_sec: f32,
    window_sec: Option<f32>,
    sync_off: Option<f32>,
    stdout_raw: bool,
) -> Result<(), Box<dyn Error>> {
    let (samples, sr) = load_wav_mono_f32(in_path)?;
    let mut cfg = cfg.clone();
    cfg.fs = sr as f32;
    let start = (start_sec.max(0.0) * cfg.fs).round() as usize;
    if start >= samples.len() {
        return Err("start offset is beyond end of wav".into());
    }
    let end = if let Some(sec) = window_sec {
        let len = (sec.max(0.0) * cfg.fs).round() as usize;
        start.saturating_add(len).min(samples.len())
    } else {
        samples.len()
    };
    if end <= start {
        return Err("window is empty".into());
    }
    let window = &samples[start..end];
    let payload = match sync_off {
        Some(sync_off) => decode_single_packet_passband_with_sync(window, &cfg, sync_off),
        None => decode_single_packet_passband(window, &cfg),
    };
    let payload = match payload {
        Some(payload) => payload,
        None => {
            let bin_dump = match sync_off {
                Some(v) => dump_passband_bins_with_sync(window, &cfg, v),
                None => dump_passband_bins(window, &cfg),
            };
            if let Some(dump) = bin_dump {
                let path = "/tmp/ofdm_decode_bins.csv";
                save_bin_dump_csv(path, &dump)?;
                eprintln!("Wrote decode bin dump: {path}");
            }
            let diag = match sync_off {
                Some(v) => diagnose_passband_window_with_sync(window, &cfg, v),
                None => diagnose_passband_window(window, &cfg),
            };
            return Err(format!(
                "decode failed: start={start_sec:.3}s window={:.3}s {} train_rms={:.4} evm[train/pilot/data]=[{:.3}/{:.3}/{:.3}]",
                (end - start) as f32 / cfg.fs,
                match sync_off {
                    Some(v) => format!("sync_off={v:.3}"),
                    None => format!("sync_off={}", diag.sync_off),
                },
                diag.train_rms,
                diag.train_recon_evm,
                diag.pilot_residual_evm,
                diag.post_eq_evm
            )
            .into());
        }
    };
    print_decode_payload(
        &payload,
        stdout_raw,
        Some((start_sec, (end - start) as f32 / cfg.fs)),
    )
}

pub(crate) fn cmd_spectrogram(cmd: SpectrogramCmd) -> Result<(), Box<dyn Error>> {
    let (samples, sr) = load_wav_mono_f32(Path::new(&cmd.in_wav))?;
    let opts = SpectrogramOptions {
        nfft: cmd.spectrogram_nfft,
        hop: cmd.spectrogram_hop,
        window: cmd.spectrogram_window.into(),
    };
    let default_opts = SpectrogramOptions::default();
    let out = Path::new(&cmd.out_png);
    if opts.nfft == default_opts.nfft
        && opts.hop == default_opts.hop
        && opts.window == default_opts.window
    {
        save_spectrogram_png(out, &samples, sr as f32)?;
    } else {
        save_spectrogram_png_with_options(out, &samples, sr as f32, opts)?;
    }
    println!("Wrote spectrogram: {}", out.display());
    Ok(())
}

pub(crate) fn cmd_scan(cmd: ScanCmd) -> Result<(), Box<dyn Error>> {
    let (samples, sr) = load_wav_mono_f32(Path::new(&cmd.in_wav))?;
    let mut cfg = cfg_for_decode_profile(&cmd.common, cmd.profile)?;
    cfg.fs = sr as f32;

    let wav_sec = samples.len() as f32 / cfg.fs;
    let start_sec = cmd.start_sec.max(0.0);
    let end_sec = cmd.end_sec.unwrap_or(wav_sec).min(wav_sec);
    let window_len = (cmd.window_sec.max(0.05) * cfg.fs).round() as usize;
    let step = ((cmd.step_ms.max(1.0) * 1e-3) * cfg.fs).round() as usize;
    let sync_step = cmd.sync_step.max(1) as usize;
    if end_sec <= start_sec {
        return Err("scan end must be after scan start".into());
    }
    if window_len == 0 || step == 0 {
        return Err("invalid scan parameters".into());
    }

    #[derive(Clone)]
    struct Candidate {
        start_sec: f32,
        sync_off: usize,
        train_rms: f32,
        train_evm: f32,
        pilot_evm: f32,
        data_evm: f32,
        decoded: bool,
        payload_len: Option<usize>,
    }

    let mut best: Vec<Candidate> = Vec::new();
    let start_idx = (start_sec * cfg.fs).round() as usize;
    let end_idx = (end_sec * cfg.fs).round() as usize;
    let mut pos = start_idx;
    while pos <= end_idx && pos + window_len <= samples.len() {
        let window = &samples[pos..pos + window_len];
        for sync_off in (cmd.sync_min.max(0) as usize..=cmd.sync_max.max(cmd.sync_min) as usize)
            .step_by(sync_step)
        {
            let diag = diagnose_passband_window_with_sync(window, &cfg, sync_off as f32);
            best.push(Candidate {
                start_sec: pos as f32 / cfg.fs,
                sync_off,
                train_rms: diag.train_rms,
                train_evm: diag.train_recon_evm,
                pilot_evm: diag.pilot_residual_evm,
                data_evm: diag.post_eq_evm,
                decoded: diag.decoded,
                payload_len: diag.decoded_payload_len,
            });
        }
        pos += step;
    }

    best.sort_by(|a, b| {
        a.data_evm
            .partial_cmp(&b.data_evm)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| {
                b.train_rms
                    .partial_cmp(&a.train_rms)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .then_with(|| {
                a.pilot_evm
                    .partial_cmp(&b.pilot_evm)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
    });

    println!(
        "Scanned wav={} sr={}Hz window={:.3}s step={:.1}ms sync=[{}..{}] step={} total={}",
        cmd.in_wav,
        sr,
        cmd.window_sec,
        cmd.step_ms,
        cmd.sync_min,
        cmd.sync_max,
        cmd.sync_step,
        best.len()
    );
    for (i, cand) in best.into_iter().take(cmd.top_k).enumerate() {
        println!(
            "#{} start={:.3}s sync_off={} train_rms={:.4} evm[train/pilot/data]=[{:.3}/{:.3}/{:.3}] decoded={} payload_len={}",
            i + 1,
            cand.start_sec,
            cand.sync_off,
            cand.train_rms,
            cand.train_evm,
            cand.pilot_evm,
            cand.data_evm,
            cand.decoded,
            cand.payload_len
                .map(|n| n.to_string())
                .unwrap_or_else(|| "-".to_string())
        );
    }
    Ok(())
}

pub(crate) fn decode_cfg_for_cmd(cmd: &DecodeCmd) -> Result<OfdmConfig, Box<dyn Error>> {
    cfg_for_decode_profile(&cmd.common, cmd.profile)
}
