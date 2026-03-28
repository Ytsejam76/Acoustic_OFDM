// Copyright (c) 2026 Elias S. G. Carotti

use std::error::Error;
use std::io::Write;
use std::path::Path;

mod audio;
mod cli_args;
mod live_profile;
mod logging;
mod mic_roundtrip;
mod rx;
mod tx;

use acoustic_ofdm::{
    decode_single_packet_passband, decode_single_packet_passband_with_sync,
    diagnose_passband_window, diagnose_passband_window_with_sync, dump_passband_bins,
    dump_passband_bins_with_sync, encode_single_packet_passband,
    encode_single_packet_passband_body, load_wav_mono_f32, save_spectrogram_png,
    save_spectrogram_png_with_options, save_wav_mono_i16, OfdmConfig, PassbandBinDump,
    SpectrogramOptions,
};
use clap::Parser;
use live_profile::{rx_default_log_file, tx_default_log_file};

use crate::cli_args::*;
use crate::logging::init_logging;
use crate::mic_roundtrip::cmd_mic_roundtrip;
use crate::rx::cmd_rx;
use crate::tx::cmd_tx;

fn cmd_encode(out_path: &Path, payload: &[u8], cfg: &OfdmConfig) -> Result<(), Box<dyn Error>> {
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

fn cmd_encode_body(out_path: &Path, payload: &[u8], cfg: &OfdmConfig) -> Result<(), Box<dyn Error>> {
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

fn cmd_decode(in_path: &Path, cfg: &OfdmConfig, stdout_raw: bool) -> Result<(), Box<dyn Error>> {
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
        println!("Decoded {} bytes", payload.len());
        println!("HEX: {hex}");
        println!("UTF8(lossy): {}", String::from_utf8_lossy(&payload));
    }
    Ok(())
}

fn cmd_decode_window(
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

fn cmd_spectrogram(cmd: SpectrogramCmd) -> Result<(), Box<dyn Error>> {
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

fn cmd_scan(cmd: ScanCmd) -> Result<(), Box<dyn Error>> {
    let (samples, sr) = load_wav_mono_f32(Path::new(&cmd.in_wav))?;
    let mut cfg = OfdmConfig::default();
    apply_common_cfg(&mut cfg, &cmd.common)?;
    let rx_profile = RxCmd {
        common: CommonCfgArgs::default(),
        profile: cmd.profile,
        duration_sec: None,
        mic_gain: None,
        in_hp_hz: 12_000.0,
        in_lp_hz: 19_000.0,
        no_input_filter: true,
        dump_wav: None,
        spectrogram: false,
        spectrogram_path: "/tmp/rx_spectrogram.png".to_string(),
        spectrogram_nfft: 512,
        spectrogram_hop: 128,
        spectrogram_window: SpectrogramWindowArg::Hann,
        oracle: false,
        stdout: false,
        verbose: false,
        log_level: None,
        log_file: None,
    };
    apply_rx_profile_cfg(&mut cfg, &rx_profile);
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
            let cand = Candidate {
                start_sec: pos as f32 / cfg.fs,
                sync_off,
                train_rms: diag.train_rms,
                train_evm: diag.train_recon_evm,
                pilot_evm: diag.pilot_residual_evm,
                data_evm: diag.post_eq_evm,
                decoded: diag.decoded,
                payload_len: diag.decoded_payload_len,
            };
            best.push(cand);
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

    let total = best.len();
    println!(
        "Scanned wav={} sr={}Hz window={:.3}s step={:.1}ms sync=[{}..{}] step={} total={}",
        cmd.in_wav,
        sr,
        cmd.window_sec,
        cmd.step_ms,
        cmd.sync_min,
        cmd.sync_max,
        cmd.sync_step,
        total
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

fn main() -> Result<(), Box<dyn Error>> {
    let cli = Cli::parse();

    match cli.command {
        Commands::Encode(cmd) => {
            let mut cfg = OfdmConfig::default();
            apply_common_cfg(&mut cfg, &cmd.common)?;
            let payload = payload_from_arg_or_stdin(&cmd.payload_text)?;
            cmd_encode(Path::new(&cmd.out_wav), &payload, &cfg)?;
        }
        Commands::EncodeBody(cmd) => {
            let mut cfg = OfdmConfig::default();
            apply_common_cfg(&mut cfg, &cmd.common)?;
            let payload = payload_from_arg_or_stdin(&cmd.payload_text)?;
            cmd_encode_body(Path::new(&cmd.out_wav), &payload, &cfg)?;
        }
        Commands::Decode(cmd) => {
            let mut cfg = OfdmConfig::default();
            apply_common_cfg(&mut cfg, &cmd.common)?;
            let rx_profile = RxCmd {
                common: CommonCfgArgs::default(),
                profile: cmd.profile,
                duration_sec: None,
                mic_gain: None,
                in_hp_hz: 12_000.0,
                in_lp_hz: 19_000.0,
                no_input_filter: true,
                dump_wav: None,
                spectrogram: false,
                spectrogram_path: "/tmp/rx_spectrogram.png".to_string(),
                spectrogram_nfft: 512,
                spectrogram_hop: 128,
                spectrogram_window: SpectrogramWindowArg::Hann,
                oracle: false,
                stdout: false,
                verbose: false,
                log_level: None,
                log_file: None,
            };
            apply_rx_profile_cfg(&mut cfg, &rx_profile);
            if let Some(start_sec) = cmd.start_sec {
                cmd_decode_window(
                    Path::new(&cmd.in_wav),
                    &cfg,
                    start_sec,
                    cmd.window_sec,
                    cmd.sync_off,
                    cmd.stdout,
                )?;
            } else {
                cmd_decode(Path::new(&cmd.in_wav), &cfg, cmd.stdout)?;
            }
        }
        Commands::MicRoundtrip(cmd) => {
            let mut cfg = OfdmConfig::default();
            apply_common_cfg(&mut cfg, &cmd.common)?;
            apply_mic_roundtrip_profile_cfg(&mut cfg, &cmd);
            let log_path = cmd
                .log_file
                .clone()
                .or_else(|| Some("acoustic_ofdm_mic_roundtrip.log".to_string()));
            init_logging(
                log_path.as_deref(),
                resolved_log_level(cmd.log_level, cmd.verbose),
            )?;
            if let Some(path) = &log_path {
                info_line!("Log file: {path}");
            }
            let opts = mic_roundtrip_audio_opts(&cmd);
            let payload = if cmd.oracle {
                ORACLE_PAYLOAD.to_vec()
            } else {
                let arg = cmd
                    .payload_text
                    .as_deref()
                    .ok_or("mic-roundtrip requires <payload_text|stdin>")?;
                payload_from_arg_or_stdin(arg)?
            };
            if payload.is_empty() {
                return Err("payload must not be empty".into());
            }
            cmd_mic_roundtrip(&payload, &cfg, &opts, cmd.stdout)?;
        }
        Commands::Roundtrip(cmd) => {
            let mut cfg = OfdmConfig::default();
            apply_common_cfg(&mut cfg, &cmd.common)?;
            let payload = payload_from_arg_or_stdin(&cmd.payload_text)?;
            cmd_encode(Path::new(&cmd.wav_path), &payload, &cfg)?;
            cmd_decode(Path::new(&cmd.wav_path), &cfg, cmd.stdout)?;
        }
        Commands::Scan(cmd) => {
            cmd_scan(cmd)?;
        }
        Commands::Spectrogram(cmd) => {
            cmd_spectrogram(cmd)?;
        }
        Commands::Rx(cmd) => {
            let mut cfg = OfdmConfig::default();
            apply_common_cfg(&mut cfg, &cmd.common)?;
            apply_rx_profile_cfg(&mut cfg, &cmd);
            let log_path = cmd
                .log_file
                .clone()
                .or_else(|| rx_default_log_file(cmd.profile));
            init_logging(
                log_path.as_deref(),
                resolved_log_level(cmd.log_level, cmd.verbose),
            )?;
            if let Some(path) = &log_path {
                info_line!("Log file: {path}");
            }
            let opts = rx_audio_opts(&cmd);
            cmd_rx(&cfg, &opts, cmd.stdout)?;
        }
        Commands::Tx(cmd) => {
            let mut cfg = OfdmConfig::default();
            apply_common_cfg(&mut cfg, &cmd.common)?;
            apply_tx_profile_cfg(&mut cfg, &cmd);
            let log_path = cmd
                .log_file
                .clone()
                .or_else(|| tx_default_log_file(cmd.profile));
            init_logging(
                log_path.as_deref(),
                resolved_log_level(cmd.log_level, cmd.verbose),
            )?;
            if let Some(path) = &log_path {
                info_line!("Log file: {path}");
            }
            let opts = tx_audio_opts(&cmd);
            let payload = if cmd.oracle {
                ORACLE_PAYLOAD.to_vec()
            } else {
                let arg = cmd
                    .payload_text
                    .as_deref()
                    .ok_or("tx requires <payload_text|stdin>")?;
                payload_from_arg_or_stdin(arg)?
            };
            if payload.is_empty() {
                return Err("payload must not be empty".into());
            }
            cmd_tx(&payload, &cfg, &opts)?;
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::{
        apply_mic_roundtrip_profile_cfg, apply_rx_profile_cfg, apply_tx_profile_cfg,
        mic_roundtrip_audio_opts, payload_from_arg, rx_audio_opts, tx_audio_opts, Cli, Commands,
        WakePreambleArg, ORACLE_PAYLOAD,
    };
    use crate::live_profile::LiveProfileArg;
    use crate::live_profile::{rx_default_log_file, tx_default_log_file};
    use acoustic_ofdm::{OfdmConfig, WakePreamble};
    use clap::Parser;

    #[test]
    fn payload_from_arg_utf8_bytes() {
        let p = payload_from_arg("ciao-OFDM");
        assert_eq!(p, b"ciao-OFDM");
    }

    #[test]
    fn parse_encode_common_options() {
        let cli = Cli::try_parse_from([
            "acoustic_ofdm_cli",
            "encode",
            "--base-freq-hz",
            "2500",
            "--wake-preamble",
            "gold",
            "out.wav",
            "hello",
        ])
        .expect("parse failed");
        match cli.command {
            Commands::Encode(cmd) => {
                assert_eq!(cmd.common.base_freq_hz, Some(2500.0));
                assert_eq!(cmd.common.wake_preamble, Some(WakePreambleArg::Gold));
                assert_eq!(cmd.out_wav, "out.wav");
                assert_eq!(cmd.payload_text, "hello");
            }
            _ => panic!("expected encode"),
        }
    }

    #[test]
    fn parse_decode_stdout_override() {
        let cli = Cli::try_parse_from([
            "acoustic_ofdm_cli",
            "decode",
            "--stdout",
            "--profile",
            "live-debug",
            "--start-sec",
            "1.25",
            "--window-sec",
            "1.40",
            "in.wav",
        ])
        .expect("parse failed");
        match cli.command {
            Commands::Decode(cmd) => {
                assert!(cmd.stdout);
                assert_eq!(cmd.profile, crate::live_profile::LiveProfileArg::LiveDebug);
                assert_eq!(cmd.start_sec, Some(1.25));
                assert_eq!(cmd.window_sec, Some(1.40));
                assert_eq!(cmd.sync_off, None);
                assert_eq!(cmd.in_wav, "in.wav");
            }
            _ => panic!("expected decode"),
        }
    }

    #[test]
    fn parse_spectrogram_options() {
        let cli = Cli::try_parse_from([
            "acoustic_ofdm_cli",
            "spectrogram",
            "--in-wav",
            "/tmp/in.wav",
            "--out-png",
            "/tmp/out.png",
            "--spectrogram-nfft",
            "1024",
            "--spectrogram-hop",
            "256",
            "--spectrogram-window",
            "blackman",
        ])
        .expect("parse failed");
        match cli.command {
            Commands::Spectrogram(cmd) => {
                assert_eq!(cmd.in_wav, "/tmp/in.wav");
                assert_eq!(cmd.out_png, "/tmp/out.png");
                assert_eq!(cmd.spectrogram_nfft, 1024);
                assert_eq!(cmd.spectrogram_hop, 256);
                assert_eq!(cmd.spectrogram_window, crate::cli_args::SpectrogramWindowArg::Blackman);
            }
            _ => panic!("expected spectrogram"),
        }
    }

    #[test]
    fn parse_rx_options() {
        let cli = Cli::try_parse_from([
            "acoustic_ofdm_cli",
            "rx",
            "--profile",
            "standard",
            "--duration-sec",
            "2.5",
            "--mic-gain",
            "0.8",
            "--wake-preamble",
            "tone",
            "--dump-wav",
            "/tmp/rx.wav",
            "--spectrogram",
            "--spectrogram-path",
            "/tmp/rx_spec.png",
            "--verbose",
        ])
        .expect("parse failed");
        match cli.command {
            Commands::Rx(cmd) => {
                assert_eq!(cmd.profile, LiveProfileArg::Standard);
                assert_eq!(cmd.duration_sec, Some(2.5));
                assert_eq!(cmd.mic_gain, Some(0.8));
                assert_eq!(cmd.common.wake_preamble, Some(WakePreambleArg::Tone));
                assert_eq!(cmd.dump_wav.as_deref(), Some("/tmp/rx.wav"));
                assert!(cmd.spectrogram);
                assert_eq!(cmd.spectrogram_path, "/tmp/rx_spec.png");
                assert!(cmd.verbose);
                assert!(!cmd.stdout);
            }
            _ => panic!("expected rx"),
        }
    }

    #[test]
    fn parse_tx_options() {
        let cli = Cli::try_parse_from([
            "acoustic_ofdm_cli",
            "tx",
            "--profile",
            "standard",
            "--spk-gain",
            "0.7",
            "--repeats",
            "4",
            "--wake-preamble",
            "gold",
            "--verbose",
            "hello",
        ])
        .expect("parse failed");
        match cli.command {
            Commands::Tx(cmd) => {
                assert_eq!(cmd.profile, LiveProfileArg::Standard);
                assert_eq!(cmd.spk_gain, Some(0.7));
                assert_eq!(cmd.repeats, Some(4));
                assert_eq!(cmd.common.wake_preamble, Some(WakePreambleArg::Gold));
                assert!(cmd.verbose);
                assert_eq!(cmd.payload_text.as_deref(), Some("hello"));
            }
            _ => panic!("expected tx"),
        }
    }

    #[test]
    fn parse_tx_oracle_options() {
        let cli = Cli::try_parse_from(["acoustic_ofdm_cli", "tx", "--oracle", "--repeats", "2"])
            .expect("parse failed");
        match cli.command {
            Commands::Tx(cmd) => {
                assert!(cmd.oracle);
                assert_eq!(cmd.repeats, Some(2));
                assert!(cmd.payload_text.is_none());
                assert_eq!(ORACLE_PAYLOAD, b"ACOUSTIC-OFDM-ORACLE");
            }
            _ => panic!("expected tx"),
        }
    }

    #[test]
    fn parse_mic_roundtrip_options() {
        let cli = Cli::try_parse_from([
            "acoustic_ofdm_cli",
            "mic-roundtrip",
            "--profile",
            "live-debug",
            "--mic-gain",
            "0.8",
            "--spk-gain",
            "0.7",
            "--repeats",
            "4",
            "--oracle",
        ])
        .expect("parse failed");
        match cli.command {
            Commands::MicRoundtrip(cmd) => {
                assert_eq!(cmd.profile, LiveProfileArg::LiveDebug);
                assert_eq!(cmd.mic_gain, Some(0.8));
                assert_eq!(cmd.spk_gain, Some(0.7));
                assert_eq!(cmd.repeats, Some(4));
                assert!(cmd.oracle);
            }
            _ => panic!("expected mic-roundtrip"),
        }
    }

    #[test]
    fn tx_live_debug_profile_defaults() {
        let cli = Cli::try_parse_from(["acoustic_ofdm_cli", "tx", "hello"]).expect("parse failed");
        match cli.command {
            Commands::Tx(cmd) => {
                assert_eq!(cmd.profile, LiveProfileArg::LiveDebug);
                let mut cfg = OfdmConfig::default();
                apply_tx_profile_cfg(&mut cfg, &cmd);
                let opts = tx_audio_opts(&cmd);
                assert_eq!(cfg.wake_preamble, WakePreamble::Tone);
                assert!((opts.spk_gain - 0.2).abs() < 1e-6);
                assert!((opts.pre_delay_sec - 0.5).abs() < 1e-6);
                assert_eq!(opts.repeats, 5);
                assert!(opts.oracle);
                assert_eq!(
                    tx_default_log_file(cmd.profile).as_deref(),
                    Some("acoustic_ofdm_tx.log")
                );
            }
            _ => panic!("expected tx"),
        }
    }

    #[test]
    fn rx_live_debug_profile_defaults() {
        let cli = Cli::try_parse_from(["acoustic_ofdm_cli", "rx"]).expect("parse failed");
        match cli.command {
            Commands::Rx(cmd) => {
                assert_eq!(cmd.profile, LiveProfileArg::LiveDebug);
                let mut cfg = OfdmConfig::default();
                apply_rx_profile_cfg(&mut cfg, &cmd);
                let opts = rx_audio_opts(&cmd);
                assert_eq!(cfg.wake_preamble, WakePreamble::Tone);
                assert!((opts.duration_sec - 5.0).abs() < 1e-6);
                assert!((opts.mic_gain - 0.2).abs() < 1e-6);
                assert_eq!(opts.dump_wav.as_deref(), Some("/tmp/rx_capture.wav"));
                assert!(opts.spectrogram);
                assert!(opts.oracle);
                assert!(opts.verbose);
                assert_eq!(
                    rx_default_log_file(cmd.profile).as_deref(),
                    Some("acoustic_ofdm_rx.log")
                );
            }
            _ => panic!("expected rx"),
        }
    }

    #[test]
    fn mic_roundtrip_live_debug_profile_defaults() {
        let cli =
            Cli::try_parse_from(["acoustic_ofdm_cli", "mic-roundtrip", "--oracle"]).expect("parse failed");
        match cli.command {
            Commands::MicRoundtrip(cmd) => {
                assert_eq!(cmd.profile, LiveProfileArg::LiveDebug);
                let mut cfg = OfdmConfig::default();
                apply_mic_roundtrip_profile_cfg(&mut cfg, &cmd);
                let opts = mic_roundtrip_audio_opts(&cmd);
                assert_eq!(cfg.wake_preamble, WakePreamble::Tone);
                assert!((opts.duration_sec - 10.0).abs() < 1e-6);
                assert!((opts.mic_gain - 0.2).abs() < 1e-6);
                assert!((opts.spk_gain - 0.2).abs() < 1e-6);
                assert_eq!(opts.repeats, 5);
                assert!(opts.oracle);
            }
            _ => panic!("expected mic-roundtrip"),
        }
    }
}

// vim: set ts=4 sw=4 et:
