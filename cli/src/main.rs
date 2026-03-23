// Copyright (c) 2026 Elias S. G. Carotti

use std::error::Error;
use std::io::Write;
use std::path::Path;

mod audio;
mod cli_args;
mod live_profile;
mod logging;
mod rx;
mod rx_support;
mod tx;

use acoustic_ofdm::{
    decode_single_packet_passband, encode_single_packet_passband, load_wav_mono_f32,
    save_wav_mono_i16, OfdmConfig,
};
use clap::Parser;
use live_profile::{rx_default_log_file, tx_default_log_file};

use crate::cli_args::*;
use crate::logging::init_logging;
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

fn cmd_decode(in_path: &Path, cfg: &OfdmConfig, stdout_raw: bool) -> Result<(), Box<dyn Error>> {
    let (samples, sr) = load_wav_mono_f32(in_path)?;
    if sr != cfg.fs as u32 {
        return Err(format!("sample-rate mismatch: wav={sr} cfg={}", cfg.fs as u32).into());
    }
    let payload = decode_single_packet_passband(&samples, cfg).ok_or("decode failed")?;
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

fn main() -> Result<(), Box<dyn Error>> {
    let cli = Cli::parse();

    match cli.command {
        Commands::Encode(cmd) => {
            let mut cfg = OfdmConfig::default();
            apply_common_cfg(&mut cfg, &cmd.common)?;
            let payload = payload_from_arg_or_stdin(&cmd.payload_text)?;
            cmd_encode(Path::new(&cmd.out_wav), &payload, &cfg)?;
        }
        Commands::Decode(cmd) => {
            let mut cfg = OfdmConfig::default();
            apply_common_cfg(&mut cfg, &cmd.common)?;
            cmd_decode(Path::new(&cmd.in_wav), &cfg, cmd.stdout)?;
        }
        Commands::Roundtrip(cmd) => {
            let mut cfg = OfdmConfig::default();
            apply_common_cfg(&mut cfg, &cmd.common)?;
            let payload = payload_from_arg_or_stdin(&cmd.payload_text)?;
            cmd_encode(Path::new(&cmd.wav_path), &payload, &cfg)?;
            cmd_decode(Path::new(&cmd.wav_path), &cfg, cmd.stdout)?;
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
        apply_rx_profile_cfg, apply_tx_profile_cfg, payload_from_arg, rx_audio_opts, tx_audio_opts,
        Cli, Commands, WakePreambleArg, ORACLE_PAYLOAD,
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
        let cli = Cli::try_parse_from(["acoustic_ofdm_cli", "decode", "--stdout", "in.wav"])
            .expect("parse failed");
        match cli.command {
            Commands::Decode(cmd) => {
                assert!(cmd.stdout);
                assert_eq!(cmd.in_wav, "in.wav");
            }
            _ => panic!("expected decode"),
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
    fn tx_live_debug_profile_defaults() {
        let cli = Cli::try_parse_from(["acoustic_ofdm_cli", "tx", "hello"]).expect("parse failed");
        match cli.command {
            Commands::Tx(cmd) => {
                assert_eq!(cmd.profile, LiveProfileArg::LiveDebug);
                let mut cfg = OfdmConfig::default();
                apply_tx_profile_cfg(&mut cfg, &cmd);
                let opts = tx_audio_opts(&cmd);
                assert_eq!(cfg.wake_preamble, WakePreamble::Gold);
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
                assert_eq!(cfg.wake_preamble, WakePreamble::Gold);
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

}

// vim: set ts=4 sw=4 et:
