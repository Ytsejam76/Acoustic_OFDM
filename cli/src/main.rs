// Copyright (c) 2026 Elias S. G. Carotti

use std::error::Error;
use std::path::Path;

use acoustic_ofdm::{
    decode_single_packet_passband,
    encode_single_packet_passband,
    load_wav_mono_f32,
    save_wav_mono_i16,
    OfdmConfig,
};

/// Parses payload text into bytes.
///
/// Parameters:
/// - `s`: payload string.
/// Returns:
/// - `Vec<u8>`: UTF-8 bytes.
fn payload_from_arg(s: &str) -> Vec<u8> {
    s.as_bytes().to_vec()
}

/// Prints command usage and exits with failure.
///
/// Parameters:
/// - none.
/// Returns:
/// - `!`: never returns.
fn usage_and_exit() -> ! {
    eprintln!("Usage:");
    eprintln!("  acoustic_ofdm_cli encode [--base-freq-hz HZ] <out.wav> <payload_text>");
    eprintln!("  acoustic_ofdm_cli decode [--base-freq-hz HZ] <in.wav>");
    eprintln!("  acoustic_ofdm_cli roundtrip [--base-freq-hz HZ] <wav_path> <payload_text>");
    std::process::exit(2);
}

/// Applies CLI option overrides to modem configuration.
///
/// Parameters:
/// - `cfg`: modem configuration to mutate.
/// - `args`: command arguments excluding executable and subcommand.
/// Returns:
/// - `Result<Vec<String>, Box<dyn Error>>`: positional arguments after options.
fn apply_cli_overrides(cfg: &mut OfdmConfig, args: &[String]) -> Result<Vec<String>, Box<dyn Error>> {
    let mut pos = Vec::new();
    let mut i = 0usize;
    while i < args.len() {
        match args[i].as_str() {
            "--base-freq-hz" => {
                if i + 1 >= args.len() {
                    return Err("--base-freq-hz requires a value".into());
                }
                let hz: f32 = args[i + 1].parse()?;
                if !hz.is_finite() || hz <= 0.0 {
                    return Err("base frequency must be a positive finite number".into());
                }
                cfg.base_freq_hz = Some(hz);
                i += 2;
            }
            s if s.starts_with("--") => {
                return Err(format!("unknown option: {s}").into());
            }
            _ => {
                pos.push(args[i].clone());
                i += 1;
            }
        }
    }
    Ok(pos)
}

/// Runs `encode` command.
///
/// Parameters:
/// - `out_path`: destination WAV path.
/// - `payload`: payload bytes.
/// - `cfg`: modem configuration.
/// Returns:
/// - `Result<(), Box<dyn Error>>`: `Ok(())` on success.
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

/// Runs `decode` command.
///
/// Parameters:
/// - `in_path`: source WAV path.
/// - `cfg`: modem configuration.
/// Returns:
/// - `Result<(), Box<dyn Error>>`: `Ok(())` on success.
fn cmd_decode(in_path: &Path, cfg: &OfdmConfig) -> Result<(), Box<dyn Error>> {
    let (samples, sr) = load_wav_mono_f32(in_path)?;
    if sr != cfg.fs as u32 {
        return Err(format!("sample-rate mismatch: wav={} cfg={}", sr, cfg.fs as u32).into());
    }
    let payload = decode_single_packet_passband(&samples, cfg).ok_or("decode failed")?;
    let hex = payload
        .iter()
        .map(|b| format!("{b:02X}"))
        .collect::<Vec<_>>()
        .join(" ");
    println!("Decoded {} bytes", payload.len());
    println!("HEX: {hex}");
    println!("UTF8(lossy): {}", String::from_utf8_lossy(&payload));
    Ok(())
}

/// Program entry point.
///
/// Parameters:
/// - none.
/// Returns:
/// - `Result<(), Box<dyn Error>>`: process result.
fn main() -> Result<(), Box<dyn Error>> {
    let mut cfg = OfdmConfig::default();
    let mut args = std::env::args();
    let _exe = args.next();
    let Some(cmd) = args.next() else {
        usage_and_exit();
    };
    let rest: Vec<String> = args.collect();
    let pos = apply_cli_overrides(&mut cfg, &rest)?;

    match cmd.as_str() {
        "encode" => {
            let Some(out) = pos.first() else {
                usage_and_exit();
            };
            let Some(payload_str) = pos.get(1) else {
                usage_and_exit();
            };
            if pos.len() != 2 {
                usage_and_exit();
            }
            cmd_encode(Path::new(&out), &payload_from_arg(&payload_str), &cfg)?;
        }
        "decode" => {
            let Some(inp) = pos.first() else {
                usage_and_exit();
            };
            if pos.len() != 1 {
                usage_and_exit();
            }
            cmd_decode(Path::new(&inp), &cfg)?;
        }
        "roundtrip" => {
            let Some(path) = pos.first() else {
                usage_and_exit();
            };
            let Some(payload_str) = pos.get(1) else {
                usage_and_exit();
            };
            if pos.len() != 2 {
                usage_and_exit();
            }
            let payload = payload_from_arg(&payload_str);
            cmd_encode(Path::new(&path), &payload, &cfg)?;
            cmd_decode(Path::new(&path), &cfg)?;
        }
        _ => usage_and_exit(),
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::{apply_cli_overrides, payload_from_arg};
    use acoustic_ofdm::OfdmConfig;

    /// Ensures payload parser keeps UTF-8 bytes exactly.
    ///
    /// Parameters:
    /// - none.
    /// Returns:
    /// - none.
    #[test]
    fn payload_from_arg_utf8_bytes() {
        let p = payload_from_arg("ciao-OFDM");
        assert_eq!(p, b"ciao-OFDM");
    }

    /// Ensures CLI base frequency option is parsed and applied.
    ///
    /// Parameters:
    /// - none.
    /// Returns:
    /// - none.
    #[test]
    fn parse_base_freq_override() {
        let mut cfg = OfdmConfig::default();
        let args = vec![
            "--base-freq-hz".to_string(),
            "2500".to_string(),
            "in.wav".to_string(),
        ];
        let pos = apply_cli_overrides(&mut cfg, &args).expect("parse failed");
        assert_eq!(pos, vec!["in.wav".to_string()]);
        assert_eq!(cfg.base_freq_hz, Some(2500.0));
    }
}

// vim: set ts=4 sw=4 et:
