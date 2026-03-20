// Copyright (c) 2026 Elias S. G. Carotti

use std::error::Error;
use std::io::{Read, Write};
use std::path::Path;
use std::sync::{Arc, Mutex};
use std::time::Duration;

use acoustic_ofdm::{
    decode_single_packet_passband,
    encode_single_packet_passband,
    load_wav_mono_f32,
    save_wav_mono_i16,
    OfdmConfig,
};
use cpal::traits::{DeviceTrait, HostTrait, StreamTrait};
use cpal::{SampleFormat, Stream, StreamConfig};
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use ringbuf::{traits::*, HeapRb};

/// Parses payload text into bytes.
///
/// Parameters:
/// - `s`: payload string.
/// Returns:
/// - `Vec<u8>`: UTF-8 bytes.
fn payload_from_arg(s: &str) -> Vec<u8> {
    s.as_bytes().to_vec()
}

/// Parses payload bytes from argument or stdin.
///
/// Parameters:
/// - `s`: payload argument, or `-` to read from stdin bytes.
/// Returns:
/// - `Result<Vec<u8>, Box<dyn Error>>`: payload bytes.
fn payload_from_arg_or_stdin(s: &str) -> Result<Vec<u8>, Box<dyn Error>> {
    if s == "-" {
        let mut buf = Vec::new();
        std::io::stdin().lock().read_to_end(&mut buf)?;
        Ok(buf)
    } else {
        Ok(payload_from_arg(s))
    }
}

/// Prints command usage and exits.
///
/// Parameters:
/// - `code`: process exit code.
/// Returns:
/// - `!`: never returns.
fn usage_and_exit(code: i32) -> ! {
    eprintln!("Usage:");
    eprintln!("  acoustic_ofdm_cli encode [--base-freq-hz HZ] <out.wav> <payload_text|stdin>");
    eprintln!("  acoustic_ofdm_cli decode [--base-freq-hz HZ] [--stdout] <in.wav>");
    eprintln!("  acoustic_ofdm_cli roundtrip [--base-freq-hz HZ] [--stdout] <wav_path> <payload_text|stdin>");
    eprintln!(
        "  acoustic_ofdm_cli codec-loop [--base-freq-hz HZ] [--iterations N] [--snr-db DB]"
    );
    eprintln!("      [--echo MS:GAIN] [--rand-echo-count N] [--rand-echo-max-ms MS]");
    eprintln!("      [--rand-echo-gain-min G] [--rand-echo-gain-max G] [--seed N]");
    eprintln!("      [--no-channel] <payload_text|stdin>");
    eprintln!("  acoustic_ofdm_cli tx [--base-freq-hz HZ] [--spk-gain GAIN] <payload_text|stdin>");
    eprintln!("  acoustic_ofdm_cli rx [--base-freq-hz HZ] [--duration-sec SEC] [--mic-gain GAIN] [--stdout] [--verbose]");
    eprintln!();
    eprintln!("Common options:");
    eprintln!("  --base-freq-hz HZ   OFDM base subcarrier frequency in Hz (encode/decode/roundtrip).");
    eprintln!("  --stdout            Write decoded payload as raw bytes to stdout (decode/roundtrip).");
    eprintln!();
    eprintln!("Codec-loop channel options:");
    eprintln!("  --snr-db DB         AWGN SNR in dB (default: 30).");
    eprintln!("  --echo MS:GAIN      Add one echo tap (repeatable), delay in ms and linear gain.");
    eprintln!("  --rand-echo-count N Random echoes per iteration (default: 0).");
    eprintln!("  --rand-echo-max-ms  Max random echo delay in ms (default: 3.0).");
    eprintln!("  --rand-echo-gain-min G  Min random echo gain (default: 0.05).");
    eprintln!("  --rand-echo-gain-max G  Max random echo gain (default: 0.35).");
    eprintln!("  --seed N            RNG seed for reproducible channel realizations.");
    eprintln!("  --no-channel        Disable channel model (deterministic encode/decode only).");
    eprintln!();
    eprintln!("RX/TX options:");
    eprintln!("  --duration-sec SEC  Stream duration in seconds (default: 10).");
    eprintln!("  --mic-gain GAIN     Microphone gain multiplier (default: 1.0).");
    eprintln!("  --spk-gain GAIN     Speaker gain multiplier (default: 1.0).");
    eprintln!("  --verbose           Print extra diagnostics (especially for rx).");
    eprintln!();
    eprintln!("Examples:");
    eprintln!("  acoustic_ofdm_cli encode /tmp/pkt.wav \"hello-ofdm\"");
    eprintln!("  printf \"hello\" | acoustic_ofdm_cli encode /tmp/pkt.wav -");
    eprintln!("  acoustic_ofdm_cli decode --stdout /tmp/pkt.wav > payload.bin");
    eprintln!("  acoustic_ofdm_cli roundtrip --base-freq-hz 2000 /tmp/pkt.wav \"hello\"");
    eprintln!("  acoustic_ofdm_cli codec-loop --iterations 50 --snr-db 16 \"hello\"");
    eprintln!("  acoustic_ofdm_cli codec-loop --echo 1.0:0.4 --echo 2.2:0.2 \"hello\"");
    eprintln!("  acoustic_ofdm_cli codec-loop --rand-echo-count 3 --rand-echo-max-ms 2.5 \"hello\"");
    eprintln!("  acoustic_ofdm_cli tx --spk-gain 0.8 \"hello\"");
    eprintln!("  acoustic_ofdm_cli rx --duration-sec 6 --verbose");
    std::process::exit(code);
}

/// Live audio options.
///
/// Parameters:
/// - none.
/// Returns:
/// - `AudioOpts`: parsed options with defaults.
#[derive(Clone, Debug)]
struct AudioOpts {
    duration_sec: f32,
    mic_gain: f32,
    spk_gain: f32,
    verbose: bool,
}

/// Channel options used by `codec-loop`.
///
/// Parameters:
/// - none.
/// Returns:
/// - `ChannelOpts`: parsed or default channel settings.
#[derive(Clone, Debug)]
struct ChannelOpts {
    snr_db: Option<f32>,
    echoes: Vec<(usize, f32)>,
    rand_echo_count: usize,
    rand_echo_max_ms: f32,
    rand_echo_gain_min: f32,
    rand_echo_gain_max: f32,
    seed: Option<u64>,
}

impl Default for ChannelOpts {
    fn default() -> Self {
        Self {
            snr_db: Some(30.0),
            echoes: Vec::new(),
            rand_echo_count: 0,
            rand_echo_max_ms: 3.0,
            rand_echo_gain_min: 0.05,
            rand_echo_gain_max: 0.35,
            seed: None,
        }
    }
}

/// Applies CLI option overrides to modem configuration.
///
/// Parameters:
/// - `cfg`: modem configuration to mutate.
/// - `args`: command arguments excluding executable and subcommand.
/// Returns:
/// - `Result<(Vec<String>, bool), Box<dyn Error>>`: positional args and stdout flag.
fn apply_cli_overrides(
    cfg: &mut OfdmConfig,
    args: &[String],
) -> Result<(Vec<String>, bool), Box<dyn Error>> {
    let mut pos = Vec::new();
    let mut stdout_raw = false;
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
            "--stdout" => {
                stdout_raw = true;
                i += 1;
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
    Ok((pos, stdout_raw))
}

/// Parses options for the `tx` command.
///
/// Parameters:
/// - `cfg`: modem configuration (mutable for common overrides).
/// - `args`: command args for `tx`.
/// Returns:
/// - `Result<(AudioOpts, Vec<u8>), Box<dyn Error>>`: parsed options and
///   payload.
fn parse_tx_args(
    cfg: &mut OfdmConfig,
    args: &[String],
) -> Result<(AudioOpts, Vec<u8>), Box<dyn Error>> {
    let mut opts = AudioOpts {
        duration_sec: 10.0,
        mic_gain: 1.0,
        spk_gain: 1.0,
        verbose: false,
    };
    let mut payload_arg: Option<String> = None;
    let mut i = 0usize;
    while i < args.len() {
        match args[i].as_str() {
            "--help" | "-h" => usage_and_exit(0),
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
            "--mic-gain" => {
                if i + 1 >= args.len() {
                    return Err("--mic-gain requires a value".into());
                }
                opts.mic_gain = args[i + 1].parse()?;
                i += 2;
            }
            "--spk-gain" => {
                if i + 1 >= args.len() {
                    return Err("--spk-gain requires a value".into());
                }
                opts.spk_gain = args[i + 1].parse()?;
                i += 2;
            }
            "--verbose" => {
                opts.verbose = true;
                i += 1;
            }
            s if s.starts_with("--") => return Err(format!("unknown option: {s}").into()),
            x => {
                if payload_arg.is_some() {
                    return Err("tx expects exactly one payload argument".into());
                }
                payload_arg = Some(x.to_string());
                i += 1;
            }
        }
    }
    let payload_s = payload_arg.ok_or("tx requires <payload_text|stdin>")?;
    let payload = payload_from_arg_or_stdin(&payload_s)?;
    if payload.is_empty() {
        return Err("payload must not be empty".into());
    }
    if !opts.spk_gain.is_finite() || opts.spk_gain < 0.0 {
        return Err("spk-gain must be >= 0".into());
    }
    Ok((opts, payload))
}

/// Parses options for the `rx` command.
///
/// Parameters:
/// - `cfg`: modem configuration (mutable for common overrides).
/// - `args`: command args for `rx`.
/// Returns:
/// - `Result<(AudioOpts, bool), Box<dyn Error>>`: parsed options and stdout flag.
fn parse_rx_args(cfg: &mut OfdmConfig, args: &[String]) -> Result<(AudioOpts, bool), Box<dyn Error>> {
    let mut opts = AudioOpts {
        duration_sec: 10.0,
        mic_gain: 1.0,
        spk_gain: 1.0,
        verbose: false,
    };
    let mut stdout_raw = false;
    let mut i = 0usize;
    while i < args.len() {
        match args[i].as_str() {
            "--help" | "-h" => usage_and_exit(0),
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
            "--stdout" => {
                stdout_raw = true;
                i += 1;
            }
            "--duration-sec" => {
                if i + 1 >= args.len() {
                    return Err("--duration-sec requires a value".into());
                }
                opts.duration_sec = args[i + 1].parse()?;
                i += 2;
            }
            "--mic-gain" => {
                if i + 1 >= args.len() {
                    return Err("--mic-gain requires a value".into());
                }
                opts.mic_gain = args[i + 1].parse()?;
                i += 2;
            }
            "--verbose" => {
                opts.verbose = true;
                i += 1;
            }
            s if s.starts_with("--") => return Err(format!("unknown option: {s}").into()),
            _ => return Err("rx does not accept positional arguments".into()),
        }
    }
    if !opts.duration_sec.is_finite() || opts.duration_sec <= 0.0 {
        return Err("duration must be > 0".into());
    }
    if !opts.mic_gain.is_finite() || opts.mic_gain < 0.0 {
        return Err("mic-gain must be >= 0".into());
    }
    Ok((opts, stdout_raw))
}

/// Parses options for the `codec-loop` command.
///
/// Parameters:
/// - `cfg`: modem configuration (mutable for common overrides).
/// - `args`: command args for the `codec-loop` subcommand.
/// Returns:
/// - `Result<(Vec<u8>, usize, ChannelOpts), Box<dyn Error>>`: payload bytes,
///   iterations and channel options.
fn parse_codec_loop_args(
    cfg: &mut OfdmConfig,
    args: &[String],
) -> Result<(Vec<u8>, usize, ChannelOpts), Box<dyn Error>> {
    let mut iterations: usize = 20;
    let mut payload_arg: Option<String> = None;
    let mut ch = ChannelOpts::default();
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
            "--iterations" => {
                if i + 1 >= args.len() {
                    return Err("--iterations requires a value".into());
                }
                iterations = args[i + 1].parse()?;
                i += 2;
            }
            "--snr-db" => {
                if i + 1 >= args.len() {
                    return Err("--snr-db requires a value".into());
                }
                let snr_db: f32 = args[i + 1].parse()?;
                if !snr_db.is_finite() {
                    return Err("snr-db must be finite".into());
                }
                ch.snr_db = Some(snr_db);
                i += 2;
            }
            "--seed" => {
                if i + 1 >= args.len() {
                    return Err("--seed requires a value".into());
                }
                ch.seed = Some(args[i + 1].parse()?);
                i += 2;
            }
            "--echo" => {
                if i + 1 >= args.len() {
                    return Err("--echo requires MS:GAIN".into());
                }
                let parts: Vec<&str> = args[i + 1].split(':').collect();
                if parts.len() != 2 {
                    return Err("--echo must be in MS:GAIN format".into());
                }
                let delay_ms: f32 = parts[0].parse()?;
                let gain: f32 = parts[1].parse()?;
                if !delay_ms.is_finite() || delay_ms < 0.0 {
                    return Err("echo delay must be finite and >= 0 ms".into());
                }
                if !gain.is_finite() {
                    return Err("echo gain must be finite".into());
                }
                let delay_samples = ((delay_ms * cfg.fs) / 1000.0).round() as usize;
                ch.echoes.push((delay_samples, gain));
                i += 2;
            }
            "--rand-echo-count" => {
                if i + 1 >= args.len() {
                    return Err("--rand-echo-count requires a value".into());
                }
                ch.rand_echo_count = args[i + 1].parse()?;
                i += 2;
            }
            "--rand-echo-max-ms" => {
                if i + 1 >= args.len() {
                    return Err("--rand-echo-max-ms requires a value".into());
                }
                ch.rand_echo_max_ms = args[i + 1].parse()?;
                i += 2;
            }
            "--rand-echo-gain-min" => {
                if i + 1 >= args.len() {
                    return Err("--rand-echo-gain-min requires a value".into());
                }
                ch.rand_echo_gain_min = args[i + 1].parse()?;
                i += 2;
            }
            "--rand-echo-gain-max" => {
                if i + 1 >= args.len() {
                    return Err("--rand-echo-gain-max requires a value".into());
                }
                ch.rand_echo_gain_max = args[i + 1].parse()?;
                i += 2;
            }
            "--no-channel" => {
                ch.snr_db = None;
                ch.echoes.clear();
                ch.rand_echo_count = 0;
                i += 1;
            }
            s if s.starts_with("--") => return Err(format!("unknown option: {s}").into()),
            x => {
                if payload_arg.is_some() {
                    return Err("codec-loop expects exactly one payload argument".into());
                }
                payload_arg = Some(x.to_string());
                i += 1;
            }
        }
    }
    let payload_s = payload_arg.ok_or("codec-loop requires <payload_text|stdin>")?;
    let payload = payload_from_arg_or_stdin(&payload_s)?;
    if payload.is_empty() {
        return Err("payload must not be empty".into());
    }
    if iterations == 0 {
        return Err("iterations must be > 0".into());
    }
    if !ch.rand_echo_max_ms.is_finite() || ch.rand_echo_max_ms < 0.0 {
        return Err("rand-echo-max-ms must be finite and >= 0".into());
    }
    if !ch.rand_echo_gain_min.is_finite()
        || !ch.rand_echo_gain_max.is_finite()
        || ch.rand_echo_gain_min > ch.rand_echo_gain_max
    {
        return Err("random echo gain bounds are invalid".into());
    }
    Ok((payload, iterations, ch))
}

/// Converts i16 sample to normalized f32 [-1, 1].
///
/// Parameters:
/// - `x`: signed integer sample.
/// Returns:
/// - `f32`: normalized sample.
fn i16_to_f32(x: i16) -> f32 {
    (x as f32) / (i16::MAX as f32)
}

/// Converts u16 sample to normalized f32 [-1, 1].
///
/// Parameters:
/// - `x`: unsigned integer sample.
/// Returns:
/// - `f32`: normalized sample.
fn u16_to_f32(x: u16) -> f32 {
    (x as f32) / (u16::MAX as f32) * 2.0 - 1.0
}

/// Converts normalized f32 [-1, 1] to i16.
///
/// Parameters:
/// - `x`: normalized sample.
/// Returns:
/// - `i16`: quantized sample.
fn f32_to_i16(x: f32) -> i16 {
    let y = x.clamp(-1.0, 1.0);
    (y * (i16::MAX as f32)) as i16
}

/// Converts normalized f32 [-1, 1] to u16.
///
/// Parameters:
/// - `x`: normalized sample.
/// Returns:
/// - `u16`: quantized sample.
fn f32_to_u16(x: f32) -> u16 {
    let y = x.clamp(-1.0, 1.0);
    (((y + 1.0) * 0.5) * (u16::MAX as f32)) as u16
}

/// Builds input stream and pushes mono frames into ring buffer.
///
/// Parameters:
/// - `dev`: input device.
/// - `cfg`: stream configuration.
/// - `fmt`: input sample format.
/// - `gain`: output gain.
/// - `prod`: ring-buffer producer.
/// Returns:
/// - `Result<Stream, Box<dyn Error>>`: input stream.
fn build_input_stream(
    dev: &cpal::Device,
    cfg: &StreamConfig,
    fmt: SampleFormat,
    gain: f32,
    prod: Arc<Mutex<ringbuf::HeapProd<f32>>>,
) -> Result<Stream, Box<dyn Error>> {
    let ch = cfg.channels as usize;
    let err_fn = |e| eprintln!("input stream error: {e}");
    let s = match fmt {
        SampleFormat::F32 => dev.build_input_stream(
            cfg,
            move |data: &[f32], _| {
                if let Ok(mut p) = prod.lock() {
                    for fr in data.chunks(ch) {
                        let _ = p.try_push(fr.first().copied().unwrap_or(0.0).mul_add(gain, 0.0).clamp(-1.0, 1.0));
                    }
                }
            },
            err_fn,
            None,
        )?,
        SampleFormat::I16 => dev.build_input_stream(
            cfg,
            move |data: &[i16], _| {
                if let Ok(mut p) = prod.lock() {
                    for fr in data.chunks(ch) {
                        let s = i16_to_f32(*fr.first().unwrap_or(&0)) * gain;
                        let _ = p.try_push(s.clamp(-1.0, 1.0));
                    }
                }
            },
            err_fn,
            None,
        )?,
        SampleFormat::U16 => dev.build_input_stream(
            cfg,
            move |data: &[u16], _| {
                if let Ok(mut p) = prod.lock() {
                    for fr in data.chunks(ch) {
                        let s = u16_to_f32(*fr.first().unwrap_or(&0)) * gain;
                        let _ = p.try_push(s.clamp(-1.0, 1.0));
                    }
                }
            },
            err_fn,
            None,
        )?,
        _ => return Err(format!("unsupported input sample format: {fmt:?}").into()),
    };
    Ok(s)
}

/// Builds output stream and pulls mono frames from ring buffer.
///
/// Parameters:
/// - `dev`: output device.
/// - `cfg`: stream configuration.
/// - `fmt`: output sample format.
/// - `cons`: ring-buffer consumer.
/// Returns:
/// - `Result<Stream, Box<dyn Error>>`: output stream.
fn build_output_stream(
    dev: &cpal::Device,
    cfg: &StreamConfig,
    fmt: SampleFormat,
    cons: Arc<Mutex<ringbuf::HeapCons<f32>>>,
    gain: f32,
) -> Result<Stream, Box<dyn Error>> {
    let ch = cfg.channels as usize;
    let err_fn = |e| eprintln!("output stream error: {e}");
    let s = match fmt {
        SampleFormat::F32 => dev.build_output_stream(
            cfg,
            move |data: &mut [f32], _| {
                if let Ok(mut c) = cons.lock() {
                    for fr in data.chunks_mut(ch) {
                        let v = (c.try_pop().unwrap_or(0.0) * gain).clamp(-1.0, 1.0);
                        for y in fr {
                            *y = v;
                        }
                    }
                }
            },
            err_fn,
            None,
        )?,
        SampleFormat::I16 => dev.build_output_stream(
            cfg,
            move |data: &mut [i16], _| {
                if let Ok(mut c) = cons.lock() {
                    for fr in data.chunks_mut(ch) {
                        let v = f32_to_i16((c.try_pop().unwrap_or(0.0) * gain).clamp(-1.0, 1.0));
                        for y in fr {
                            *y = v;
                        }
                    }
                }
            },
            err_fn,
            None,
        )?,
        SampleFormat::U16 => dev.build_output_stream(
            cfg,
            move |data: &mut [u16], _| {
                if let Ok(mut c) = cons.lock() {
                    for fr in data.chunks_mut(ch) {
                        let v = f32_to_u16((c.try_pop().unwrap_or(0.0) * gain).clamp(-1.0, 1.0));
                        for y in fr {
                            *y = v;
                        }
                    }
                }
            },
            err_fn,
            None,
        )?,
        _ => return Err(format!("unsupported output sample format: {fmt:?}").into()),
    };
    Ok(s)
}

/// Generates one zero-mean unit-variance Gaussian random sample.
///
/// Parameters:
/// - `rng`: random generator.
/// Returns:
/// - `f32`: Gaussian random variable.
fn gaussian_sample(rng: &mut StdRng) -> f32 {
    let u1 = rng.random::<f32>().clamp(1e-7, 1.0);
    let u2 = rng.random::<f32>();
    (-2.0 * u1.ln()).sqrt() * (2.0 * std::f32::consts::PI * u2).cos()
}

/// Applies simple multipath and AWGN channel to a waveform.
///
/// Parameters:
/// - `tx`: transmit samples.
/// - `ch`: channel options.
/// - `rng`: random generator.
/// Returns:
/// - `Vec<f32>`: channel output waveform.
fn apply_channel(tx: &[f32], ch: &ChannelOpts, rng: &mut StdRng) -> Vec<f32> {
    let mut y = tx.to_vec();
    for &(delay_samp, gain) in &ch.echoes {
        if delay_samp == 0 {
            for (yi, &x) in y.iter_mut().zip(tx.iter()) {
                *yi += gain * x;
            }
            continue;
        }
        if y.len() < tx.len() + delay_samp {
            y.resize(tx.len() + delay_samp, 0.0);
        }
        for (i, &x) in tx.iter().enumerate() {
            y[i + delay_samp] += gain * x;
        }
    }
    if let Some(snr_db) = ch.snr_db {
        let sig_pow = y.iter().map(|v| v * v).sum::<f32>() / (y.len().max(1) as f32);
        let snr_lin = 10.0f32.powf(snr_db / 10.0);
        let noise_var = if snr_lin > 0.0 { sig_pow / snr_lin } else { sig_pow };
        let sigma = noise_var.sqrt();
        for v in &mut y {
            *v += sigma * gaussian_sample(rng);
        }
    }
    y
}

/// Runs encode/decode self-test loop and prints summary stats.
///
/// Parameters:
/// - `payload`: payload bytes to test.
/// - `iterations`: number of encode/decode runs.
/// - `cfg`: modem configuration.
/// - `ch`: channel options.
/// Returns:
/// - `Result<(), Box<dyn Error>>`: `Ok(())` on success.
fn cmd_codec_loop(
    payload: &[u8],
    iterations: usize,
    cfg: &OfdmConfig,
    ch: &ChannelOpts,
) -> Result<(), Box<dyn Error>> {
    if payload.len() > cfg.packet_payload_bytes {
        return Err(format!(
            "payload too long for single-packet app: {} > {}",
            payload.len(),
            cfg.packet_payload_bytes
        )
        .into());
    }
    let seed = ch.seed.unwrap_or_else(|| rand::rng().random::<u64>());
    let mut rng = StdRng::seed_from_u64(seed);
    let mut ok = 0usize;
    println!("codec-loop channel:");
    println!("  seed       : {}", seed);
    println!(
        "  snr_db     : {}",
        ch.snr_db
            .map(|v| format!("{v:.2}"))
            .unwrap_or_else(|| "none".to_string())
    );
    println!("  fixed_echoes: {}", ch.echoes.len());
    for (idx, &(d, g)) in ch.echoes.iter().enumerate() {
        println!("    [{}] delay={} samples gain={:.3}", idx, d, g);
    }
    println!(
        "  random_echoes: {} per iter (max_delay_ms={:.2}, gain=[{:.3},{:.3}])",
        ch.rand_echo_count,
        ch.rand_echo_max_ms,
        ch.rand_echo_gain_min,
        ch.rand_echo_gain_max
    );
    println!();
    for i in 0..iterations {
        let mut iter_ch = ch.clone();
        if ch.rand_echo_count > 0 {
            let max_delay_samples = ((ch.rand_echo_max_ms * cfg.fs) / 1000.0).round() as usize;
            for _ in 0..ch.rand_echo_count {
                let d = if max_delay_samples == 0 {
                    0
                } else {
                    rng.random_range(0..=max_delay_samples)
                };
                let g = if (ch.rand_echo_gain_max - ch.rand_echo_gain_min).abs() < f32::EPSILON {
                    ch.rand_echo_gain_min
                } else {
                    rng.random_range(ch.rand_echo_gain_min..=ch.rand_echo_gain_max)
                };
                iter_ch.echoes.push((d, g));
            }
        }
        let tx = encode_single_packet_passband(payload, cfg);
        let rx = apply_channel(&tx, &iter_ch, &mut rng);
        let dec = decode_single_packet_passband(&rx, cfg);
        let pass = dec.as_deref() == Some(payload);
        if pass {
            ok += 1;
        }
        println!(
            "[{}/{}] {}",
            i + 1,
            iterations,
            if pass { "OK" } else { "FAIL" }
        );
    }
    println!();
    println!("codec-loop summary:");
    println!("  iterations : {}", iterations);
    println!("  ok         : {}", ok);
    println!("  fail       : {}", iterations - ok);
    println!("  success    : {:.1}%", 100.0 * (ok as f32) / (iterations as f32));
    Ok(())
}

/// Runs `tx`: OFDM-encodes payload and transmits through speaker.
///
/// Parameters:
/// - `payload`: payload bytes to transmit.
/// - `cfg`: modem configuration.
/// - `opts`: TX options.
/// Returns:
/// - `Result<(), Box<dyn Error>>`: `Ok(())` on success.
fn cmd_tx(payload: &[u8], cfg: &OfdmConfig, opts: &AudioOpts) -> Result<(), Box<dyn Error>> {
    if payload.len() > cfg.packet_payload_bytes {
        return Err(format!(
            "payload too long for single-packet app: {} > {}",
            payload.len(),
            cfg.packet_payload_bytes
        )
        .into());
    }

    let host = cpal::default_host();
    let out_dev = host.default_output_device().ok_or("no output device")?;
    let out_cfg = out_dev.default_output_config()?;

    let mut cfg_rt = cfg.clone();
    cfg_rt.fs = out_cfg.sample_rate().0 as f32;
    let tx = encode_single_packet_passband(payload, &cfg_rt);

    let out_cap = tx.len() + 4096;
    let out_rb = HeapRb::<f32>::new(out_cap.max(4096));
    let (mut out_prod, out_cons) = out_rb.split();
    for &s in &tx {
        let _ = out_prod.try_push(s);
    }
    let out_cons = Arc::new(Mutex::new(out_cons));

    let out_stream = build_output_stream(
        &out_dev,
        &out_cfg.clone().into(),
        out_cfg.sample_format(),
        out_cons,
        opts.spk_gain,
    )?;

    println!("Output device: {}", out_dev.name()?);
    println!("Stream config: {} Hz, out {:?}", out_cfg.sample_rate().0, out_cfg.sample_format());
    println!("Transmit samples: {}", tx.len());
    if opts.verbose {
        let tx_dur = (tx.len() as f32) / cfg_rt.fs;
        let peak = tx
            .iter()
            .fold(0.0f32, |m, &v| if v.abs() > m { v.abs() } else { m });
        println!("TX diagnostics: duration={:.3}s peak={:.3} spk_gain={:.3}", tx_dur, peak, opts.spk_gain);
    }
    out_stream.play()?;
    let play_sec = (tx.len() as f32 / cfg_rt.fs) + 0.25;
    std::thread::sleep(Duration::from_secs_f32(play_sec.max(0.25)));
    drop(out_stream);
    println!("Transmit done.");
    Ok(())
}

/// Computes basic signal diagnostics.
///
/// Parameters:
/// - `x`: input waveform.
/// Returns:
/// - `(f32, f32, usize)`: RMS, peak absolute value, and index of first sample
///   above 2% full-scale (or `x.len()` when none).
fn signal_diag(x: &[f32]) -> (f32, f32, usize) {
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

/// Runs `rx`: captures from microphone and attempts OFDM sync+decode.
///
/// Parameters:
/// - `cfg`: modem configuration.
/// - `opts`: RX options.
/// - `stdout_raw`: print decoded bytes to stdout when true.
/// Returns:
/// - `Result<(), Box<dyn Error>>`: `Ok(())` on success.
fn cmd_rx(cfg: &OfdmConfig, opts: &AudioOpts, stdout_raw: bool) -> Result<(), Box<dyn Error>> {
    let host = cpal::default_host();
    let in_dev = host.default_input_device().ok_or("no input device")?;
    let in_cfg = in_dev.default_input_config()?;
    let mut cfg_rt = cfg.clone();
    cfg_rt.fs = in_cfg.sample_rate().0 as f32;

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

    println!("Input device : {}", in_dev.name()?);
    println!("Stream config: {} Hz, in {:?}", in_cfg.sample_rate().0, in_cfg.sample_format());
    println!("Listening for {:.2}s ...", opts.duration_sec);
    in_stream.play()?;
    std::thread::sleep(Duration::from_secs_f32(opts.duration_sec));
    drop(in_stream);

    let mut rx = Vec::<f32>::new();
    while let Some(s) = in_cons.try_pop() {
        rx.push(s);
    }
    println!("Captured samples: {}", rx.len());
    let (rms, peak, first_loud) = signal_diag(&rx);
    if opts.verbose {
        println!(
            "RX diagnostics: rms={:.5} peak={:.5} first_loud_sample={} ({:.3}s)",
            rms,
            peak,
            first_loud,
            (first_loud as f32) / cfg_rt.fs
        );
    }

    let mut dec = decode_single_packet_passband(&rx, &cfg_rt);
    if dec.is_none() {
        let coarse_step = ((cfg_rt.fs * 0.02).round() as usize).max(1);
        let max_off = rx
            .len()
            .saturating_sub(((cfg_rt.fs * 0.15).round() as usize).max(1));
        if opts.verbose {
            println!(
                "Decode retry sweep: step={} samples (~{:.1} ms), max_offset={}",
                coarse_step,
                1000.0 * (coarse_step as f32) / cfg_rt.fs,
                max_off
            );
        }
        let mut attempts = 0usize;
        for off in (0..=max_off).step_by(coarse_step) {
            attempts += 1;
            if let Some(bytes) = decode_single_packet_passband(&rx[off..], &cfg_rt) {
                if opts.verbose {
                    println!("Decode recovered at offset {} samples ({:.3}s)", off, (off as f32) / cfg_rt.fs);
                }
                dec = Some(bytes);
                break;
            }
        }
        if opts.verbose {
            println!("Decode attempts: {}", attempts + 1);
        }
    }

    match dec {
        Some(bytes) => {
            println!("Decode: OK ({} bytes)", bytes.len());
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
                println!("HEX: {hex}");
                println!("UTF8(lossy): {}", String::from_utf8_lossy(&bytes));
            }
        }
        None => {
            println!("Decode: FAIL (sync/CFO/equalization/CRC path)");
        }
    }
    Ok(())
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
fn cmd_decode(in_path: &Path, cfg: &OfdmConfig, stdout_raw: bool) -> Result<(), Box<dyn Error>> {
    let (samples, sr) = load_wav_mono_f32(in_path)?;
    if sr != cfg.fs as u32 {
        return Err(format!("sample-rate mismatch: wav={} cfg={}", sr, cfg.fs as u32).into());
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

/// Program entry point.
///
/// Parameters:
/// - none.
/// Returns:
/// - `Result<(), Box<dyn Error>>`: process result.
fn main() -> Result<(), Box<dyn Error>> {
    let mut args = std::env::args();
    let _exe = args.next();
    let Some(cmd) = args.next() else {
        usage_and_exit(2);
    };
    if cmd == "--help" || cmd == "-h" {
        usage_and_exit(0);
    }
    let rest: Vec<String> = args.collect();

    match cmd.as_str() {
        "encode" => {
            let mut cfg = OfdmConfig::default();
            let (pos, _stdout_raw) = apply_cli_overrides(&mut cfg, &rest)?;
            let Some(out) = pos.first() else {
                usage_and_exit(2);
            };
            let Some(payload_str) = pos.get(1) else {
                usage_and_exit(2);
            };
            if pos.len() != 2 {
                usage_and_exit(2);
            }
            let payload = payload_from_arg_or_stdin(payload_str)?;
            cmd_encode(Path::new(out), &payload, &cfg)?;
        }
        "decode" => {
            let mut cfg = OfdmConfig::default();
            let (pos, stdout_raw) = apply_cli_overrides(&mut cfg, &rest)?;
            let Some(inp) = pos.first() else {
                usage_and_exit(2);
            };
            if pos.len() != 1 {
                usage_and_exit(2);
            }
            cmd_decode(Path::new(inp), &cfg, stdout_raw)?;
        }
        "roundtrip" => {
            let mut cfg = OfdmConfig::default();
            let (pos, stdout_raw) = apply_cli_overrides(&mut cfg, &rest)?;
            let Some(path) = pos.first() else {
                usage_and_exit(2);
            };
            let Some(payload_str) = pos.get(1) else {
                usage_and_exit(2);
            };
            if pos.len() != 2 {
                usage_and_exit(2);
            }
            let payload = payload_from_arg_or_stdin(payload_str)?;
            cmd_encode(Path::new(path), &payload, &cfg)?;
            cmd_decode(Path::new(path), &cfg, stdout_raw)?;
        }
        "codec-loop" => {
            let mut cfg = OfdmConfig::default();
            let (payload, iterations, ch) = parse_codec_loop_args(&mut cfg, &rest)?;
            cmd_codec_loop(&payload, iterations, &cfg, &ch)?;
        }
        "rx" => {
            let mut cfg = OfdmConfig::default();
            let (opts, stdout_raw) = parse_rx_args(&mut cfg, &rest)?;
            cmd_rx(&cfg, &opts, stdout_raw)?;
        }
        "tx" => {
            let mut cfg = OfdmConfig::default();
            let (opts, payload) = parse_tx_args(&mut cfg, &rest)?;
            cmd_tx(&payload, &cfg, &opts)?;
        }
        _ => usage_and_exit(2),
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::{
        apply_cli_overrides, parse_codec_loop_args, parse_rx_args, parse_tx_args, payload_from_arg,
    };
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
        let (pos, stdout_raw) = apply_cli_overrides(&mut cfg, &args).expect("parse failed");
        assert_eq!(pos, vec!["in.wav".to_string()]);
        assert_eq!(cfg.base_freq_hz, Some(2500.0));
        assert!(!stdout_raw);
    }

    /// Ensures stdout flag is parsed as a common option.
    ///
    /// Parameters:
    /// - none.
    /// Returns:
    /// - none.
    #[test]
    fn parse_stdout_override() {
        let mut cfg = OfdmConfig::default();
        let args = vec!["--stdout".to_string(), "in.wav".to_string()];
        let (pos, stdout_raw) = apply_cli_overrides(&mut cfg, &args).expect("parse failed");
        assert_eq!(pos, vec!["in.wav".to_string()]);
        assert!(stdout_raw);
    }

    /// Ensures RX option parser accepts duration and mic gain.
    ///
    /// Parameters:
    /// - none.
    /// Returns:
    /// - none.
    #[test]
    fn parse_rx_options() {
        let mut cfg = OfdmConfig::default();
        let args = vec![
            "--duration-sec".to_string(),
            "2.5".to_string(),
            "--mic-gain".to_string(),
            "0.8".to_string(),
            "--verbose".to_string(),
        ];
        let (o, stdout_raw) = parse_rx_args(&mut cfg, &args).expect("parse failed");
        assert!((o.duration_sec - 2.5).abs() < 1e-6);
        assert!((o.mic_gain - 0.8).abs() < 1e-6);
        assert!(o.verbose);
        assert!(!stdout_raw);
    }

    /// Ensures TX option parser accepts payload and speaker gain.
    ///
    /// Parameters:
    /// - none.
    /// Returns:
    /// - none.
    #[test]
    fn parse_tx_options() {
        let mut cfg = OfdmConfig::default();
        let args = vec![
            "--spk-gain".to_string(),
            "0.7".to_string(),
            "--verbose".to_string(),
            "hello".to_string(),
        ];
        let (o, payload) = parse_tx_args(&mut cfg, &args).expect("parse failed");
        assert!((o.spk_gain - 0.7).abs() < 1e-6);
        assert!(o.verbose);
        assert_eq!(payload, b"hello");
    }

    /// Ensures codec-loop parser handles iterations and payload.
    ///
    /// Parameters:
    /// - none.
    /// Returns:
    /// - none.
    #[test]
    fn parse_codec_loop_options() {
        let mut cfg = OfdmConfig::default();
        let args = vec![
            "--base-freq-hz".to_string(),
            "2200".to_string(),
            "--iterations".to_string(),
            "7".to_string(),
            "--snr-db".to_string(),
            "18".to_string(),
            "--echo".to_string(),
            "1.5:0.3".to_string(),
            "--rand-echo-count".to_string(),
            "2".to_string(),
            "--rand-echo-max-ms".to_string(),
            "2.0".to_string(),
            "--rand-echo-gain-min".to_string(),
            "0.1".to_string(),
            "--rand-echo-gain-max".to_string(),
            "0.25".to_string(),
            "--seed".to_string(),
            "123".to_string(),
            "hello".to_string(),
        ];
        let (p, n, ch) = parse_codec_loop_args(&mut cfg, &args).expect("parse failed");
        assert_eq!(n, 7);
        assert_eq!(p, b"hello");
        assert_eq!(cfg.base_freq_hz, Some(2200.0));
        assert_eq!(ch.snr_db, Some(18.0));
        assert_eq!(ch.seed, Some(123));
        assert_eq!(ch.echoes.len(), 1);
        assert!(ch.echoes[0].0 > 0);
        assert!((ch.echoes[0].1 - 0.3).abs() < 1e-6);
        assert_eq!(ch.rand_echo_count, 2);
        assert!((ch.rand_echo_max_ms - 2.0).abs() < 1e-6);
        assert!((ch.rand_echo_gain_min - 0.1).abs() < 1e-6);
        assert!((ch.rand_echo_gain_max - 0.25).abs() < 1e-6);
    }
}

// vim: set ts=4 sw=4 et:
