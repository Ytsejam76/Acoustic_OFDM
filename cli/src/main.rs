// Copyright (c) 2026 Elias S. G. Carotti

use std::error::Error;
use std::io::{Read, Write};
use std::path::Path;
use std::sync::{Arc, Mutex};
use std::time::{Duration, Instant};

mod live_profile;
mod logging;

use acoustic_ofdm::{
    diagnose_passband_window,
    dump_passband_constellation,
    dump_passband_sync_metric,
    decode_single_packet_passband,
    encode_single_packet_passband,
    load_wav_mono_f32,
    save_spectrogram_png,
    save_wav_mono_i16,
    OfdmConfig,
    WakePreamble,
};
use clap::{Args, Parser, Subcommand, ValueEnum};
use cpal::traits::{DeviceTrait, HostTrait, StreamTrait};
use cpal::{SampleFormat, Stream, StreamConfig};
use live_profile::{rx_default_log_file, tx_default_log_file, LiveProfileArg};
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use ringbuf::{traits::*, HeapRb};
use logging::init_logging;

#[derive(Copy, Clone, Debug, Eq, PartialEq, ValueEnum)]
enum WakePreambleArg {
    Gold,
    Pn,
    Chirp,
    Tone,
}

impl From<WakePreambleArg> for WakePreamble {
    fn from(value: WakePreambleArg) -> Self {
        match value {
            WakePreambleArg::Gold => WakePreamble::Gold,
            WakePreambleArg::Pn => WakePreamble::Pn,
            WakePreambleArg::Chirp => WakePreamble::Chirp,
            WakePreambleArg::Tone => WakePreamble::Tone,
        }
    }
}

#[derive(Debug, Clone, Args, Default)]
struct CommonCfgArgs {
    #[arg(long)]
    base_freq_hz: Option<f32>,
    #[arg(long, value_enum)]
    wake_preamble: Option<WakePreambleArg>,
}

#[derive(Debug, Parser)]
#[command(name = "acoustic_ofdm_cli")]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Debug, Subcommand)]
enum Commands {
    Encode(EncodeCmd),
    Decode(DecodeCmd),
    Roundtrip(RoundtripCmd),
    #[command(name = "codec-loop")]
    CodecLoop(CodecLoopCmd),
    Rx(RxCmd),
    Tx(TxCmd),
}

#[derive(Debug, Args)]
struct EncodeCmd {
    #[command(flatten)]
    common: CommonCfgArgs,
    out_wav: String,
    payload_text: String,
}

#[derive(Debug, Args)]
struct DecodeCmd {
    #[command(flatten)]
    common: CommonCfgArgs,
    #[arg(long)]
    stdout: bool,
    in_wav: String,
}

#[derive(Debug, Args)]
struct RoundtripCmd {
    #[command(flatten)]
    common: CommonCfgArgs,
    #[arg(long)]
    stdout: bool,
    wav_path: String,
    payload_text: String,
}

#[derive(Debug, Args)]
struct CodecLoopCmd {
    #[command(flatten)]
    common: CommonCfgArgs,
    #[arg(long, default_value_t = 20)]
    iterations: usize,
    #[arg(long)]
    snr_db: Option<f32>,
    #[arg(long = "echo")]
    echoes: Vec<String>,
    #[arg(long, default_value_t = 0)]
    rand_echo_count: usize,
    #[arg(long, default_value_t = 3.0)]
    rand_echo_max_ms: f32,
    #[arg(long, default_value_t = 0.05)]
    rand_echo_gain_min: f32,
    #[arg(long, default_value_t = 0.35)]
    rand_echo_gain_max: f32,
    #[arg(long)]
    seed: Option<u64>,
    #[arg(long)]
    no_channel: bool,
    payload_text: String,
}

#[derive(Debug, Args)]
struct TxCmd {
    #[command(flatten)]
    common: CommonCfgArgs,
    #[arg(long, value_enum, default_value_t = LiveProfileArg::LiveDebug)]
    profile: LiveProfileArg,
    #[arg(long)]
    spk_gain: Option<f32>,
    #[arg(long)]
    pre_delay_sec: Option<f32>,
    #[arg(long)]
    repeats: Option<usize>,
    #[arg(long)]
    gap_sec: Option<f32>,
    #[arg(long)]
    oracle: bool,
    #[arg(long)]
    verbose: bool,
    #[arg(long)]
    log_file: Option<String>,
    payload_text: Option<String>,
}

#[derive(Debug, Args)]
struct RxCmd {
    #[command(flatten)]
    common: CommonCfgArgs,
    #[arg(long, value_enum, default_value_t = LiveProfileArg::LiveDebug)]
    profile: LiveProfileArg,
    #[arg(long)]
    duration_sec: Option<f32>,
    #[arg(long)]
    mic_gain: Option<f32>,
    #[arg(long, default_value_t = 12_000.0)]
    in_hp_hz: f32,
    #[arg(long, default_value_t = 19_000.0)]
    in_lp_hz: f32,
    #[arg(long)]
    no_input_filter: bool,
    #[arg(long)]
    dump_wav: Option<String>,
    #[arg(long)]
    spectrogram: bool,
    #[arg(long, default_value = "/tmp/rx_spectrogram.png")]
    spectrogram_path: String,
    #[arg(long)]
    oracle: bool,
    #[arg(long)]
    stdout: bool,
    #[arg(long)]
    verbose: bool,
    #[arg(long)]
    log_file: Option<String>,
}

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
    pre_delay_sec: f32,
    repeats: usize,
    gap_sec: f32,
    input_filter: bool,
    input_hp_hz: f32,
    input_lp_hz: f32,
    dump_wav: Option<String>,
    spectrogram: bool,
    spectrogram_path: String,
    oracle: bool,
    verbose: bool,
}

const ORACLE_PAYLOAD: &[u8] = b"ACOUSTIC-OFDM-ORACLE";

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

fn apply_common_cfg(cfg: &mut OfdmConfig, common: &CommonCfgArgs) -> Result<(), Box<dyn Error>> {
    if let Some(hz) = common.base_freq_hz {
        if !hz.is_finite() || hz <= 0.0 {
            return Err("base frequency must be a positive finite number".into());
        }
        cfg.base_freq_hz = Some(hz);
    }
    if let Some(w) = common.wake_preamble {
        cfg.wake_preamble = w.into();
    }
    Ok(())
}

fn apply_tx_profile_cfg(cfg: &mut OfdmConfig, cmd: &TxCmd) {
    if cmd.common.wake_preamble.is_none() {
        cfg.wake_preamble = match cmd.profile {
            LiveProfileArg::Standard => cfg.wake_preamble,
            LiveProfileArg::LiveDebug => WakePreamble::Gold,
        };
    }
}

fn apply_rx_profile_cfg(cfg: &mut OfdmConfig, cmd: &RxCmd) {
    if cmd.common.wake_preamble.is_none() {
        cfg.wake_preamble = match cmd.profile {
            LiveProfileArg::Standard => cfg.wake_preamble,
            LiveProfileArg::LiveDebug => WakePreamble::Gold,
        };
    }
}


fn rx_audio_opts(cmd: &RxCmd) -> AudioOpts {
    let (duration_sec, mic_gain, dump_wav, spectrogram, spectrogram_path, oracle, verbose) =
        match cmd.profile {
            LiveProfileArg::Standard => (
                10.0,
                1.0,
                None,
                false,
                "/tmp/rx_spectrogram.png".to_string(),
                false,
                false,
            ),
            LiveProfileArg::LiveDebug => (
                5.0,
                0.2,
                Some("/tmp/rx_capture.wav".to_string()),
                true,
                "/tmp/rx_spectrogram.png".to_string(),
                true,
                true,
            ),
        };
    AudioOpts {
        duration_sec: cmd.duration_sec.unwrap_or(duration_sec),
        mic_gain: cmd.mic_gain.unwrap_or(mic_gain),
        spk_gain: 1.0,
        pre_delay_sec: 0.2,
        repeats: 3,
        gap_sec: 0.35,
        input_filter: false,
        input_hp_hz: cmd.in_hp_hz,
        input_lp_hz: cmd.in_lp_hz,
        dump_wav: cmd.dump_wav.clone().or(dump_wav),
        spectrogram: cmd.spectrogram || spectrogram,
        spectrogram_path: if cmd.spectrogram_path != "/tmp/rx_spectrogram.png" {
            cmd.spectrogram_path.clone()
        } else {
            spectrogram_path
        },
        oracle: cmd.oracle || oracle,
        verbose: cmd.verbose || verbose,
    }
}

fn tx_audio_opts(cmd: &TxCmd) -> AudioOpts {
    let (spk_gain, pre_delay_sec, repeats, gap_sec, oracle, verbose) = match cmd.profile {
        LiveProfileArg::Standard => (1.0, 0.2, 3, 0.35, false, false),
        LiveProfileArg::LiveDebug => (0.2, 0.5, 5, 0.35, true, false),
    };
    AudioOpts {
        duration_sec: 10.0,
        mic_gain: 1.0,
        spk_gain: cmd.spk_gain.unwrap_or(spk_gain),
        pre_delay_sec: cmd.pre_delay_sec.unwrap_or(pre_delay_sec),
        repeats: cmd.repeats.unwrap_or(repeats),
        gap_sec: cmd.gap_sec.unwrap_or(gap_sec),
        input_filter: false,
        input_hp_hz: 12_000.0,
        input_lp_hz: 19_000.0,
        dump_wav: None,
        spectrogram: false,
        spectrogram_path: "/tmp/rx_spectrogram.png".to_string(),
        oracle: cmd.oracle || oracle,
        verbose: cmd.verbose || verbose,
    }
}

fn parse_echo_specs(specs: &[String], fs: f32) -> Result<Vec<(usize, f32)>, Box<dyn Error>> {
    let mut echoes = Vec::new();
    for spec in specs {
        let parts: Vec<&str> = spec.split(':').collect();
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
        let delay_samples = ((delay_ms * fs) / 1000.0).round() as usize;
        echoes.push((delay_samples, gain));
    }
    Ok(echoes)
}

fn codec_loop_channel_opts(cmd: &CodecLoopCmd, fs: f32) -> Result<ChannelOpts, Box<dyn Error>> {
    let mut ch = ChannelOpts {
        snr_db: cmd.snr_db.or(Some(30.0)),
        echoes: parse_echo_specs(&cmd.echoes, fs)?,
        rand_echo_count: cmd.rand_echo_count,
        rand_echo_max_ms: cmd.rand_echo_max_ms,
        rand_echo_gain_min: cmd.rand_echo_gain_min,
        rand_echo_gain_max: cmd.rand_echo_gain_max,
        seed: cmd.seed,
    };
    if cmd.no_channel {
        ch.snr_db = None;
        ch.echoes.clear();
        ch.rand_echo_count = 0;
    }
    if cmd.iterations == 0 {
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
    Ok(ch)
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

fn sinc(x: f32) -> f32 {
    if x.abs() < 1e-6 {
        1.0
    } else {
        (std::f32::consts::PI * x).sin() / (std::f32::consts::PI * x)
    }
}

fn fir_bandpass(len: usize, f_lo_hz: f32, f_hi_hz: f32, fs: f32) -> Vec<f32> {
    let len = len.max(3) | 1;
    let m = (len - 1) as f32 * 0.5;
    let fl = (f_lo_hz / fs).clamp(0.0, 0.49);
    let fh = (f_hi_hz / fs).clamp((fl + 1.0 / fs).min(0.49), 0.49);
    let mut h = Vec::with_capacity(len);
    for n in 0..len {
        let x = (n as f32) - m;
        let ideal = 2.0 * fh * sinc(2.0 * fh * x) - 2.0 * fl * sinc(2.0 * fl * x);
        let w = 0.54 - 0.46 * (2.0 * std::f32::consts::PI * (n as f32) / ((len - 1) as f32)).cos();
        h.push(ideal * w);
    }
    let sum = h.iter().sum::<f32>().abs().max(1e-9);
    for v in &mut h {
        *v /= sum;
    }
    h
}

fn fir_filter(x: &[f32], h: &[f32]) -> Vec<f32> {
    let mut y = vec![0.0f32; x.len()];
    for n in 0..x.len() {
        let mut acc = 0.0f32;
        let kmax = (n + 1).min(h.len());
        for k in 0..kmax {
            acc += x[n - k] * h[k];
        }
        y[n] = acc;
    }
    y
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
    cfg_rt.sync_half_len = ((0.25 * cfg_rt.fs * 0.5).round() as usize).max(64);
    cfg_rt.use_pilots = Some(true);
    let tx = encode_single_packet_passband(payload, &cfg_rt);

    let pre_n = (opts.pre_delay_sec * cfg_rt.fs).round().max(0.0) as usize;
    let gap_n = (opts.gap_sec * cfg_rt.fs).round().max(0.0) as usize;
    let total_burst = opts
        .repeats
        .saturating_mul(tx.len())
        .saturating_add(opts.repeats.saturating_sub(1).saturating_mul(gap_n));
    let total_n = pre_n.saturating_add(total_burst);
    let out_cap = total_n + 4096;
    let out_rb = HeapRb::<f32>::new(out_cap.max(4096));
    let (mut out_prod, out_cons) = out_rb.split();
    for _ in 0..pre_n {
        let _ = out_prod.try_push(0.0);
    }
    for rep in 0..opts.repeats {
        for &s in &tx {
            let _ = out_prod.try_push(s);
        }
        if rep + 1 < opts.repeats {
            for _ in 0..gap_n {
                let _ = out_prod.try_push(0.0);
            }
        }
    }
    let out_cons = Arc::new(Mutex::new(out_cons));

    let out_stream = build_output_stream(
        &out_dev,
        &out_cfg.clone().into(),
        out_cfg.sample_format(),
        out_cons,
        opts.spk_gain,
    )?;

    info_line!("Output device: {}", out_dev.name()?);
    info_line!("Stream config: {} Hz, out {:?}", out_cfg.sample_rate().0, out_cfg.sample_format());
    info_line!("Wake preamble: {}", cfg_rt.wake_preamble.as_str());
    if opts.oracle {
        info_line!("Oracle mode: enabled ({} bytes)", payload.len());
    }
    info_line!("Transmit samples: {}", tx.len());
    if opts.verbose {
        let tx_dur = (tx.len() as f32) / cfg_rt.fs;
        let peak = tx.iter().fold(0.0f32, |m, &v| if v.abs() > m { v.abs() } else { m });
        info_line!(
            "TX diagnostics: duration={:.3}s peak={:.3} spk_gain={:.3} repeats={} pre_delay={:.2}s gap={:.2}s",
            tx_dur, peak, opts.spk_gain, opts.repeats, opts.pre_delay_sec, opts.gap_sec
        );
    }
    out_stream.play()?;
    if opts.verbose {
        for i in 0..opts.repeats {
            info_line!("TX burst {}/{}", i + 1, opts.repeats);
        }
    }
    let play_sec = (total_n as f32 / cfg_rt.fs) + 0.25;
    std::thread::sleep(Duration::from_secs_f32(play_sec.max(0.25)));
    drop(out_stream);
    info_line!("Transmit done.");
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

/// Estimates clipping severity in a captured waveform.
///
/// Parameters:
/// - `x`: input waveform.
/// Returns:
/// - `(usize, f32)`: number of near-full-scale samples and their fraction.
fn clipping_diag(x: &[f32]) -> (usize, f32) {
    if x.is_empty() {
        return (0, 0.0);
    }
    let clipped = x.iter().filter(|&&s| s.abs() >= 0.995).count();
    (clipped, (clipped as f32) / (x.len() as f32))
}

/// Filters input samples for sync detection only.
///
/// Parameters:
/// - `x`: input samples.
/// - `fs`: sample rate in Hz.
/// - `hp_hz`: high-pass cutoff in Hz.
/// - `lp_hz`: low-pass cutoff in Hz.
/// - `enabled`: when false, returns input copy unchanged.
/// Returns:
/// - `Vec<f32>`: filtered samples for wake/sync detection.
fn filter_for_sync_detection(
    x: &[f32],
    fs: f32,
    hp_hz: f32,
    lp_hz: f32,
    enabled: bool,
) -> Vec<f32> {
    if !enabled || x.is_empty() {
        return x.to_vec();
    }
    let hp = hp_hz.clamp(10.0, 0.45 * fs);
    let lp = lp_hz.clamp((hp + 10.0).min(0.49 * fs), 0.49 * fs);
    let h = fir_bandpass(129, hp, lp, fs);
    fir_filter(x, &h)
}

/// Generates a deterministic bipolar PN sequence.
///
/// Parameters:
/// - `n`: number of chips.
/// Returns:
/// - `Vec<f32>`: PN chips in `{-1, +1}`.
fn pn_sequence(n: usize) -> Vec<f32> {
    let mut state: u16 = 0x01FF;
    let mut out = Vec::with_capacity(n);
    for _ in 0..n {
        let bit = (state & 1) as u8;
        out.push(if bit == 0 { -1.0 } else { 1.0 });
        let fb = ((state >> 8) ^ (state >> 4)) & 1;
        state = (state >> 1) | (fb << 8);
    }
    out
}

/// Generates a deterministic bipolar Gold-like sequence.
///
/// Parameters:
/// - `n`: number of chips.
/// Returns:
/// - `Vec<f32>`: Gold chips in `{-1, +1}`.
fn gold_sequence(n: usize) -> Vec<f32> {
    let mut s1: u16 = 0x01FF;
    let mut s2: u16 = 0x0155;
    let mut out = Vec::with_capacity(n);
    for _ in 0..n {
        let b1 = (s1 & 1) as u8;
        let b2 = (s2 & 1) as u8;
        out.push(if (b1 ^ b2) == 0 { -1.0 } else { 1.0 });

        let fb1 = ((s1 >> 8) ^ (s1 >> 4)) & 1;
        let fb2 = ((s2 >> 8) ^ (s2 >> 7) ^ (s2 >> 4) ^ (s2 >> 1)) & 1;
        s1 = (s1 >> 1) | (fb1 << 8);
        s2 = (s2 >> 1) | (fb2 << 8);
    }
    out
}

/// Builds the wake preamble reference used by passband packets.
///
/// Parameters:
/// - `cfg`: modem configuration.
/// Returns:
/// - `Vec<f32>`: wake waveform samples.
fn make_wake_ref(cfg: &OfdmConfig) -> Vec<f32> {
    let n = (cfg.wake_ms * 1e-3 * cfg.fs) as usize;
    let ramp = ((0.001 * cfg.fs) as usize).min(n / 4);
    let pn = pn_sequence(n);
    let gold = gold_sequence(n);
    let mut out = vec![0.0f32; n];
    for i in 0..n {
        let t = i as f32 / cfg.fs;
        let w = match cfg.wake_preamble {
            WakePreamble::Tone => (2.0 * std::f32::consts::PI * cfg.wake_freq * t).sin(),
            WakePreamble::Chirp => {
                let tmax = ((n - 1) as f32 / cfg.fs).max(1.0 / cfg.fs);
                let k = (cfg.sync_chirp_f1 - cfg.sync_chirp_f0) / tmax;
                (2.0 * std::f32::consts::PI * (cfg.sync_chirp_f0 * t + 0.5 * k * t * t)).sin()
            }
            WakePreamble::Pn => {
                let chip = pn[i];
                chip * (2.0 * std::f32::consts::PI * cfg.wake_freq * t).sin()
            }
            WakePreamble::Gold => {
                let chip = gold[i];
                chip * (2.0 * std::f32::consts::PI * cfg.wake_freq * t).sin()
            }
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

/// Finds likely wake start indices using normalized correlation.
///
/// Parameters:
/// - `rx`: captured samples.
/// - `wake`: wake reference samples.
/// - `step`: search step in samples.
/// - `top_k`: number of candidates returned.
/// - `packet_len`: expected packet length in samples after wake start.
/// Returns:
/// - `Vec<(usize, f32)>`: `(start_index, score)` sorted by descending score.
fn wake_candidates(
    rx: &[f32],
    wake: &[f32],
    step: usize,
    top_k: usize,
    packet_len: usize,
) -> Vec<(usize, f32)> {
    if rx.len() < wake.len() || wake.is_empty() || step == 0 || top_k == 0 {
        return Vec::new();
    }
    let n = rx.len();
    let m = wake.len();
    let mut pref = vec![0.0f32; n + 1];
    for (i, &x) in rx.iter().enumerate() {
        pref[i + 1] = pref[i] + x * x;
    }
    let w_energy = wake.iter().map(|v| v * v).sum::<f32>().max(1e-12);
    let mut scored: Vec<(usize, f32)> = Vec::new();
    let min_sep = (m / 2).max(1);
    let last = n - m;
    for i in (0..=last).step_by(step) {
        let e = (pref[i + m] - pref[i]).max(1e-12);
        let mut dot = 0.0f32;
        for k in 0..m {
            dot += rx[i + k] * wake[k];
        }
        let corr = dot.abs() / (e.sqrt() * w_energy.sqrt());

        let pkt_end = i.saturating_add(packet_len).min(n);
        let post_e = (pref[pkt_end] - pref[i + m]).max(1e-12);
        let post_len = pkt_end.saturating_sub(i + m).max(1);
        let post_rms = (post_e / (post_len as f32)).sqrt();
        let wake_rms = (e / (m as f32)).sqrt();
        let energy_ratio = (post_rms / wake_rms.max(1e-6)).clamp(0.0, 4.0);

        // Favor candidates followed by sustained packet energy.
        let time_bias = 1.0 - 0.15 * ((i as f32) / (n as f32));
        let score = corr * (0.35 + 0.65 * energy_ratio) * time_bias.max(0.5);
        scored.push((i, score));
    }
    scored.sort_by(|a, b| b.1.total_cmp(&a.1));
    let mut filtered: Vec<(usize, f32)> = Vec::new();
    for (idx, sc) in scored {
        if filtered.iter().any(|(j, _)| idx.abs_diff(*j) < min_sep) {
            continue;
        }
        filtered.push((idx, sc));
        if filtered.len() == top_k {
            break;
        }
    }
    filtered
}

fn sample_linear(x: &[f32], pos: f32) -> f32 {
    if x.is_empty() || pos < 0.0 {
        return 0.0;
    }
    let i0 = pos.floor() as usize;
    if i0 >= x.len() {
        return 0.0;
    }
    let i1 = (i0 + 1).min(x.len() - 1);
    let a = pos - (i0 as f32);
    x[i0] * (1.0 - a) + x[i1] * a
}

fn fractional_wake_score(rx: &[f32], wake: &[f32], start: f32) -> f32 {
    if wake.is_empty() {
        return 0.0;
    }
    let w_energy = wake.iter().map(|v| v * v).sum::<f32>().max(1e-12);
    let mut dot = 0.0f32;
    let mut e = 0.0f32;
    for (k, &wk) in wake.iter().enumerate() {
        let s = sample_linear(rx, start + (k as f32));
        dot += s * wk;
        e += s * s;
    }
    dot.abs() / (e.sqrt().max(1e-12) * w_energy.sqrt())
}

fn refine_wake_candidates_fractional(
    rx: &[f32],
    wake: &[f32],
    cands: &[(usize, f32)],
) -> Vec<(usize, f32)> {
    let mut refined = Vec::with_capacity(cands.len());
    for (idx, base_score) in cands {
        let mut best_pos = *idx as f32;
        let mut best_score = *base_score;
        for di in -2..=2 {
            for frac_q in 0..4 {
                let pos = (*idx as f32) + (di as f32) + 0.25 * (frac_q as f32);
                if pos < 0.0 {
                    continue;
                }
                let score = fractional_wake_score(rx, wake, pos);
                if score > best_score {
                    best_score = score;
                    best_pos = pos;
                }
            }
        }
        refined.push((best_pos.round().max(0.0) as usize, best_score));
    }
    refined.sort_by(|a, b| b.1.total_cmp(&a.1));
    let mut dedup = Vec::with_capacity(refined.len());
    for (idx, sc) in refined {
        if dedup.iter().any(|(j, _)| idx.abs_diff(*j) < 4) {
            continue;
        }
        dedup.push((idx, sc));
    }
    dedup
}

#[derive(Clone, Copy, Debug)]
struct ActiveRegion {
    start: usize,
    end: usize,
    mean_rms: f32,
    peak_rms: f32,
}

fn burst_active_regions(x: &[f32], fs: f32) -> Vec<ActiveRegion> {
    if x.is_empty() {
        return Vec::new();
    }
    let win = ((0.010 * fs).round() as usize).max(1);
    let hop = ((0.002 * fs).round() as usize).max(1);
    if x.len() < win {
        let rms = (x.iter().map(|v| v * v).sum::<f32>() / (x.len() as f32)).sqrt();
        return vec![ActiveRegion {
            start: 0,
            end: x.len(),
            mean_rms: rms,
            peak_rms: rms,
        }];
    }
    let mut sum = x[..win].iter().map(|v| v * v).sum::<f32>();
    let mut env = Vec::<(usize, f32)>::new();
    let mut start = 0usize;
    loop {
        env.push((start, (sum / (win as f32)).sqrt()));
        if start + hop + win > x.len() {
            break;
        }
        for k in 0..hop {
            sum += x[start + win + k] * x[start + win + k] - x[start + k] * x[start + k];
        }
        start += hop;
    }
    let mut vals = env.iter().map(|(_, e)| *e).collect::<Vec<_>>();
    vals.sort_by(|a, b| a.total_cmp(b));
    let noise = vals[vals.len() / 5].max(1e-4);
    let th = (2.5 * noise).max(noise + 0.015);
    let min_run = ((0.025 * fs).round() as usize).max(hop);
    let pre = ((0.020 * fs).round() as usize).max(1);
    let post = ((0.180 * fs).round() as usize).max(1);
    let merge_gap = ((0.250 * fs).round() as usize).max(1);
    let mut runs = Vec::<ActiveRegion>::new();
    let mut cur: Option<(usize, usize, f32, f32, usize)> = None;
    for (s, e) in env {
        let active = e >= th;
        match (cur, active) {
            (None, true) => cur = Some((s, s + win, e, e, 1)),
            (Some((a, _b, sum_rms, peak_rms, nframes)), true) => {
                cur = Some((a, s + win, sum_rms + e, peak_rms.max(e), nframes + 1));
            }
            (Some((a, b, sum_rms, peak_rms, nframes)), false) => {
                if b.saturating_sub(a) >= min_run {
                    runs.push(ActiveRegion {
                        start: a.saturating_sub(pre),
                        end: (b + post).min(x.len()),
                        mean_rms: sum_rms / (nframes as f32),
                        peak_rms,
                    });
                }
                cur = None;
            }
            (None, false) => {}
        }
    }
    if let Some((a, b, sum_rms, peak_rms, nframes)) = cur {
        if b.saturating_sub(a) >= min_run {
            runs.push(ActiveRegion {
                start: a.saturating_sub(pre),
                end: (b + post).min(x.len()),
                mean_rms: sum_rms / (nframes as f32),
                peak_rms,
            });
        }
    }
    let mut merged = Vec::<ActiveRegion>::new();
    for r in runs {
        if let Some(last) = merged.last_mut() {
            if r.start <= last.end.saturating_add(merge_gap) {
                let last_len = last.end.saturating_sub(last.start).max(1) as f32;
                let r_len = r.end.saturating_sub(r.start).max(1) as f32;
                last.end = last.end.max(r.end);
                last.mean_rms = (last.mean_rms * last_len + r.mean_rms * r_len) / (last_len + r_len);
                last.peak_rms = last.peak_rms.max(r.peak_rms);
                continue;
            }
        }
        merged.push(r);
    }
    merged.sort_by(|a, b| {
        let sa = a.peak_rms * (0.5 + a.mean_rms) * ((a.end - a.start) as f32).sqrt();
        let sb = b.peak_rms * (0.5 + b.mean_rms) * ((b.end - b.start) as f32).sqrt();
        sb.total_cmp(&sa).then_with(|| a.start.cmp(&b.start))
    });
    merged
}

fn filter_candidates_by_regions(cands: &[(usize, f32)], regions: &[ActiveRegion]) -> Vec<(usize, f32)> {
    if regions.is_empty() {
        return Vec::new();
    }
    let keep_regions = regions.len().min(3);
    cands.iter()
        .copied()
        .filter_map(|(idx, score)| {
            regions
                .iter()
                .take(keep_regions)
                .find(|r| idx >= r.start && idx < r.end)
                .map(|r| {
                    let region_boost = (1.0 + 2.0 * r.mean_rms + 1.5 * r.peak_rms).max(1.0);
                    (idx, score * region_boost)
                })
        })
        .collect::<Vec<_>>()
        .tap_mut(|v| v.sort_by(|a, b| b.1.total_cmp(&a.1)))
}

trait TapMut: Sized {
    fn tap_mut<F: FnOnce(&mut Self)>(mut self, f: F) -> Self {
        f(&mut self);
        self
    }
}

impl<T> TapMut for T {}

/// Estimates one packet waveform length in passband samples.
///
/// Parameters:
/// - `cfg`: modem configuration.
/// Returns:
/// - `usize`: estimated packet sample length including wake and guard.
fn estimated_packet_len_samples(cfg: &OfdmConfig) -> usize {
    let bps = cfg.modulation.bits_per_symbol();
    let used_bins = cfg.used_bins.len();
    let pilots_on = cfg.use_pilots.unwrap_or(matches!(cfg.modulation, acoustic_ofdm::Modulation::Qpsk));
    let pilot_count = if pilots_on {
        cfg.num_pilots
            .unwrap_or(cfg.pilot_bins.len())
            .min(cfg.pilot_bins.len())
            .min(used_bins)
    } else {
        0
    };
    let n_data_carriers = used_bins.saturating_sub(pilot_count).max(1);
    let max_payload_bytes = cfg.packet_payload_bytes + 16;
    let max_bits = max_payload_bytes * 8;
    let bits_per_ofdm = n_data_carriers * bps;
    let n_data_ofdm = max_bits.div_ceil(bits_per_ofdm) + 2;
    let baseband_len = 2 * cfg.sync_half_len + (1 + n_data_ofdm) * (cfg.nfft + cfg.ncp);
    let wake_len = (cfg.wake_ms * 1e-3 * cfg.fs) as usize;
    let guard_len = (cfg.wake_guard_ms * 1e-3 * cfg.fs) as usize;
    wake_len + guard_len + baseband_len
}

fn ranked_offset_hypotheses(
    rx: &[f32],
    seeds: &[(usize, f32, f32, bool)],
    est_pkt: usize,
    pad: usize,
    cfg: &OfdmConfig,
) -> Vec<(usize, acoustic_ofdm::PassbandDiagnostics)> {
    let back = ((cfg.fs * 0.008).round() as isize).max(1);
    let fwd = ((cfg.fs * 0.020).round() as isize).max(1);
    let step = ((cfg.fs * 0.0005).round() as isize).max(1);
    let mut scored = Vec::new();
    for (idx, _, _, _) in seeds.iter().take(3) {
        for dj in (-back..=fwd).step_by(step as usize) {
            let off_i = *idx as isize + dj;
            if off_i < 0 {
                continue;
            }
            let off = off_i as usize;
            if off >= rx.len() {
                continue;
            }
            let end = off.saturating_add(est_pkt + pad).min(rx.len());
            if end <= off + cfg.nfft + cfg.ncp {
                continue;
            }
            let diag = diagnose_passband_window(&rx[off..end], cfg);
            let score = diagnostic_candidate_score(&diag);
            scored.push((off, score, diag));
        }
    }
    scored.sort_by(|a, b| {
        let pa = plausible_candidate(&a.2);
        let pb = plausible_candidate(&b.2);
        pb.cmp(&pa)
            .then_with(|| b.1.total_cmp(&a.1))
            .then_with(|| a.0.cmp(&b.0))
    });
    let mut dedup = Vec::new();
    let min_sep = ((cfg.fs * 0.0015).round() as usize).max(1);
    for (off, _score, diag) in scored {
        if dedup.iter().any(|(j, _)| off.abs_diff(*j) < min_sep) {
            continue;
        }
        dedup.push((off, diag));
        if dedup.len() >= 12 {
            break;
        }
    }
    dedup
}

/// Attempts a fast real-time decode on a short rolling buffer.
///
/// Parameters:
/// - `rx`: rolling capture buffer.
/// - `cfg`: modem configuration.
/// Returns:
/// - `Option<Vec<u8>>`: decoded payload when successful.
fn quick_realtime_decode(rx_raw: &[f32], rx_sync: &[f32], cfg: &OfdmConfig) -> Option<Vec<u8>> {
    let est_pkt = estimated_packet_len_samples(cfg);
    let pad = ((0.050 * cfg.fs).round() as usize).max(1);
    if rx_raw.len() < est_pkt + pad || rx_sync.len() < est_pkt + pad {
        return None;
    }
    let wake = make_wake_ref(cfg);
    let step = ((cfg.fs * 0.001).round() as usize).max(1);
    let gate = filter_for_sync_detection(rx_raw, cfg.fs, 12_000.0, 19_000.0, true);
    let regions = burst_active_regions(&gate, cfg.fs);
    let cands0 = filter_candidates_by_regions(&wake_candidates(rx_sync, &wake, step, 6, est_pkt), &regions);
    let cands = refine_wake_candidates_fractional(rx_sync, &wake, &cands0);
    for (idx, _) in cands.iter().take(4) {
        let off = *idx;
        let end = off.saturating_add(est_pkt + pad).min(rx_raw.len());
        if end > off + cfg.nfft + cfg.ncp {
            if let Some(bytes) = decode_single_packet_passband(&rx_raw[off..end], cfg) {
                return Some(bytes);
            }
        }
    }
    let ranked_cands: Vec<(usize, f32, f32, bool)> = cands
        .iter()
        .take(3)
        .filter_map(|(idx, wake_score)| {
            let end = idx.saturating_add(est_pkt + pad).min(rx_raw.len());
            (end > *idx + cfg.nfft + cfg.ncp).then(|| {
                let diag = diagnose_passband_window(&rx_raw[*idx..end], cfg);
                (*idx, *wake_score, diagnostic_candidate_score(&diag), plausible_candidate(&diag))
            })
        })
        .collect();
    for (off, _) in ranked_offset_hypotheses(rx_raw, &ranked_cands, est_pkt, pad, cfg) {
        let end = off.saturating_add(est_pkt + pad).min(rx_raw.len());
        if let Some(bytes) = decode_single_packet_passband(&rx_raw[off..end], cfg) {
            return Some(bytes);
        }
    }
    None
}

fn print_passband_diagnostics(label: &str, pkt_audio: &[f32], cfg: &OfdmConfig) {
    let d = diagnose_passband_window(pkt_audio, cfg);
    println!(
        "{}: enough={} sync_off={} cfo={:.1}Hz train_rms={:.4} hest[min/mean/max]=[{:.3}/{:.3}/{:.3}] evm[train/pilot/data]=[{:.3}/{:.3}/{:.3}] decoded={}",
        label,
        d.enough_samples,
        d.sync_off,
        d.cfo_hz,
        d.train_rms,
        d.hest_mag_min,
        d.hest_mag_mean,
        d.hest_mag_max,
        d.train_recon_evm,
        d.pilot_residual_evm,
        d.post_eq_evm,
        d.decoded
    );
}


fn diagnostic_candidate_score(diag: &acoustic_ofdm::PassbandDiagnostics) -> f32 {
    if !diag.enough_samples {
        return -1e9;
    }
    let sync_penalty = 0.0001 * (diag.sync_off as f32);
    let cfo_penalty = 0.02 * diag.cfo_hz.abs();
    let train_bonus = 120.0 * diag.train_rms;
    let hest_bonus = 5.0 * diag.hest_mag_mean - 0.6 * (diag.hest_mag_max - diag.hest_mag_min);
    let evm_penalty = 4.0 * diag.train_recon_evm + 6.0 * diag.pilot_residual_evm + 3.0 * diag.post_eq_evm;
    let decoded_bonus = if diag.decoded { 1000.0 } else { 0.0 };
    decoded_bonus + train_bonus + hest_bonus - sync_penalty - cfo_penalty - evm_penalty
}

fn plausible_candidate(diag: &acoustic_ofdm::PassbandDiagnostics) -> bool {
    diag.enough_samples
        && diag.train_rms >= 0.01
        && diag.hest_mag_mean >= 0.15
        && diag.hest_mag_max >= 0.4
        && (diag.train_recon_evm == 0.0 || diag.train_recon_evm <= 0.75)
        && (diag.pilot_residual_evm == 0.0 || diag.pilot_residual_evm <= 1.00)
        && (diag.post_eq_evm == 0.0 || diag.post_eq_evm <= 1.00)
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
    info_line!("Stream config: {} Hz, in {:?}", in_cfg.sample_rate().0, in_cfg.sample_format());
    info_line!("RX detector: wake-correlation-v2");
    info_line!("Wake preamble: {}", cfg_rt.wake_preamble.as_str());
    if opts.oracle {
        info_line!("Oracle mode: enabled (expect {} bytes)", ORACLE_PAYLOAD.len());
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
    let mut last_rt_try = Instant::now();
    while t0.elapsed() < Duration::from_secs_f32(opts.duration_sec) {
        while let Some(s) = in_cons.try_pop() {
            rx.push(s);
        }
        if last_rt_try.elapsed() >= Duration::from_millis(300) {
            let tail_span = ((cfg_rt.fs * 3.0).round() as usize).max(1);
            let st = rx.len().saturating_sub(tail_span);
            if let Some(bytes) = quick_realtime_decode(&rx[st..], &rx[st..], &cfg_rt) {
                println!(
                    "Realtime decode: OK at t={:.3}s ({} bytes)",
                    t0.elapsed().as_secs_f32(),
                    bytes.len()
                );
                if opts.oracle {
                    println!(
                        "Oracle verdict: {}",
                        if bytes.as_slice() == ORACLE_PAYLOAD { "match" } else { "mismatch" }
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
                    println!("HEX: {hex}");
                    println!("UTF8(lossy): {}", String::from_utf8_lossy(&bytes));
                }
                drop(in_stream);
                return Ok(());
            }
            if opts.verbose {
                println!(
                    "Realtime: checked t={:.3}s captured={}",
                    t0.elapsed().as_secs_f32(),
                    rx.len()
                );
            }
            last_rt_try = Instant::now();
        }
        std::thread::sleep(Duration::from_millis(10));
    }
    drop(in_stream);

    while let Some(s) = in_cons.try_pop() {
        rx.push(s);
    }
    info_line!("Captured samples: {}", rx.len());
    if let Some(path) = &opts.dump_wav {
        save_wav_mono_i16(Path::new(path), &rx, cfg_rt.fs.round() as u32)?;
        info_line!("Saved RX capture: {}", path);
        if opts.spectrogram {
            let spec_path = Path::new(&opts.spectrogram_path);
            save_spectrogram_png(spec_path, &rx, cfg_rt.fs)?;
            info_line!("Saved spectrogram PNG: {}", spec_path.display());
        }
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
        println!(
            "RX diagnostics: rms={:.5} peak={:.5} first_loud_sample={} ({:.3}s)",
            rms,
            peak,
            first_loud,
            (first_loud as f32) / cfg_rt.fs
        );
        println!(
            "RX clipping: {} samples ({:.2}%) at |x| >= 0.995",
            clipped,
            100.0 * clipped_frac
        );
        if peak < 0.01 {
            warn_line!("RX warning: very low capture level; increase speaker volume or mic gain.");
        }
        if clipped_frac >= 0.001 {
            warn_line!("RX warning: capture is clipping; reduce speaker volume, mic gain, or disable AGC.");
        } else if clipped_frac >= 0.0001 {
            warn_line!("RX warning: capture is close to clipping.");
        }
    }

    let est_pkt = estimated_packet_len_samples(&cfg_rt);
    let pad = ((0.050 * cfg_rt.fs).round() as usize).max(1);
    if opts.verbose {
        println!(
            "Decode window: est_packet={} samples ({:.3}s), pad={} samples",
            est_pkt,
            (est_pkt as f32) / cfg_rt.fs,
            pad
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
        let mut ranked_cands: Vec<(usize, f32, f32, bool)> = Vec::new();
        for (idx, wake_score) in cands.iter().copied().take(8) {
            let end = idx.saturating_add(est_pkt + pad).min(rx.len());
            if end <= idx + cfg_rt.nfft + cfg_rt.ncp {
                continue;
            }
            let diag = diagnose_passband_window(&rx[idx..end], &cfg_rt);
            let diag_score = diagnostic_candidate_score(&diag);
            let plausible = plausible_candidate(&diag);
            ranked_cands.push((idx, wake_score, diag_score, plausible));
        }
        ranked_cands.sort_by(|a, b| {
            b.3.cmp(&a.3)
                .then_with(|| b.2.total_cmp(&a.2))
                .then_with(|| b.1.total_cmp(&a.1))
                .then_with(|| a.0.cmp(&b.0))
        });
        if opts.verbose {
            println!(
                "Wake search: {} candidates (step={} samples, active_regions={})",
                cands.len(),
                coarse_step,
                active_regions.len()
            );
            for (i, r) in active_regions.iter().take(5).enumerate() {
                println!(
                    "  region {:2}: [{:.3}s, {:.3}s] mean_rms={:.4} peak_rms={:.4}",
                    i + 1,
                    (r.start as f32) / cfg_rt.fs,
                    (r.end as f32) / cfg_rt.fs,
                    r.mean_rms,
                    r.peak_rms
                );
            }
            for (i, (idx, sc)) in cands.iter().enumerate() {
                println!(
                    "  cand {:2}: idx={} t={:.3}s score={:.4}",
                    i + 1,
                    idx,
                    (*idx as f32) / cfg_rt.fs,
                    sc
                );
            }
            for (i, (idx, wake_score, diag_score, plausible)) in ranked_cands.iter().take(3).enumerate() {
                println!(
                    "  rank {:2}: idx={} t={:.3}s wake_score={:.4} diag_score={:.4} plausible={}",
                    i + 1,
                    idx,
                    (*idx as f32) / cfg_rt.fs,
                    wake_score,
                    diag_score,
                    plausible
                );
            }
            for (i, (idx, _, _, _)) in ranked_cands.iter().take(3).enumerate() {
                let off = *idx;
                let end = off.saturating_add(est_pkt + pad).min(rx.len());
                if end > off + cfg_rt.nfft + cfg_rt.ncp {
                    print_passband_diagnostics(&format!("  cand {:2} diag", i + 1), &rx[off..end], &cfg_rt);
                }
            }
            if let Some((idx, _, _, _)) = ranked_cands.first() {
                let off = *idx;
                let end = off.saturating_add(est_pkt + pad).min(rx.len());
                if end > off + cfg_rt.nfft + cfg_rt.ncp {
                    if let Some(dump) = dump_passband_constellation(&rx[off..end], &cfg_rt) {
                        let pre_path = Path::new("/tmp/ofdm_constellation_pre_eq.csv");
                        let post_path = Path::new("/tmp/ofdm_constellation_post_eq.csv");
                        {
                            let mut out = std::io::BufWriter::new(std::fs::File::create(pre_path)?);
                            writeln!(out, "re,im")?;
                            for z in &dump.pre_eq {
                                writeln!(out, "{},{}", z.re, z.im)?;
                            }
                            out.flush()?;
                        }
                        {
                            let mut out = std::io::BufWriter::new(std::fs::File::create(post_path)?);
                            writeln!(out, "re,im")?;
                            for z in &dump.post_eq {
                                writeln!(out, "{},{}", z.re, z.im)?;
                            }
                            out.flush()?;
                        }
                        println!("Saved constellation CSV: {}", pre_path.display());
                        println!("Saved constellation CSV: {}", post_path.display());
                    }
                    if let Some(sync_dump) = dump_passband_sync_metric(&rx[off..end], &cfg_rt) {
                        let sync_path = Path::new("/tmp/ofdm_sync_metric.csv");
                        let mut out = std::io::BufWriter::new(std::fs::File::create(sync_path)?);
                        writeln!(out, "offset,metric")?;
                        for (i, m) in sync_dump.metrics.iter().enumerate() {
                            writeln!(out, "{},{}", i, m)?;
                        }
                        out.flush()?;
                        println!(
                            "Saved sync metric CSV: {} (coarse_sync_off={}, refined_sync_off={})",
                            sync_path.display(),
                            sync_dump.coarse_sync_off,
                            sync_dump.refined_sync_off
                        );
                    }
                }
            }
        }
        let refined = ranked_offset_hypotheses(&rx, &ranked_cands, est_pkt, pad, &cfg_rt);
        if opts.verbose {
            println!(
                "Metric refinement: {} shortlisted offsets from top {} seeds",
                refined.len(),
                ranked_cands.len().min(3)
            );
            for (i, (off, diag)) in refined.iter().take(6).enumerate() {
                println!(
                    "  refine {:2}: off={} t={:.3}s sync_off={} train_rms={:.4} evm[train/pilot/data]=[{:.3}/{:.3}/{:.3}]",
                    i + 1,
                    off,
                    (*off as f32) / cfg_rt.fs,
                    diag.sync_off,
                    diag.train_rms,
                    diag.train_recon_evm,
                    diag.pilot_residual_evm,
                    diag.post_eq_evm
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
                    println!(
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
        println!("Decode attempts: {}", attempts);
    }

    match dec {
        Some(bytes) => {
            println!("Decode: OK ({} bytes)", bytes.len());
            if opts.oracle {
                println!(
                    "Oracle verdict: {}",
                    if bytes.as_slice() == ORACLE_PAYLOAD { "match" } else { "mismatch" }
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
        Commands::CodecLoop(cmd) => {
            let mut cfg = OfdmConfig::default();
            apply_common_cfg(&mut cfg, &cmd.common)?;
            let payload = payload_from_arg_or_stdin(&cmd.payload_text)?;
            if payload.is_empty() {
                return Err("payload must not be empty".into());
            }
            let ch = codec_loop_channel_opts(&cmd, cfg.fs)?;
            cmd_codec_loop(&payload, cmd.iterations, &cfg, &ch)?;
        }
        Commands::Rx(cmd) => {
            let mut cfg = OfdmConfig::default();
            apply_common_cfg(&mut cfg, &cmd.common)?;
            apply_rx_profile_cfg(&mut cfg, &cmd);
            let log_path = cmd
                .log_file
                .clone()
                .or_else(|| rx_default_log_file(cmd.profile));
            init_logging(log_path.as_deref(), cmd.verbose)?;
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
            init_logging(log_path.as_deref(), cmd.verbose)?;
            let opts = tx_audio_opts(&cmd);
            let payload = if cmd.oracle {
                ORACLE_PAYLOAD.to_vec()
            } else {
                let arg = cmd.payload_text.as_deref().ok_or("tx requires <payload_text|stdin>")?;
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
        Cli, Commands, LiveProfileArg, ORACLE_PAYLOAD, WakePreambleArg,
    };
    use acoustic_ofdm::{OfdmConfig, WakePreamble};
    use clap::Parser;

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

    /// Ensures common encode options are parsed by clap.
    ///
    /// Parameters:
    /// - none.
    /// Returns:
    /// - none.
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

    /// Ensures decode stdout option is parsed by clap.
    ///
    /// Parameters:
    /// - none.
    /// Returns:
    /// - none.
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

    /// Ensures RX option parser accepts duration, mic gain, and spectrogram path.
    ///
    /// Parameters:
    /// - none.
    /// Returns:
    /// - none.
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

    /// Ensures TX option parser accepts payload and speaker settings.
    ///
    /// Parameters:
    /// - none.
    /// Returns:
    /// - none.
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

    /// Ensures TX oracle mode does not require an explicit payload.
    ///
    /// Parameters:
    /// - none.
    /// Returns:
    /// - none.
    #[test]
    fn parse_tx_oracle_options() {
        let cli =
            Cli::try_parse_from(["acoustic_ofdm_cli", "tx", "--oracle", "--repeats", "2"])
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

    /// Ensures default TX live-debug profile resolves to the expected lab values.
    ///
    /// Parameters:
    /// - none.
    /// Returns:
    /// - none.
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
            }
            _ => panic!("expected tx"),
        }
    }

    /// Ensures default RX live-debug profile resolves to the expected lab values.
    ///
    /// Parameters:
    /// - none.
    /// Returns:
    /// - none.
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
            }
            _ => panic!("expected rx"),
        }
    }

    /// Ensures codec-loop parser handles channel options and payload.
    ///
    /// Parameters:
    /// - none.
    /// Returns:
    /// - none.
    #[test]
    fn parse_codec_loop_options() {
        let cli = Cli::try_parse_from([
            "acoustic_ofdm_cli",
            "codec-loop",
            "--base-freq-hz",
            "2200",
            "--iterations",
            "7",
            "--snr-db",
            "18",
            "--echo",
            "1.5:0.3",
            "--rand-echo-count",
            "2",
            "--rand-echo-max-ms",
            "2.0",
            "--rand-echo-gain-min",
            "0.1",
            "--rand-echo-gain-max",
            "0.25",
            "--seed",
            "123",
            "hello",
        ])
        .expect("parse failed");
        match cli.command {
            Commands::CodecLoop(cmd) => {
                assert_eq!(cmd.common.base_freq_hz, Some(2200.0));
                assert_eq!(cmd.iterations, 7);
                assert_eq!(cmd.snr_db, Some(18.0));
                assert_eq!(cmd.echoes, vec!["1.5:0.3"]);
                assert_eq!(cmd.rand_echo_count, 2);
                assert!((cmd.rand_echo_max_ms - 2.0).abs() < 1e-6);
                assert!((cmd.rand_echo_gain_min - 0.1).abs() < 1e-6);
                assert!((cmd.rand_echo_gain_max - 0.25).abs() < 1e-6);
                assert_eq!(cmd.seed, Some(123));
                assert_eq!(cmd.payload_text, "hello");
            }
            _ => panic!("expected codec-loop"),
        }
    }
}

// vim: set ts=4 sw=4 et:
