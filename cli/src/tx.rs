// Copyright (c) 2026 Elias S. G. Carotti

use std::error::Error;
use std::time::Duration;

use acoustic_ofdm::{encode_single_packet_passband, save_wav_mono_i16, OfdmConfig};
use cpal::traits::{DeviceTrait, HostTrait, StreamTrait};

use crate::audio::build_output_stream;
use crate::cli_args::AudioOpts;
use crate::info_line;

const CAL_TONE_SEC: f32 = 0.5;
const CAL_TONE_GAIN: f32 = 0.2;

fn calibration_frequencies(fc: f32, fs: f32) -> [f32; 3] {
    let nyquist = 0.5 * fs;
    let low = (fc - 3_000.0).clamp(1_000.0, 0.9 * nyquist);
    let mid = fc.clamp(1_000.0, 0.9 * nyquist);
    let high = (fc + 3_000.0).clamp(1_000.0, 0.9 * nyquist);
    [low, mid, high]
}

#[derive(Clone, Debug)]
pub(crate) struct TxPlan {
    pub(crate) fs: f32,
    pub(crate) packet: Vec<f32>,
    pub(crate) scheduled: Vec<f32>,
    pub(crate) burst_starts_sec: Vec<f32>,
}

pub(crate) fn build_tx_plan(
    payload: &[u8],
    cfg: &OfdmConfig,
    opts: &AudioOpts,
    fs: f32,
) -> Result<TxPlan, Box<dyn Error>> {
    if payload.len() > cfg.packet_payload_bytes {
        return Err(format!(
            "payload too long for single-packet app: {} > {}",
            payload.len(),
            cfg.packet_payload_bytes
        )
        .into());
    }

    let mut cfg_rt = cfg.clone();
    cfg_rt.fs = fs;
    cfg_rt.use_pilots = Some(true);
    let packet = encode_single_packet_passband(payload, &cfg_rt);

    let pre_n = (opts.pre_delay_sec * fs).round().max(0.0) as usize;
    let gap_n = (opts.gap_sec * fs).round().max(0.0) as usize;
    let cal_n = (CAL_TONE_SEC * fs).round().max(0.0) as usize;
    let total_burst = opts
        .repeats
        .saturating_mul(packet.len())
        .saturating_add(opts.repeats.saturating_sub(1).saturating_mul(gap_n));
    let total_n = cal_n.saturating_add(pre_n).saturating_add(total_burst);

    let mut scheduled = Vec::with_capacity(total_n);
    let mut burst_starts_sec = Vec::with_capacity(opts.repeats);
    if cal_n > 0 {
        let cal_freqs = calibration_frequencies(cfg_rt.fc, fs);
        let ramp = ((0.01 * fs).round() as usize).max(1).min(cal_n / 4);
        let seg_len = (cal_n / 3).max(1);
        for n in 0..cal_n {
            let t = n as f32 / fs;
            let env = if n < ramp {
                n as f32 / ramp as f32
            } else if n + ramp >= cal_n {
                (cal_n - 1 - n) as f32 / ramp as f32
            } else {
                1.0
            };
            let f = if n < seg_len {
                cal_freqs[0]
            } else if n < 2 * seg_len {
                cal_freqs[1]
            } else {
                cal_freqs[2]
            };
            scheduled.push(CAL_TONE_GAIN * env * (2.0 * std::f32::consts::PI * f * t).sin());
        }
    }
    scheduled.extend(std::iter::repeat_n(0.0, pre_n));
    for rep in 0..opts.repeats {
        burst_starts_sec.push((scheduled.len() as f32) / fs);
        scheduled.extend_from_slice(&packet);
        if rep + 1 < opts.repeats {
            scheduled.extend(std::iter::repeat_n(0.0, gap_n));
        }
    }

    if (opts.spk_gain - 1.0).abs() > f32::EPSILON {
        for s in &mut scheduled {
            *s *= opts.spk_gain;
        }
    }

    Ok(TxPlan {
        fs,
        packet,
        scheduled,
        burst_starts_sec,
    })
}

pub(crate) fn cmd_tx(
    payload: &[u8],
    cfg: &OfdmConfig,
    opts: &AudioOpts,
) -> Result<(), Box<dyn Error>> {
    let host = cpal::default_host();
    let out_dev = host.default_output_device().ok_or("no output device")?;
    let out_cfg = out_dev.default_output_config()?;
    let plan = build_tx_plan(payload, cfg, opts, out_cfg.sample_rate().0 as f32)?;

    if let Some(path) = &opts.dump_wav {
        save_wav_mono_i16(
            std::path::Path::new(path),
            &plan.scheduled,
            plan.fs.round() as u32,
        )?;
        info_line!("Saved TX WAV: {path}");
    }

    let out_stream = build_output_stream(
        &out_dev,
        &out_cfg.clone().into(),
        out_cfg.sample_format(),
        plan.scheduled.clone(),
    )?;

    info_line!("Output device: {}", out_dev.name()?);
    info_line!(
        "Stream config: {} Hz, {} ch, out {:?}",
        out_cfg.sample_rate().0,
        out_cfg.channels(),
        out_cfg.sample_format()
    );
    info_line!("Wake preamble: {}", cfg.wake_preamble.as_str());
    if opts.oracle {
        info_line!("Oracle mode: enabled ({} bytes)", payload.len());
    }
    info_line!("Transmit samples: {}", plan.packet.len());
    if opts.verbose {
        let cal_freqs = calibration_frequencies(cfg.fc, plan.fs);
        let peak = plan
            .packet
            .iter()
            .fold(0.0f32, |m, &v| if v.abs() > m { v.abs() } else { m });
        info_line!(
            "TX diagnostics: duration={:.3}s peak={peak:.3} cal_tones_seq={:.2}s@[{:.1},{:.1},{:.1}]Hz spk_gain={:.3} repeats={} pre_delay={:.2}s gap={:.2}s",
            (plan.packet.len() as f32) / plan.fs,
            CAL_TONE_SEC,
            cal_freqs[0],
            cal_freqs[1],
            cal_freqs[2],
            opts.spk_gain,
            opts.repeats,
            opts.pre_delay_sec,
            opts.gap_sec,
        );
    }
    out_stream.play()?;
    if opts.verbose {
        for i in 0..opts.repeats {
            info_line!("TX burst {}/{}", i + 1, opts.repeats);
        }
    }
    let play_sec = (plan.scheduled.len() as f32 / plan.fs) + 0.25;
    std::thread::sleep(Duration::from_secs_f32(play_sec.max(0.25)));
    drop(out_stream);
    info_line!("Transmit done.");
    Ok(())
}
