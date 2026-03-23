// Copyright (c) 2026 Elias S. G. Carotti

use std::error::Error;
use std::time::Duration;

use acoustic_ofdm::{encode_single_packet_passband, OfdmConfig};
use cpal::traits::{DeviceTrait, HostTrait, StreamTrait};

use crate::audio::build_output_stream;
use crate::cli_args::AudioOpts;
use crate::info_line;

pub(crate) fn cmd_tx(payload: &[u8], cfg: &OfdmConfig, opts: &AudioOpts) -> Result<(), Box<dyn Error>> {
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
    let mut tx_samples = Vec::with_capacity(total_n);
    tx_samples.extend(std::iter::repeat_n(0.0, pre_n));
    for rep in 0..opts.repeats {
        tx_samples.extend_from_slice(&tx);
        if rep + 1 < opts.repeats {
            tx_samples.extend(std::iter::repeat_n(0.0, gap_n));
        }
    }

    let out_stream = build_output_stream(
        &out_dev,
        &out_cfg.clone().into(),
        out_cfg.sample_format(),
        tx_samples,
    )?;

    info_line!("Output device: {}", out_dev.name()?);
    info_line!(
        "Stream config: {} Hz, out {:?}",
        out_cfg.sample_rate().0,
        out_cfg.sample_format()
    );
    info_line!("Wake preamble: {}", cfg_rt.wake_preamble.as_str());
    if opts.oracle {
        info_line!("Oracle mode: enabled ({} bytes)", payload.len());
    }
    info_line!("Transmit samples: {}", tx.len());
    if opts.verbose {
        let peak = tx.iter().fold(0.0f32, |m, &v| if v.abs() > m { v.abs() } else { m });
        info_line!(
            "TX diagnostics: duration={:.3}s peak={peak:.3} spk_gain={:.3} repeats={} pre_delay={:.2}s gap={:.2}s",
            (tx.len() as f32) / cfg_rt.fs,
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
    let play_sec = (total_n as f32 / cfg_rt.fs) + 0.25;
    std::thread::sleep(Duration::from_secs_f32(play_sec.max(0.25)));
    drop(out_stream);
    info_line!("Transmit done.");
    Ok(())
}
