// Copyright (c) 2026 Elias S. G. Carotti

use std::error::Error;
use std::path::Path;
use std::sync::{Arc, Mutex};
use std::time::{Duration, Instant};

use acoustic_ofdm::{
    save_spectrogram_png, save_spectrogram_png_with_options, save_wav_mono_i16, OfdmConfig,
};
use cpal::traits::{DeviceTrait, HostTrait, StreamTrait};
use ringbuf::{traits::*, HeapRb};

use crate::audio::{build_input_stream, clipping_diag, signal_diag};
use crate::cli_args::{AudioOpts, ORACLE_PAYLOAD};
use crate::{debug_line, info_line, warn_line};

pub(crate) fn cmd_rx(
    cfg: &OfdmConfig,
    opts: &AudioOpts,
    _stdout_raw: bool,
) -> Result<(), Box<dyn Error>> {
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

    info_line!("Input device : {}", in_dev.name()?);
    info_line!(
        "Stream config: {} Hz, {} ch, in {:?}",
        in_cfg.sample_rate().0,
        in_cfg.channels(),
        in_cfg.sample_format()
    );
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
    let deadline = Instant::now() + Duration::from_secs_f32(opts.duration_sec);
    while Instant::now() < deadline {
        while let Some(s) = in_cons.try_pop() {
            rx.push(s);
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

    if opts.verbose {
        let (rms, peak, first_loud) = signal_diag(&rx);
        let (clipped, clipped_frac) = clipping_diag(&rx);
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

    Ok(())
}
