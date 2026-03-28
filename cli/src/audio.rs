// Copyright (c) 2026 Elias S. G. Carotti

use std::error::Error;
use std::sync::{Arc, Mutex};

use cpal::traits::DeviceTrait;
use cpal::{SampleFormat, Stream, StreamConfig};
use ringbuf::traits::Producer;
use ringbuf::HeapProd;

pub(crate) fn i16_to_f32(x: i16) -> f32 {
    (x as f32) / 32768.0
}
pub(crate) fn u16_to_f32(x: u16) -> f32 {
    ((x as f32) - 32768.0) / 32768.0
}
pub(crate) fn f32_to_i16(x: f32) -> i16 {
    (x.clamp(-1.0, 1.0) * 32767.0).round() as i16
}
pub(crate) fn f32_to_u16(x: f32) -> u16 {
    ((x.clamp(-1.0, 1.0) * 32767.0) + 32768.0).round() as u16
}

pub(crate) fn signal_diag(x: &[f32]) -> (f32, f32, usize) {
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

pub(crate) fn clipping_diag(x: &[f32]) -> (usize, f32) {
    if x.is_empty() {
        return (0, 0.0);
    }
    let clipped = x.iter().filter(|&&s| s.abs() >= 0.995).count();
    (clipped, (clipped as f32) / (x.len() as f32))
}

pub(crate) fn build_input_stream(
    dev: &cpal::Device,
    config: &StreamConfig,
    fmt: SampleFormat,
    gain: f32,
    prod: Arc<Mutex<HeapProd<f32>>>,
) -> Result<Stream, Box<dyn Error>> {
    let err_fn = |e| eprintln!("input stream error: {e}");
    let channels = usize::from(config.channels.max(1));
    let stream = match fmt {
        SampleFormat::I16 => dev.build_input_stream(
            config,
            move |data: &[i16], _| {
                if let Ok(mut p) = prod.lock() {
                    for frame in data.chunks_exact(channels) {
                        let mono = frame.iter().map(|&s| i16_to_f32(s)).sum::<f32>() / (channels as f32);
                        let _ = p.try_push(mono * gain);
                    }
                }
            },
            err_fn,
            None,
        )?,
        SampleFormat::U16 => dev.build_input_stream(
            config,
            move |data: &[u16], _| {
                if let Ok(mut p) = prod.lock() {
                    for frame in data.chunks_exact(channels) {
                        let mono = frame.iter().map(|&s| u16_to_f32(s)).sum::<f32>() / (channels as f32);
                        let _ = p.try_push(mono * gain);
                    }
                }
            },
            err_fn,
            None,
        )?,
        SampleFormat::F32 => dev.build_input_stream(
            config,
            move |data: &[f32], _| {
                if let Ok(mut p) = prod.lock() {
                    for frame in data.chunks_exact(channels) {
                        let mono = frame.iter().copied().sum::<f32>() / (channels as f32);
                        let _ = p.try_push(mono * gain);
                    }
                }
            },
            err_fn,
            None,
        )?,
        _ => return Err("unsupported input sample format".into()),
    };
    Ok(stream)
}

pub(crate) fn build_output_stream(
    dev: &cpal::Device,
    config: &StreamConfig,
    fmt: SampleFormat,
    mut tx: Vec<f32>,
) -> Result<Stream, Box<dyn Error>> {
    let err_fn = |e| eprintln!("output stream error: {e}");
    let channels = usize::from(config.channels.max(1));
    let mut idx = 0usize;
    if tx.is_empty() {
        tx.push(0.0);
    }
    let stream = match fmt {
        SampleFormat::I16 => dev.build_output_stream(
            config,
            move |data: &mut [i16], _| {
                for frame in data.chunks_exact_mut(channels) {
                    let v = if idx < tx.len() { tx[idx] } else { 0.0 };
                    let y = f32_to_i16(v);
                    for s in frame {
                        *s = y;
                    }
                    idx += 1;
                }
            },
            err_fn,
            None,
        )?,
        SampleFormat::U16 => dev.build_output_stream(
            config,
            move |data: &mut [u16], _| {
                for frame in data.chunks_exact_mut(channels) {
                    let v = if idx < tx.len() { tx[idx] } else { 0.0 };
                    let y = f32_to_u16(v);
                    for s in frame {
                        *s = y;
                    }
                    idx += 1;
                }
            },
            err_fn,
            None,
        )?,
        SampleFormat::F32 => dev.build_output_stream(
            config,
            move |data: &mut [f32], _| {
                for frame in data.chunks_exact_mut(channels) {
                    let y = if idx < tx.len() { tx[idx] } else { 0.0 };
                    for s in frame {
                        *s = y;
                    }
                    idx += 1;
                }
            },
            err_fn,
            None,
        )?,
        _ => return Err("unsupported output sample format".into()),
    };
    Ok(stream)
}
