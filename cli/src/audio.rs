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

fn sinc(x: f32) -> f32 {
    if x.abs() < 1e-8 {
        1.0
    } else {
        (std::f32::consts::PI * x).sin() / (std::f32::consts::PI * x)
    }
}

pub(crate) fn fir_bandpass(len: usize, f_lo_hz: f32, f_hi_hz: f32, fs: f32) -> Vec<f32> {
    let m = (len - 1) as f32 / 2.0;
    let mut h = vec![0.0f32; len];
    let lo = f_lo_hz / fs;
    let hi = f_hi_hz / fs;
    for (n, hn) in h.iter_mut().enumerate() {
        let k = n as f32 - m;
        let ideal = 2.0 * hi * sinc(2.0 * hi * k) - 2.0 * lo * sinc(2.0 * lo * k);
        let win = 0.54 - 0.46 * (2.0 * std::f32::consts::PI * n as f32 / (len as f32 - 1.0)).cos();
        *hn = ideal * win;
    }
    let sum = h.iter().sum::<f32>().abs().max(1e-12);
    for hn in &mut h {
        *hn /= sum;
    }
    h
}

pub(crate) fn fir_filter(x: &[f32], h: &[f32]) -> Vec<f32> {
    let mut y = vec![0.0f32; x.len()];
    for n in 0..x.len() {
        let mut acc = 0.0f32;
        let kmax = h.len().min(n + 1);
        for k in 0..kmax {
            acc += h[k] * x[n - k];
        }
        y[n] = acc;
    }
    y
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
