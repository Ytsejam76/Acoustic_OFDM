// Copyright (c) 2026 Elias S. G. Carotti

use std::error::Error;
use std::path::Path;

use image::{GrayImage, Luma};
use rustfft::{num_complex::Complex32, FftPlanner};

fn hann_window(n: usize) -> Vec<f32> {
    if n <= 1 {
        return vec![1.0; n.max(1)];
    }
    (0..n)
        .map(|i| 0.5 - 0.5 * (2.0 * std::f32::consts::PI * (i as f32) / ((n - 1) as f32)).cos())
        .collect()
}

/// Saves a simple grayscale spectrogram PNG for a mono waveform.
///
/// Parameters:
/// - `path`: output PNG path.
/// - `x`: input mono waveform.
/// - `fs`: sample rate in Hz.
/// Returns:
/// - `Result<(), Box<dyn Error>>`: `Ok(())` when the PNG is written.
pub fn save_spectrogram_png(path: &Path, x: &[f32], fs: f32) -> Result<(), Box<dyn Error>> {
    if x.is_empty() {
        return Ok(());
    }
    let nfft = 512usize;
    let hop = 128usize;
    if x.len() < nfft {
        return Ok(());
    }
    let window = hann_window(nfft);
    let mut planner = FftPlanner::<f32>::new();
    let fft = planner.plan_fft_forward(nfft);
    let n_frames = 1 + (x.len() - nfft) / hop;
    let n_bins = nfft / 2;
    let mut spec = vec![0.0f32; n_frames * n_bins];
    let mut frame = vec![Complex32::new(0.0, 0.0); nfft];
    let mut max_db = -120.0f32;
    let mut min_db = 0.0f32;
    for t in 0..n_frames {
        let s0 = t * hop;
        for i in 0..nfft {
            frame[i] = Complex32::new(x[s0 + i] * window[i], 0.0);
        }
        fft.process(&mut frame);
        for k in 0..n_bins {
            let p = frame[k].norm_sqr() / (nfft as f32);
            let db = 10.0 * p.max(1e-12).log10();
            spec[t * n_bins + k] = db;
            max_db = max_db.max(db);
            min_db = min_db.min(db);
        }
    }
    let floor_db = (max_db - 70.0).max(min_db);
    let width = n_frames as u32;
    let height = n_bins as u32;
    let mut img = GrayImage::new(width, height);
    for t in 0..n_frames {
        for k in 0..n_bins {
            let db = spec[t * n_bins + k];
            let norm = ((db - floor_db) / (max_db - floor_db).max(1e-6)).clamp(0.0, 1.0);
            let y = (n_bins - 1 - k) as u32;
            img.put_pixel(t as u32, y, Luma([(255.0 * norm) as u8]));
        }
    }
    img.save(path)?;
    let _ = fs;
    Ok(())
}

// vim: set ts=4 sw=4 et:
