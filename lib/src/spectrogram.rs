// Copyright (c) 2026 Elias S. G. Carotti

use std::error::Error;
use std::path::Path;

use plotters::prelude::*;
use rustfft::{num_complex::Complex32, FftPlanner};

fn hann_window(n: usize) -> Vec<f32> {
    if n <= 1 {
        return vec![1.0; n.max(1)];
    }
    (0..n)
        .map(|i| 0.5 - 0.5 * (2.0 * std::f32::consts::PI * (i as f32) / ((n - 1) as f32)).cos())
        .collect()
}

fn viridis_like(t: f32) -> RGBColor {
    let t = t.clamp(0.0, 1.0);
    let anchors = [
        (0.0, (13.0, 8.0, 135.0)),
        (0.25, (59.0, 82.0, 139.0)),
        (0.5, (33.0, 145.0, 140.0)),
        (0.75, (94.0, 201.0, 98.0)),
        (1.0, (253.0, 231.0, 37.0)),
    ];
    for w in anchors.windows(2) {
        let (t0, c0) = w[0];
        let (t1, c1) = w[1];
        if t <= t1 {
            let a = ((t - t0) / (t1 - t0)).clamp(0.0, 1.0);
            let lerp = |x0: f32, x1: f32| (x0 + a * (x1 - x0)).round() as u8;
            return RGBColor(
                lerp(c0.0, c1.0),
                lerp(c0.1, c1.1),
                lerp(c0.2, c1.2),
            );
        }
    }
    RGBColor(253, 231, 37)
}

/// Saves a human-readable spectrogram PNG for a mono waveform.
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
    let t_max = ((n_frames - 1) * hop + nfft) as f32 / fs;
    let f_max_khz = fs / 2000.0;

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

    let root = BitMapBackend::new(path, (1280, 720)).into_drawing_area();
    root.fill(&RGBColor(245, 245, 240))?;
    let (main_area, legend_area) = root.split_horizontally(1180);

    let mut chart = ChartBuilder::on(&main_area)
        .caption("RX Spectrogram", ("sans-serif", 28).into_font())
        .margin(20)
        .x_label_area_size(45)
        .y_label_area_size(60)
        .build_cartesian_2d(0f32..t_max, 0f32..f_max_khz)?;

    chart
        .configure_mesh()
        .x_desc("Time [s]")
        .y_desc("Frequency [kHz]")
        .axis_desc_style(("sans-serif", 20))
        .label_style(("sans-serif", 16))
        .light_line_style(RGBAColor(0, 0, 0, 0.08))
        .bold_line_style(RGBAColor(0, 0, 0, 0.18))
        .draw()?;

    let dt = hop as f32 / fs;
    let df_khz = fs / (nfft as f32) / 1000.0;
    let mut cells = Vec::with_capacity(n_frames * n_bins);
    for t in 0..n_frames {
        for k in 0..n_bins {
            let db = spec[t * n_bins + k];
            let norm = ((db - floor_db) / (max_db - floor_db).max(1e-6)).clamp(0.0, 1.0);
            let color = viridis_like(norm).filled();
            let x0 = t as f32 * dt;
            let x1 = (t as f32 + 1.0) * dt;
            let y0 = k as f32 * df_khz;
            let y1 = (k as f32 + 1.0) * df_khz;
            cells.push(Rectangle::new([(x0, y0), (x1, y1)], color));
        }
    }
    chart.draw_series(cells)?;

    legend_area.fill(&RGBColor(245, 245, 240))?;
    let mut legend = ChartBuilder::on(&legend_area)
        .margin(30)
        .y_label_area_size(50)
        .build_cartesian_2d(0f32..1f32, floor_db..max_db)?;
    legend
        .configure_mesh()
        .disable_x_mesh()
        .disable_x_axis()
        .y_desc("Power [dB]")
        .axis_desc_style(("sans-serif", 18))
        .label_style(("sans-serif", 14))
        .draw()?;
    let n_steps = 256usize;
    legend.draw_series((0..n_steps).map(|i| {
        let y0 = floor_db + (i as f32) * (max_db - floor_db) / (n_steps as f32);
        let y1 = floor_db + ((i + 1) as f32) * (max_db - floor_db) / (n_steps as f32);
        Rectangle::new([(0.0, y0), (1.0, y1)], viridis_like((i as f32) / ((n_steps - 1) as f32)).filled())
    }))?;

    root.present()?;
    Ok(())
}

// vim: set ts=4 sw=4 et:
