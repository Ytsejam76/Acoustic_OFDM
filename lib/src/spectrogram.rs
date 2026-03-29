// Copyright (c) 2026 Elias S. G. Carotti

use std::error::Error;
use std::path::Path;

use plotters::prelude::*;
use rustfft::{num_complex::Complex32, FftPlanner};

/// Analysis window used for STFT-based spectrogram generation.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SpectrogramWindow {
    /// Hann window, good general-purpose sidelobe suppression.
    ///
    /// Reference: [Wikipedia: Hann window](https://en.wikipedia.org/wiki/Window_function#Hann_window)
    Hann,
    /// Hamming window, slightly narrower main lobe than Hann.
    ///
    /// Reference: [Wikipedia: Hamming window](https://en.wikipedia.org/wiki/Window_function#Hamming_window)
    Hamming,
    /// Blackman window, stronger sidelobe suppression at the cost of width.
    ///
    /// Reference: [Wikipedia: Blackman window](https://en.wikipedia.org/wiki/Window_function#Blackman_window)
    Blackman,
    /// Rectangular window, no tapering.
    ///
    /// Reference: [Wikipedia: Rectangular window](https://en.wikipedia.org/wiki/Window_function#Rectangular_window)
    Rect,
}

/// Configuration for spectrogram rendering.
#[derive(Clone, Copy, Debug)]
pub struct SpectrogramOptions {
    /// FFT size used for each STFT frame.
    pub nfft: usize,
    /// Hop size, in samples, between adjacent STFT frames.
    pub hop: usize,
    /// Analysis window applied before each FFT.
    pub window: SpectrogramWindow,
}

impl Default for SpectrogramOptions {
    fn default() -> Self {
        Self {
            nfft: 512,
            hop: 128,
            window: SpectrogramWindow::Hann,
        }
    }
}

fn window_samples(n: usize, kind: SpectrogramWindow) -> Vec<f32> {
    if n <= 1 {
        return vec![1.0; n.max(1)];
    }
    match kind {
        SpectrogramWindow::Rect => vec![1.0; n],
        SpectrogramWindow::Hann => (0..n)
            .map(|i| 0.5 - 0.5 * (2.0 * std::f32::consts::PI * (i as f32) / ((n - 1) as f32)).cos())
            .collect(),
        SpectrogramWindow::Hamming => (0..n)
            .map(|i| {
                0.54 - 0.46 * (2.0 * std::f32::consts::PI * (i as f32) / ((n - 1) as f32)).cos()
            })
            .collect(),
        SpectrogramWindow::Blackman => (0..n)
            .map(|i| {
                let phi = 2.0 * std::f32::consts::PI * (i as f32) / ((n - 1) as f32);
                0.42 - 0.5 * phi.cos() + 0.08 * (2.0 * phi).cos()
            })
            .collect(),
    }
}

fn jet_like(t: f32) -> RGBColor {
    let t = t.clamp(0.0, 1.0);
    let anchors = [
        (0.0, (0.0, 0.0, 131.0)),
        (0.125, (0.0, 60.0, 170.0)),
        (0.375, (5.0, 255.0, 255.0)),
        (0.625, (255.0, 255.0, 0.0)),
        (0.875, (250.0, 0.0, 0.0)),
        (1.0, (128.0, 0.0, 0.0)),
    ];
    for w in anchors.windows(2) {
        let (t0, c0) = w[0];
        let (t1, c1) = w[1];
        if t <= t1 {
            let a = ((t - t0) / (t1 - t0)).clamp(0.0, 1.0);
            let lerp = |x0: f32, x1: f32| (x0 + a * (x1 - x0)).round() as u8;
            return RGBColor(lerp(c0.0, c1.0), lerp(c0.1, c1.1), lerp(c0.2, c1.2));
        }
    }
    RGBColor(128, 0, 0)
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
    save_spectrogram_png_with_options(path, x, fs, SpectrogramOptions::default())
}

pub fn save_spectrogram_png_with_options(
    path: &Path,
    x: &[f32],
    fs: f32,
    opts: SpectrogramOptions,
) -> Result<(), Box<dyn Error>> {
    if x.is_empty() {
        return Ok(());
    }
    let nfft = opts.nfft.max(16);
    let hop = opts.hop.max(1).min(nfft);
    if x.len() < nfft {
        return Ok(());
    }

    let window = window_samples(nfft, opts.window);
    let mut planner = FftPlanner::<f32>::new();
    let fft = planner.plan_fft_forward(nfft);
    let n_frames = 1 + (x.len() - nfft) / hop;
    let n_bins = nfft / 2 + 1;
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
    let (main_area, legend_area) = root.split_horizontally(1120);

    let mut chart = ChartBuilder::on(&main_area)
        .caption("Spectrogram", ("sans-serif", 28).into_font())
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
            let color = jet_like(norm).filled();
            let x0 = t as f32 * dt;
            let x1 = (t as f32 + 1.0) * dt;
            let y0 = k as f32 * df_khz;
            let y1 = (k as f32 + 1.0) * df_khz;
            cells.push(Rectangle::new([(x0, y0), (x1, y1)], color));
        }
    }
    chart.draw_series(cells)?;

    legend_area.fill(&RGBColor(245, 245, 240))?;
    let (_lw, lh) = legend_area.dim_in_pixel();
    let lh = lh as i32;
    let bar_left = 28i32;
    let bar_right = 52i32;
    let bar_top = 36i32;
    let bar_bottom = (lh - 42).max(bar_top + 1);
    let n_steps = 256usize;
    for i in 0..n_steps {
        let y0 = bar_bottom - ((i as i32) * (bar_bottom - bar_top) / (n_steps as i32));
        let y1 = bar_bottom - (((i + 1) as i32) * (bar_bottom - bar_top) / (n_steps as i32));
        legend_area.draw(&Rectangle::new(
            [(bar_left, y1), (bar_right, y0.max(y1 + 1))],
            jet_like((i as f32) / ((n_steps - 1) as f32)).filled(),
        ))?;
    }
    legend_area.draw(&Rectangle::new(
        [(bar_left, bar_top), (bar_right, bar_bottom)],
        ShapeStyle::from(&BLACK).stroke_width(1),
    ))?;

    let tick_vals = [
        max_db,
        max_db - 20.0,
        max_db - 40.0,
        max_db - 60.0,
        floor_db,
    ];
    for &db in &tick_vals {
        let a = ((db - floor_db) / (max_db - floor_db).max(1e-6)).clamp(0.0, 1.0);
        let y = bar_bottom - ((a * ((bar_bottom - bar_top) as f32)).round() as i32);
        legend_area.draw(&PathElement::new(
            vec![(bar_right + 2, y), (bar_right + 10, y)],
            BLACK,
        ))?;
        legend_area.draw(&Text::new(
            format!("{db:.0}"),
            (bar_right + 14, y + 5),
            ("sans-serif", 15).into_font(),
        ))?;
    }
    legend_area.draw(&Text::new(
        "Power [dB]",
        (bar_right + 42, lh / 2),
        ("sans-serif", 18)
            .into_font()
            .transform(FontTransform::Rotate90),
    ))?;

    root.present()?;
    Ok(())
}

// vim: set ts=4 sw=4 et:
