// Copyright (c) 2026 Elias S. G. Carotti

use std::error::Error;
use std::path::Path;

use plotters::coord::Shift;
use plotters::prelude::*;
use rustfft::num_complex::Complex32;

fn constellation_bounds(points: &[Complex32]) -> Option<f32> {
    let max_abs = points
        .iter()
        .map(|z| z.re.abs().max(z.im.abs()))
        .fold(0.0f32, f32::max);
    if max_abs <= 0.0 {
        None
    } else {
        Some((1.2 * max_abs).max(1.25))
    }
}

fn draw_constellation_panel(
    area: &DrawingArea<BitMapBackend<'_>, Shift>,
    title: &str,
    points: &[Complex32],
    lim: f32,
) -> Result<(), Box<dyn Error>> {
    let mut chart = ChartBuilder::on(area)
        .caption(title, ("sans-serif", 24).into_font())
        .margin(16)
        .x_label_area_size(40)
        .y_label_area_size(45)
        .build_cartesian_2d(-lim..lim, -lim..lim)?;

    chart
        .configure_mesh()
        .x_desc("In-Phase")
        .y_desc("Quadrature")
        .axis_desc_style(("sans-serif", 18))
        .label_style(("sans-serif", 14))
        .light_line_style(RGBAColor(0, 0, 0, 0.08))
        .bold_line_style(RGBAColor(0, 0, 0, 0.18))
        .draw()?;

    chart.draw_series(std::iter::once(PathElement::new(
        vec![(-lim, 0.0), (lim, 0.0)],
        BLACK.mix(0.25),
    )))?;
    chart.draw_series(std::iter::once(PathElement::new(
        vec![(0.0, -lim), (0.0, lim)],
        BLACK.mix(0.25),
    )))?;
    chart.draw_series(points.iter().map(|z| {
        Circle::new(
            (z.re, z.im),
            4,
            ShapeStyle::from(&RGBColor(33, 145, 140)).stroke_width(1),
        )
    }))?;
    Ok(())
}

pub fn save_constellation_comparison_png(
    path: &Path,
    pre_eq: &[Complex32],
    post_eq: &[Complex32],
) -> Result<(), Box<dyn Error>> {
    if pre_eq.is_empty() && post_eq.is_empty() {
        return Ok(());
    }

    let lim = constellation_bounds(pre_eq)
        .into_iter()
        .chain(constellation_bounds(post_eq))
        .fold(1.25f32, f32::max);

    let root = BitMapBackend::new(path, (1280, 720)).into_drawing_area();
    root.fill(&RGBColor(245, 245, 240))?;
    let (left, right) = root.split_horizontally(640);
    draw_constellation_panel(&left, "Constellation: Pre-EQ", pre_eq, lim)?;
    draw_constellation_panel(&right, "Constellation: Post-EQ", post_eq, lim)?;
    root.present()?;
    Ok(())
}

// vim: set ts=4 sw=4 et:
