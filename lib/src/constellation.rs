// Copyright (c) 2026 Elias S. G. Carotti

use std::error::Error;
use std::path::Path;

use plotters::coord::Shift;
use plotters::prelude::*;
use rustfft::num_complex::Complex32;

use crate::debug::PassbandChannelCompareDump;

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

fn phase_deg(z: Complex32) -> f32 {
    z.arg().to_degrees()
}

pub fn save_channel_compare_png(
    path: &Path,
    dump: &PassbandChannelCompareDump,
) -> Result<(), Box<dyn Error>> {
    if dump.rows.is_empty() {
        return Ok(());
    }

    let mut bins = dump.rows.iter().map(|r| r.used_bin).collect::<Vec<_>>();
    bins.sort_unstable();
    bins.dedup();

    let mut actual = Vec::with_capacity(bins.len());
    let mut est_train = Vec::with_capacity(bins.len());
    let mut est_pilot = Vec::with_capacity(bins.len());
    for &bin in &bins {
        let rows = dump
            .rows
            .iter()
            .filter(|r| r.used_bin == bin)
            .collect::<Vec<_>>();
        let n = rows.len().max(1) as f32;
        let a = rows
            .iter()
            .fold(Complex32::new(0.0, 0.0), |acc, r| acc + r.actual_h)
            / n;
        let t = rows
            .iter()
            .fold(Complex32::new(0.0, 0.0), |acc, r| acc + r.estimated_h_train)
            / n;
        let p = rows
            .iter()
            .fold(Complex32::new(0.0, 0.0), |acc, r| acc + r.estimated_h_pilot)
            / n;
        actual.push((bin as f32, a));
        est_train.push((bin as f32, t));
        est_pilot.push((bin as f32, p));
    }

    let x_min = *bins.first().unwrap() as f32 - 0.5;
    let x_max = *bins.last().unwrap() as f32 + 0.5;
    let mag_max = actual
        .iter()
        .chain(est_train.iter())
        .chain(est_pilot.iter())
        .map(|(_, z)| z.norm())
        .fold(0.0f32, f32::max)
        .max(1.0);

    let root = BitMapBackend::new(path, (1280, 800)).into_drawing_area();
    root.fill(&RGBColor(245, 245, 240))?;
    let (top, bottom) = root.split_vertically(400);

    {
        let mut chart = ChartBuilder::on(&top)
            .caption("Channel Envelope", ("sans-serif", 24).into_font())
            .margin(16)
            .x_label_area_size(40)
            .y_label_area_size(55)
            .build_cartesian_2d(x_min..x_max, 0.0f32..(1.2 * mag_max))?;
        chart
            .configure_mesh()
            .x_desc("Used Bin")
            .y_desc("|H|")
            .axis_desc_style(("sans-serif", 18))
            .label_style(("sans-serif", 14))
            .draw()?;
        chart
            .draw_series(std::iter::once(PathElement::new(
                actual
                    .iter()
                    .map(|(b, z)| (*b, z.norm()))
                    .collect::<Vec<_>>(),
                RGBColor(214, 39, 40),
            )))?
            .label("Actual")
            .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], RGBColor(214, 39, 40)));
        chart
            .draw_series(std::iter::once(PathElement::new(
                est_train
                    .iter()
                    .map(|(b, z)| (*b, z.norm()))
                    .collect::<Vec<_>>(),
                RGBColor(31, 119, 180),
            )))?
            .label("Training Estimate")
            .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], RGBColor(31, 119, 180)));
        chart
            .draw_series(std::iter::once(PathElement::new(
                est_pilot
                    .iter()
                    .map(|(b, z)| (*b, z.norm()))
                    .collect::<Vec<_>>(),
                RGBColor(44, 160, 44),
            )))?
            .label("Pilot-Adjusted")
            .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], RGBColor(44, 160, 44)));
        chart
            .configure_series_labels()
            .border_style(BLACK.mix(0.3))
            .background_style(WHITE.mix(0.8))
            .draw()?;
    }

    {
        let mut chart = ChartBuilder::on(&bottom)
            .caption("Channel Phase", ("sans-serif", 24).into_font())
            .margin(16)
            .x_label_area_size(40)
            .y_label_area_size(55)
            .build_cartesian_2d(x_min..x_max, -180.0f32..180.0f32)?;
        chart
            .configure_mesh()
            .x_desc("Used Bin")
            .y_desc("Phase (deg)")
            .axis_desc_style(("sans-serif", 18))
            .label_style(("sans-serif", 14))
            .draw()?;
        chart.draw_series(std::iter::once(PathElement::new(
            actual
                .iter()
                .map(|(b, z)| (*b, phase_deg(*z)))
                .collect::<Vec<_>>(),
            RGBColor(214, 39, 40),
        )))?;
        chart.draw_series(std::iter::once(PathElement::new(
            est_train
                .iter()
                .map(|(b, z)| (*b, phase_deg(*z)))
                .collect::<Vec<_>>(),
            RGBColor(31, 119, 180),
        )))?;
        chart.draw_series(std::iter::once(PathElement::new(
            est_pilot
                .iter()
                .map(|(b, z)| (*b, phase_deg(*z)))
                .collect::<Vec<_>>(),
            RGBColor(44, 160, 44),
        )))?;
    }

    root.present()?;
    Ok(())
}

// vim: set ts=4 sw=4 et:
