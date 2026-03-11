#!/usr/bin/env python3
"""Render BER/PER vs SNR plots from Octave stats using seaborn."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy.io import loadmat
from scipy.interpolate import PchipInterpolator


def as_array(x):
    return np.asarray(x, dtype=float).reshape(-1)


def smooth_xy(x, y, npts=240):
    x = as_array(x)
    y = as_array(y)
    m = ~np.isnan(x) & ~np.isnan(y)
    x = x[m]
    y = y[m]
    if x.size < 3:
        return x, y
    order = np.argsort(x)
    x = x[order]
    y = y[order]
    xu, idx = np.unique(x, return_index=True)
    yu = y[idx]
    if xu.size < 3:
        return xu, yu
    xd = np.linspace(float(xu.min()), float(xu.max()), int(npts))
    yd = PchipInterpolator(xu, yu)(xd)
    return xd, yd


def parse_stats(path: Path):
    data = loadmat(path, simplify_cells=True)
    stats = data["stats"]
    if "oracle" in stats and "nonoracle" in stats:
        return stats["oracle"], stats["nonoracle"], True
    return stats, None, False


def apply_theme():
    sns.set_theme(
        style="ticks",
        context="talk",
        palette="colorblind",
        rc={
            "figure.facecolor": "#f6f7fb",
            "axes.facecolor": "#fbfcff",
            "axes.edgecolor": "#3a3f47",
            "axes.grid": True,
            "grid.color": "#c7ced9",
            "grid.linestyle": "--",
            "grid.alpha": 0.6,
            "axes.titleweight": "bold",
            "axes.titlepad": 10.0,
            "axes.labelcolor": "#1f2430",
            "xtick.color": "#2b313d",
            "ytick.color": "#2b313d",
            "legend.frameon": True,
            "legend.framealpha": 0.9,
            "legend.facecolor": "#ffffff",
        },
    )


def plot_compare(oracle, est, out_png: Path, smooth: bool = True):
    apply_theme()
    fig, axes = plt.subplots(2, 1, figsize=(11, 9), sharex=True, constrained_layout=True)

    x1 = as_array(oracle["snr_db"])
    per1 = as_array(oracle["per"])
    ber1 = as_array(oracle["ber_decoded"])

    x2 = as_array(est["snr_db"])
    per2 = as_array(est["per"])
    ber2 = as_array(est["ber_decoded"])

    if smooth:
        x1s, per1s = smooth_xy(x1, per1)
        x2s, per2s = smooth_xy(x2, per2)
        sns.lineplot(x=x1s, y=np.clip(per1s, 0, 1), linewidth=2.5, label="Oracle sync (smooth)", ax=axes[0])
        sns.lineplot(x=x2s, y=np.clip(per2s, 0, 1), linewidth=2.5, linestyle="--", label="Estimated sync (smooth)", ax=axes[0])
    else:
        sns.lineplot(x=x1, y=per1, marker="o", linewidth=2.0, label="Oracle sync", ax=axes[0])
        sns.lineplot(x=x2, y=per2, marker="s", linewidth=2.0, linestyle="--", label="Estimated sync", ax=axes[0])
    sns.scatterplot(x=x1, y=per1, s=36, marker="o", color=axes[0].lines[0].get_color(), ax=axes[0], legend=False)
    sns.scatterplot(x=x2, y=per2, s=36, marker="s", color=axes[0].lines[1].get_color(), ax=axes[0], legend=False)
    axes[0].set_title("Packet Error Rate vs SNR")
    axes[0].set_xlabel("SNR (dB)")
    axes[0].set_ylabel("PER")
    axes[0].set_ylim(0, 1)
    axes[0].legend(loc="upper right", frameon=True)

    m1 = ~np.isnan(ber1)
    m2 = ~np.isnan(ber2)
    has_positive = (np.any(m1 & (ber1 > 0)) or np.any(m2 & (ber2 > 0)))
    if np.any(m1):
        if smooth:
            x1b, y1b = smooth_xy(x1[m1], ber1[m1])
            sns.lineplot(x=x1b, y=np.maximum(y1b, 0), linewidth=2.5, label="Oracle sync (smooth)", ax=axes[1])
        else:
            sns.lineplot(x=x1[m1], y=np.maximum(ber1[m1], 0), marker="o", linewidth=2.0, label="Oracle sync", ax=axes[1])
        sns.scatterplot(x=x1[m1], y=ber1[m1], s=36, marker="o", color=axes[1].lines[-1].get_color(), ax=axes[1], legend=False)
    if np.any(m2):
        if smooth:
            x2b, y2b = smooth_xy(x2[m2], ber2[m2])
            sns.lineplot(x=x2b, y=np.maximum(y2b, 0), linewidth=2.5, linestyle="--", label="Estimated sync (smooth)", ax=axes[1])
        else:
            sns.lineplot(x=x2[m2], y=np.maximum(ber2[m2], 0), marker="s", linewidth=2.0, linestyle="--", label="Estimated sync", ax=axes[1])
        sns.scatterplot(x=x2[m2], y=ber2[m2], s=36, marker="s", color=axes[1].lines[-1].get_color(), ax=axes[1], legend=False)
    if has_positive:
        axes[1].set_yscale("log")
    else:
        axes[1].set_ylim(0, 0.05)
    axes[1].set_title("BER on Decoded Bits vs SNR")
    axes[1].set_xlabel("SNR (dB)")
    axes[1].set_ylabel("BER")
    axes[1].legend(loc="upper right", frameon=True)

    fig.savefig(out_png, dpi=180)
    plt.close(fig)


def plot_single(stats, out_png: Path, smooth: bool = True):
    apply_theme()
    fig, axes = plt.subplots(2, 1, figsize=(11, 9), sharex=True, constrained_layout=True)

    x = as_array(stats["snr_db"])
    per = as_array(stats["per"])
    ber = as_array(stats["ber_decoded"])
    dec = as_array(stats["decode_rate"])

    if smooth:
        xs1, ys1 = smooth_xy(x, per)
        xs2, ys2 = smooth_xy(x, 1 - dec)
        sns.lineplot(x=xs1, y=np.clip(ys1, 0, 1), linewidth=2.5, label="PER (smooth)", ax=axes[0])
        sns.lineplot(x=xs2, y=np.clip(ys2, 0, 1), linewidth=2.5, linestyle="--", label="1 - Decode rate (smooth)", ax=axes[0])
    else:
        sns.lineplot(x=x, y=per, marker="o", linewidth=2.0, label="PER", ax=axes[0])
        sns.lineplot(x=x, y=1 - dec, marker="s", linewidth=2.0, linestyle="--", label="1 - Decode rate", ax=axes[0])
    sns.scatterplot(x=x, y=per, s=36, marker="o", color=axes[0].lines[0].get_color(), ax=axes[0], legend=False)
    sns.scatterplot(x=x, y=1 - dec, s=36, marker="s", color=axes[0].lines[1].get_color(), ax=axes[0], legend=False)
    axes[0].set_title("Packet Error and Erasure Rates vs SNR")
    axes[0].set_xlabel("SNR (dB)")
    axes[0].set_ylabel("Rate")
    axes[0].set_ylim(0, 1)
    axes[0].legend(loc="upper right", frameon=True)

    m = ~np.isnan(ber)
    has_positive = np.any(m & (ber > 0))
    if np.any(m):
        if smooth:
            xsb, ysb = smooth_xy(x[m], ber[m])
            sns.lineplot(x=xsb, y=np.maximum(ysb, 0), linewidth=2.5, label="BER (smooth)", ax=axes[1])
        else:
            sns.lineplot(x=x[m], y=np.maximum(ber[m], 0), marker="o", linewidth=2.0, label="BER", ax=axes[1])
        sns.scatterplot(x=x[m], y=ber[m], s=36, marker="o", color=axes[1].lines[-1].get_color(), ax=axes[1], legend=False)
    if has_positive:
        axes[1].set_yscale("log")
    else:
        axes[1].set_ylim(0, 0.05)
    axes[1].set_title("BER on Decoded Bits vs SNR")
    axes[1].set_xlabel("SNR (dB)")
    axes[1].set_ylabel("BER")

    fig.savefig(out_png, dpi=180)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stats", required=True, help="Path to snr_sweep_stats.mat")
    ap.add_argument("--out", required=True, help="Output PNG path")
    ap.add_argument("--no-smooth", action="store_true", help="Disable curve smoothing")
    args = ap.parse_args()

    stats_path = Path(args.stats)
    out_png = Path(args.out)
    out_png.parent.mkdir(parents=True, exist_ok=True)

    s1, s2, is_compare = parse_stats(stats_path)
    if is_compare:
        plot_compare(s1, s2, out_png, smooth=not args.no_smooth)
    else:
        plot_single(s1, out_png, smooth=not args.no_smooth)
    print(f"Saved seaborn plot to: {out_png}")


if __name__ == "__main__":
    main()
