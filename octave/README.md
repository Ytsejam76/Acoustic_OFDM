# Octave Experiments

This directory contains the fast simulation path for the acoustic OFDM modem.
It is the preferred place to try equalizer ideas before moving them into the
Rust live pipeline.

## Why the Octave path exists

The main purpose of the Octave code is to let us explore DSP assumptions one at
a time under controlled simulated channels.

In the current equalizer work, the main question is not only whether a given
denoiser works, but also which assumptions it adds:

- per-symbol pilot residual denoising in the delay domain
- temporal smoothing across OFDM symbols
- explicit model-order selection such as MDL
- white-noise versus colored-noise observation models

The Octave implementation is used to test those assumptions in isolation before
they are translated into the Rust modem.

## Current equalizer modes

The decoder parameter struct accepts `p.equalizer_mode` with these values:

- `training-pilot`
- `pilot-denoise`
- `pilot-denoise-mdl`
- `pilot-denoise-temporal`
- `pilot-denoise-wiener`
- `pilot-denoise-wiener-psd`
- `pilot-denoise-wiener-shrink`

## Current design rationale

### `pilot-denoise`

This is the baseline delay-domain residual denoiser:

1. equalize the current symbol with the training-based baseline
2. estimate pilot residual observations
3. interpolate the residual over used bins
4. transform to delay taps
5. denoise the taps
6. reconstruct the residual correction

### `pilot-denoise-mdl`

This adds MDL-based delay-domain model-order selection. It was implemented as
an explicit experiment, but recent results suggest that per-symbol MDL is too
opinionated for the current open-air setting and can underperform.

### `pilot-denoise-temporal`

This keeps the per-symbol denoiser and adds causal temporal fusion across
symbols. It is useful as a practical temporal baseline.

### `pilot-denoise-wiener`

This is the current clean temporal Wiener experiment:

- keep the full delay-tap estimate
- do not apply MDL truncation
- do not apply per-symbol tap shrinkage before Wiener smoothing
- apply a temporal Wiener smoother across recent OFDM symbols

The point of this mode is to test a more agnostic temporal MMSE-style baseline
without imposing explicit low-order delay structure.

### `pilot-denoise-wiener-psd`

This is a separate frequency-domain Wiener experiment:

- use all unused FFT bins as measurements of residual non-signal energy
- smooth that disturbance estimate across frequency and time
- interpolate it onto the active bins
- map it into the equalized domain with `|H[k]|^2`
- apply a Wiener-like gain on the equalized subcarriers

This mode is the current baseline for "use all the bins" thinking. It avoids
delay-tap model assumptions and instead uses the unused spectrum as a
disturbance sensor field.

### `pilot-denoise-wiener-shrink`

This is the older hybrid experiment:

- apply the per-symbol tap shrinker first
- then apply temporal Wiener smoothing

It is kept only for comparison. It is not the clean MMSE baseline because the
nonlinear shrinkage step changes the interpretation of the Wiener stage.

## Plotting the delay-domain impulse response

When the equalizer reaches the residual-tap path, `ofdm_test_channel` can plot
the residual delay taps before and after denoising in the figure:

- `impulse_response_denoise_compare`

This figure is useful for inspecting what the denoiser is doing in the
delay-domain basis.

## Running one experiment

From the repository root:

```bash
octave --quiet --eval "addpath('octave'); p=struct(); p.equalizer_mode='pilot-denoise-wiener-psd'; p.modulation='QPSK'; p.use_pilots=true; p.used_bins=[2 3 4 5 6 7 8 9]; p.pilot_bins=[2 4 7 9]; p.pause_before_exit=false; p.make_plots=true; p.save_images=true; p.out_dir='output/octave_eq'; p.save_decoder_constellation=true; p.verbose=true; ofdm_test_channel(p);"
```

## Running the equalizer comparison sweep

Use the dedicated sweep helper:

```bash
octave --quiet --eval "addpath('octave'); cfg=struct(); cfg.num_trials=20; cfg.show_progress=true; cfg.make_octave_plots=true; cfg.save_plot=true; cfg.out_dir='output/octave_eq'; cfg.plot_filename='equalizer_sweep.png'; cfg.oracle_sync=true; cfg.equalizer_modes={'pilot-denoise','pilot-denoise-temporal','pilot-denoise-wiener','pilot-denoise-wiener-psd'}; cfg.base_params=struct(); cfg.base_params.modulation='QPSK'; cfg.base_params.use_pilots=true; cfg.base_params.used_bins=[2 3 4 5 6 7 8 9]; cfg.base_params.pilot_bins=[2 4 7 9]; ofdm_equalizer_sweep(cfg);"
```

This saves:

- `output/octave_eq/equalizer_sweep.png`
- `output/octave_eq/equalizer_sweep_stats.mat`

## Running against recorded interference

Recorded ambient audio can be mixed into the passband waveform offline through
`p.interference_file`. The intended local folder is:

- `octave/noise_recordings/`

That directory is ignored by Git so local captures never get committed.

For one run with a fixed chunk:

```bash
octave --quiet --eval "addpath('octave'); p=struct(); p.equalizer_mode='pilot-denoise-wiener-psd'; p.interference_file='octave/noise_recordings/interference.wav'; p.interference_ratio_db=12; p.interference_offset_seconds=1.5; p.make_plots=true; p.save_images=false; p.pause_before_exit=false; ofdm_test_channel(p);"
```

For a repeatable chunk sweep across modes:

```bash
octave --quiet --eval "addpath('octave'); cfg=struct(); cfg.num_chunks=8; cfg.interference_file='octave/noise_recordings/interference.wav'; cfg.interference_ratio_db=12; cfg.chunk_advance_samples=12000; cfg.equalizer_modes={'pilot-denoise','pilot-denoise-temporal','pilot-denoise-wiener-psd'}; ofdm_recorded_interference_sweep(cfg);"
```

This uses the same chunk offsets for every mode so the comparison stays fair.

For saved plots on the same recorded segments across multiple modes:

```bash
RUN_KIND=compare \
  NUM_SEGMENTS=3 \
  CHUNK_ADVANCE_SAMPLES=12000 \
  SAVE_IMAGES=true \
  OUT_DIR=output/octave_recorded_compare \
  bash scripts/octave_recorded_interference.sh
```

This saves one directory per `(segment, mode)` pair and also writes a
`summary.tsv` file in the chosen output root.

## Current recommendation

For exploratory work, start with:

1. `pilot-denoise`
2. `pilot-denoise-temporal`
3. `pilot-denoise-wiener-psd`
4. `pilot-denoise-wiener`

Only use `pilot-denoise-mdl` or `pilot-denoise-wiener-shrink` as comparison
experiments, not as default directions.
