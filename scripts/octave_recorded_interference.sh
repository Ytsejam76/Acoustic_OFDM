#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

MODE="${MODE:-pilot-denoise-wiener-psd}"
INTERFERENCE_FILE="${INTERFERENCE_FILE:-octave/noise_recordings/interference.wav}"
INTERFERENCE_RATIO_DB="${INTERFERENCE_RATIO_DB:-12}"
MODULATION="${MODULATION:-QPSK}"
USED_BINS="${USED_BINS:-[2 3 4 5 6 7 8 9]}"
PILOT_BINS="${PILOT_BINS:-[2 4 7 9]}"

RUN_KIND="${RUN_KIND:-single}"
OFFSET_SECONDS="${OFFSET_SECONDS:-1.5}"
OFFSET_SAMPLES="${OFFSET_SAMPLES:-}"
MAKE_PLOTS="${MAKE_PLOTS:-true}"
SAVE_IMAGES="${SAVE_IMAGES:-false}"
OUT_DIR="${OUT_DIR:-output/octave_eq}"

NUM_CHUNKS="${NUM_CHUNKS:-8}"
CHUNK_ADVANCE_SAMPLES="${CHUNK_ADVANCE_SAMPLES:-12000}"
NUM_SEGMENTS="${NUM_SEGMENTS:-3}"
START_OFFSET_SAMPLES="${START_OFFSET_SAMPLES:-0}"
SHOW_PROGRESS="${SHOW_PROGRESS:-true}"
DEFAULT_MODES="{'pilot-denoise','pilot-denoise-temporal','pilot-denoise-wiener-psd'}"
MODES="${MODES:-$DEFAULT_MODES}"

if [ ! -f "$INTERFERENCE_FILE" ]; then
  printf 'missing interference file: %s\n' "$INTERFERENCE_FILE" >&2
  exit 1
fi

bool_to_octave() {
  case "${1,,}" in
    1|true|yes|on) printf 'true' ;;
    0|false|no|off) printf 'false' ;;
    *)
      printf 'invalid boolean: %s\n' "$1" >&2
      exit 1
      ;;
  esac
}

MAKE_PLOTS_OCT="$(bool_to_octave "$MAKE_PLOTS")"
SAVE_IMAGES_OCT="$(bool_to_octave "$SAVE_IMAGES")"
SHOW_PROGRESS_OCT="$(bool_to_octave "$SHOW_PROGRESS")"

if [ "$RUN_KIND" = "single" ]; then
  OFFSET_EXPR="p.interference_offset_seconds=${OFFSET_SECONDS};"
  if [ -n "$OFFSET_SAMPLES" ]; then
    OFFSET_EXPR="p.interference_offset_samples=${OFFSET_SAMPLES};"
  fi

  octave --quiet --eval "addpath('octave'); p=struct(); p.equalizer_mode='${MODE}'; p.modulation='${MODULATION}'; p.use_pilots=true; p.used_bins=${USED_BINS}; p.pilot_bins=${PILOT_BINS}; p.interference_file='${INTERFERENCE_FILE}'; p.interference_ratio_db=${INTERFERENCE_RATIO_DB}; ${OFFSET_EXPR} p.make_plots=${MAKE_PLOTS_OCT}; p.save_images=${SAVE_IMAGES_OCT}; p.out_dir='${OUT_DIR}'; p.pause_before_exit=false; ofdm_test_channel(p);"
  exit 0
fi

if [ "$RUN_KIND" = "sweep" ]; then
  octave --quiet --eval "addpath('octave'); cfg=struct(); cfg.num_chunks=${NUM_CHUNKS}; cfg.show_progress=${SHOW_PROGRESS_OCT}; cfg.interference_file='${INTERFERENCE_FILE}'; cfg.interference_ratio_db=${INTERFERENCE_RATIO_DB}; cfg.chunk_advance_samples=${CHUNK_ADVANCE_SAMPLES}; cfg.equalizer_modes=${MODES}; ofdm_recorded_interference_sweep(cfg);"
  exit 0
fi

if [ "$RUN_KIND" = "compare" ]; then
  octave --quiet --eval "addpath('octave'); cfg=struct(); cfg.num_segments=${NUM_SEGMENTS}; cfg.show_progress=${SHOW_PROGRESS_OCT}; cfg.interference_file='${INTERFERENCE_FILE}'; cfg.interference_ratio_db=${INTERFERENCE_RATIO_DB}; cfg.start_offset_samples=${START_OFFSET_SAMPLES}; cfg.chunk_advance_samples=${CHUNK_ADVANCE_SAMPLES}; cfg.equalizer_modes=${MODES}; cfg.make_plots=${MAKE_PLOTS_OCT}; cfg.save_images=${SAVE_IMAGES_OCT}; cfg.out_dir='${OUT_DIR}'; ofdm_recorded_interference_compare(cfg);"
  exit 0
fi

printf 'unsupported RUN_KIND: %s\n' "$RUN_KIND" >&2
exit 1
