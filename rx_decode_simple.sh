#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT"

OUT_DIR="$ROOT/output"
DURATION="4"
WAV="$OUT_DIR/rx_capture.wav"
LOG="acoustic_ofdm_mic_roundtrip.log"

mkdir -p "$OUT_DIR"
rm -f \
  "$OUT_DIR/tx_roundtrip.wav" \
  "$OUT_DIR/tx_packet.wav" \
  "$OUT_DIR/tx_symbols_only.wav" \
  "$OUT_DIR/rx_capture.wav" \
  "$OUT_DIR/rx_spectrogram.png" \
  "$OUT_DIR/rx_capture_check.png" \
  "$OUT_DIR/ofdm_decode_bins.csv" \
  "$OUT_DIR/ofdm_constellation.png" \
  "$OUT_DIR/ofdm_constellation_pre_eq.csv" \
  "$OUT_DIR/ofdm_constellation_post_eq.csv" \
  "$OUT_DIR/ofdm_bins.csv" \
  "$OUT_DIR/ofdm_sync_metric.csv" \
  "$ROOT/$LOG"

cargo run -p acoustic_ofdm_cli -- \
  encode \
  "$OUT_DIR/tx_packet.wav" \
  ACOUSTIC-OFDM-ORACLE

cargo run -p acoustic_ofdm_cli -- \
  encode-body \
  "$OUT_DIR/tx_symbols_only.wav" \
  ACOUSTIC-OFDM-ORACLE

exec cargo run -p acoustic_ofdm_cli -- \
  mic-roundtrip \
  --profile live-debug \
  --duration-sec "$DURATION" \
  --mic-gain 0.6 \
  --spk-gain 0.6 \
  --repeats 1 \
  --gap-sec 0.15 \
  --pre-delay-sec 0.05 \
  --oracle \
  --spectrogram \
  --dump-tx-wav "$OUT_DIR/tx_roundtrip.wav" \
  --dump-wav "$WAV" \
  --spectrogram-path "$OUT_DIR/rx_spectrogram.png" \
  --log-level debug \
  --log-file "$LOG"
