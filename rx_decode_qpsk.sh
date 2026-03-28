#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT"

OUT_DIR="$ROOT/output"
DURATION="4"
WAV="$OUT_DIR/rx_capture_qpsk.wav"
LOG="acoustic_ofdm_mic_roundtrip_qpsk.log"

mkdir -p "$OUT_DIR"
rm -f \
  "$OUT_DIR/tx_roundtrip_qpsk.wav" \
  "$OUT_DIR/tx_packet_qpsk.wav" \
  "$OUT_DIR/tx_symbols_only_qpsk.wav" \
  "$OUT_DIR/rx_capture_qpsk.wav" \
  "$OUT_DIR/rx_spectrogram_qpsk.png" \
  "$OUT_DIR/ofdm_constellation.png" \
  "$OUT_DIR/ofdm_constellation_pre_eq.csv" \
  "$OUT_DIR/ofdm_constellation_post_eq.csv" \
  "$ROOT/$LOG"

cargo run -p acoustic_ofdm_cli -- \
  encode \
  --modulation qpsk \
  "$OUT_DIR/tx_packet_qpsk.wav" \
  ACOUSTIC-OFDM-ORACLE

cargo run -p acoustic_ofdm_cli -- \
  encode-body \
  --modulation qpsk \
  "$OUT_DIR/tx_symbols_only_qpsk.wav" \
  ACOUSTIC-OFDM-ORACLE

exec cargo run -p acoustic_ofdm_cli -- \
  mic-roundtrip \
  --profile live-debug \
  --modulation qpsk \
  --duration-sec "$DURATION" \
  --mic-gain 0.6 \
  --spk-gain 0.6 \
  --repeats 1 \
  --gap-sec 0.15 \
  --pre-delay-sec 0.05 \
  --oracle \
  --spectrogram \
  --dump-tx-wav "$OUT_DIR/tx_roundtrip_qpsk.wav" \
  --dump-wav "$WAV" \
  --spectrogram-path "$OUT_DIR/rx_spectrogram_qpsk.png" \
  --log-level debug \
  --log-file "$LOG"
