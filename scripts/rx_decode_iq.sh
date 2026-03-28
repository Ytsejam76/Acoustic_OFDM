#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

OUT_DIR="$ROOT/output"
DURATION="4"
WAV="$OUT_DIR/rx_capture_iq.wav"
LOG="acoustic_ofdm_mic_roundtrip_iq.log"

mkdir -p "$OUT_DIR"
rm -f \
  "$OUT_DIR/tx_roundtrip_iq.wav" \
  "$OUT_DIR/tx_packet_iq.wav" \
  "$OUT_DIR/tx_symbols_only_iq.wav" \
  "$OUT_DIR/rx_capture_iq.wav" \
  "$OUT_DIR/rx_spectrogram_iq.png" \
  "$OUT_DIR/ofdm_decode_bins.csv" \
  "$OUT_DIR/ofdm_constellation.png" \
  "$OUT_DIR/ofdm_constellation_pre_eq.csv" \
  "$OUT_DIR/ofdm_constellation_post_eq.csv" \
  "$OUT_DIR/ofdm_channel_compare.csv" \
  "$OUT_DIR/ofdm_channel_compare.png" \
  "$ROOT/$LOG"

cargo run -p acoustic_ofdm_cli -- \
  encode \
  --fec-mode hamming74 \
  --passband-mode iq \
  --fs-baseband 16000 \
  --nfft 2048 \
  --ncp 1024 \
  --sync-half-len 2048 \
  "$OUT_DIR/tx_packet_iq.wav" \
  ACOUSTIC-OFDM-ORACLE

cargo run -p acoustic_ofdm_cli -- \
  encode-body \
  --fec-mode hamming74 \
  --passband-mode iq \
  --fs-baseband 16000 \
  --nfft 2048 \
  --ncp 1024 \
  --sync-half-len 2048 \
  "$OUT_DIR/tx_symbols_only_iq.wav" \
  ACOUSTIC-OFDM-ORACLE

exec cargo run -p acoustic_ofdm_cli -- \
  mic-roundtrip \
  --profile live-debug \
  --fec-mode hamming74 \
  --passband-mode iq \
  --fs-baseband 16000 \
  --nfft 2048 \
  --ncp 1024 \
  --sync-half-len 2048 \
  --duration-sec "$DURATION" \
  --mic-gain 0.6 \
  --spk-gain 0.6 \
  --repeats 1 \
  --gap-sec 0.15 \
  --pre-delay-sec 0.05 \
  --oracle \
  --spectrogram \
  --dump-tx-wav "$OUT_DIR/tx_roundtrip_iq.wav" \
  --dump-wav "$WAV" \
  --spectrogram-path "$OUT_DIR/rx_spectrogram_iq.png" \
  --log-level debug \
  --log-file "$LOG"
