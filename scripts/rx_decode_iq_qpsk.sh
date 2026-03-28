#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

OUT_DIR="$ROOT/output"
DURATION="4"
WAV="$OUT_DIR/rx_capture_iq_qpsk.wav"
LOG="acoustic_ofdm_mic_roundtrip_iq_qpsk.log"

mkdir -p "$OUT_DIR"
rm -f \
  "$OUT_DIR/tx_roundtrip_iq_qpsk.wav" \
  "$OUT_DIR/tx_packet_iq_qpsk.wav" \
  "$OUT_DIR/tx_symbols_only_iq_qpsk.wav" \
  "$OUT_DIR/rx_capture_iq_qpsk.wav" \
  "$OUT_DIR/rx_spectrogram_iq_qpsk.png" \
  "$OUT_DIR/ofdm_decode_bins.csv" \
  "$OUT_DIR/ofdm_constellation.png" \
  "$OUT_DIR/ofdm_constellation_pre_eq.csv" \
  "$OUT_DIR/ofdm_constellation_post_eq.csv" \
  "$OUT_DIR/ofdm_channel_compare.csv" \
  "$OUT_DIR/ofdm_channel_compare.png" \
  "$OUT_DIR/ofdm_pre_crc_bytes.bin" \
  "$ROOT/$LOG"

cargo run -p acoustic_ofdm_cli -- \
  encode \
  --fec-mode hamming74 \
  --passband-mode iq \
  --fs-baseband 16000 \
  --nfft 2048 \
  --ncp 1024 \
  --sync-half-len 2048 \
  --modulation qpsk \
  "$OUT_DIR/tx_packet_iq_qpsk.wav" \
  ACOUSTIC-OFDM-ORACLE

cargo run -p acoustic_ofdm_cli -- \
  encode-body \
  --fec-mode hamming74 \
  --passband-mode iq \
  --fs-baseband 16000 \
  --nfft 2048 \
  --ncp 1024 \
  --sync-half-len 2048 \
  --modulation qpsk \
  "$OUT_DIR/tx_symbols_only_iq_qpsk.wav" \
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
  --modulation qpsk \
  --duration-sec "$DURATION" \
  --mic-gain 0.6 \
  --spk-gain 0.6 \
  --repeats 1 \
  --gap-sec 0.15 \
  --pre-delay-sec 0.05 \
  --oracle \
  --spectrogram \
  --dump-tx-wav "$OUT_DIR/tx_roundtrip_iq_qpsk.wav" \
  --dump-wav "$WAV" \
  --spectrogram-path "$OUT_DIR/rx_spectrogram_iq_qpsk.png" \
  --log-level debug \
  --log-file "$LOG"
