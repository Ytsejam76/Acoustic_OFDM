#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

OUT_DIR="$ROOT/output"
mkdir -p "$OUT_DIR"

exec cargo run -p acoustic_ofdm_cli -- \
  tx \
  --profile live-debug \
  --fec-mode hamming74 \
  --passband-mode iq \
  --fs-baseband 22050 \
  --nfft 1024 \
  --ncp 512 \
  --sync-half-len 1024 \
  --spk-gain 0.6 \
  --oracle \
  --pre-delay-sec 0.05 \
  --repeats 8 \
  --gap-sec 0.15 \
  --dump-wav "$OUT_DIR/tx_iq.wav"
