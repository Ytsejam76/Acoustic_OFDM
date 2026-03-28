#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT"

OUT_DIR="$ROOT/output"
mkdir -p "$OUT_DIR"

exec cargo run -p acoustic_ofdm_cli -- \
  encode-body \
  "$OUT_DIR/tx_symbols_only.wav" \
  ACOUSTIC-OFDM-ORACLE
