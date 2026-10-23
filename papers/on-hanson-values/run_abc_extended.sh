#!/usr/bin/env bash
# Batteries A+B+C (or $BATTERIES, e.g. d; for e set DESIGN=five), extended design, 3 replicates (1296 interviews both arms for abc).
# Qwen runs first (API only); Talkie second (local llama-server). Resumable via trial_id.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$ROOT"
OUT="${1:-papers/on-hanson-values/abc_extended_r3}"
BATTERIES="${BATTERIES:-abc}"
echo "=== $(date -u +%Y-%m-%dT%H:%M:%SZ) $BATTERIES ${DESIGN:-extended} r=3 out=$OUT ==="

set -a
# shellcheck disable=SC1091
source .env 2>/dev/null || true
set +a

PY=.venv/bin/python3.14
[[ -x "$PY" ]] || PY=python3.14
CAFF=()
if command -v caffeinate >/dev/null; then
  CAFF=(caffeinate -i)
fi
RUN=( "${CAFF[@]}" "$PY" papers/on-hanson-values/historical_value_experiment.py run \
  --design "${DESIGN:-extended}" --replicates 3 --batteries "$BATTERIES" --out "$OUT" --pause 0.05 )

echo "--- phase 1: qwen1930 (OpenRouter) ---"
"${RUN[@]}" --arms qwen1930

echo "--- phase 2: talkie (local llama-server) ---"
export GGML_METAL_NE11_MM_MIN=1024
export TALKIE_NGL=0
papers/on-hanson-values/run.sh &
SPID=$!
cleanup() { kill "$SPID" 2>/dev/null || true; }
trap cleanup EXIT

for _ in $(seq 1 240); do
  if curl -sS -m 3 http://127.0.0.1:8080/v1/models 2>/dev/null | rg -q talkie; then
    echo "llama-server ready pid=$SPID"
    break
  fi
  if ! kill -0 "$SPID" 2>/dev/null; then
    echo "llama-server exited during startup"
    exit 1
  fi
  sleep 2
done

"${RUN[@]}" --arms talkie

echo "--- judge + analyze ---"
"${CAFF[@]}" "$PY" papers/on-hanson-values/historical_value_experiment.py judge --out "$OUT"
"${CAFF[@]}" "$PY" papers/on-hanson-values/historical_value_experiment.py analyze --out "$OUT"

echo "=== $(date -u +%Y-%m-%dT%H:%M:%SZ) finished ==="
