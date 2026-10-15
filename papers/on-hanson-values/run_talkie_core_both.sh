#!/usr/bin/env bash
# Run Talkie interviews for core_both while keeping llama-server alive in this session.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$ROOT"
export GGML_METAL_NE11_MM_MIN=1024
export TALKIE_NGL=0
LOG="${TALKIE_RUN_LOG:-/tmp/talkie-core_both-run.log}"
exec >>"$LOG" 2>&1
echo "=== $(date -u +%Y-%m-%dT%H:%M:%SZ) starting Talkie core_both run ==="

papers/on-hanson-values/run.sh &
SPID=$!
cleanup() { kill "$SPID" 2>/dev/null || true; }
trap cleanup EXIT

for _ in $(seq 1 180); do
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

set -a
# shellcheck disable=SC1091
source .env 2>/dev/null || true
set +a

.venv/bin/python3.14 papers/on-hanson-values/historical_value_experiment.py run \
  --design core --replicates 2 --arms talkie,qwen1930 --out papers/on-hanson-values/core_both \
  --pause 0.05

.venv/bin/python3.14 papers/on-hanson-values/historical_value_experiment.py judge \
  --out papers/on-hanson-values/core_both

.venv/bin/python3.14 papers/on-hanson-values/historical_value_experiment.py analyze \
  --out papers/on-hanson-values/core_both

echo "=== $(date -u +%Y-%m-%dT%H:%M:%SZ) finished ==="
