#!/usr/bin/env bash
# Battery E (battery_e/prereg.txt): 40 domains x 5 worlds x 3 replicates x 2 arms = 1200 interviews.
# Qwen (OpenRouter) and Talkie (local llama-server) run in parallel into one output directory, in blocked
# order (complete domain/replicate blocks, cells alternating). Every $INTERIM seconds the new interviews are
# judged and the analysis is refreshed, so partial results appear in $OUT/battery_e_contrasts.csv.
# Interim results are for monitoring only; the run is completed regardless. Resumable via trial_id.
# Deviation 4 (the four extra worlds): CONDITIONS=C_explained,C_moral,D_media,D_full MERGE_WITH=<five-world dir>
# MERGED=<dir> bash run_battery_e.sh <new dir>; the analysis then runs on the merged nine-world directory.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$ROOT"
OUT="${1:-papers/on-hanson-values/e_five_r3}"
INTERIM="${INTERIM:-1800}"
echo "=== $(date -u +%Y-%m-%dT%H:%M:%SZ) battery e five r=3 blocked out=$OUT ==="

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
EXP=papers/on-hanson-values/historical_value_experiment.py
RUN=( "$PY" "$EXP" run --design five --order blocked --replicates 3 --batteries e --out "$OUT" --pause 0.05
      ${CONDITIONS:+--conditions "$CONDITIONS"} )
AOUT="$OUT"
if [[ -n "${MERGE_WITH:-}" ]]; then
  AOUT="${MERGED:?set MERGED to the merged output directory}"
fi
analyze() {
  if [[ -n "${MERGE_WITH:-}" ]]; then
    "$PY" "$EXP" merge --sources "$MERGE_WITH" "$OUT" --out "$AOUT"
  fi
  MPLBACKEND=Agg "$PY" "$EXP" analyze --out "$AOUT" --bootstrap 2000
}
mkdir -p "$OUT"

export GGML_METAL_NE11_MM_MIN=1024
export TALKIE_NGL=0
papers/on-hanson-values/run.sh > "$OUT/llama-server.log" 2>&1 &
SPID=$!
PIDS=()
cleanup() { kill "$SPID" ${PIDS[@]+"${PIDS[@]}"} 2>/dev/null || true; }
trap cleanup EXIT

for _ in $(seq 1 240); do
  if curl -sS -m 3 http://127.0.0.1:8080/v1/models 2>/dev/null | grep -q talkie; then
    echo "llama-server ready pid=$SPID"
    break
  fi
  if ! kill -0 "$SPID" 2>/dev/null; then
    echo "llama-server exited during startup"
    exit 1
  fi
  sleep 2
done

"${CAFF[@]}" "${RUN[@]}" --arms qwen1930 > "$OUT/run_qwen1930.log" 2>&1 &
PIDS+=($!)
sleep 5  # the first arm writes the manifest
"${CAFF[@]}" "${RUN[@]}" --arms talkie > "$OUT/run_talkie.log" 2>&1 &
PIDS+=($!)
echo "arms running: qwen1930 pid=${PIDS[0]} talkie pid=${PIDS[1]}"

interim() {
  "$PY" "$EXP" judge --out "$OUT" > "$OUT/judge_interim.log" 2>&1 || true
  analyze > "$OUT/analyze_interim.log" 2>&1 || true
  echo "--- interim $(date -u +%H:%M:%SZ): $(grep -c '"status": "ok"' "$OUT/interviews.jsonl" 2>/dev/null || echo 0) interviews ok"
  sed -n '/Battery E complete/,/^Planned interviews/p' "$OUT/analyze_interim.log" || true
}

while kill -0 "${PIDS[0]}" 2>/dev/null || kill -0 "${PIDS[1]}" 2>/dev/null; do
  for _ in $(seq 1 $((INTERIM / 10))); do
    sleep 10
    kill -0 "${PIDS[0]}" 2>/dev/null || kill -0 "${PIDS[1]}" 2>/dev/null || break
  done
  interim
done
wait "${PIDS[@]}" || echo "an arm exited with an error; see $OUT/run_*.log"

echo "--- final judge + analyze ---"
"${CAFF[@]}" "$PY" "$EXP" judge --out "$OUT"
analyze
echo "=== $(date -u +%Y-%m-%dT%H:%M:%SZ) finished ==="
