#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"

LLAMA_SERVER="${LLAMA_SERVER:-$HOME/tools/llama.cpp/build/bin/llama-server}"
MODEL="${TALKIE_MODEL:-$HOME/models/talkie/talkie-Q5_K_M.gguf}"
PORT="${TALKIE_PORT:-8080}"

if [[ ! -x "$LLAMA_SERVER" ]]; then
  echo "Missing llama-server at $LLAMA_SERVER" >&2
  exit 1
fi
if [[ ! -f "$MODEL" ]]; then
  echo "Missing Talkie GGUF at $MODEL" >&2
  exit 1
fi

# Avoid Metal fp16 overflow on Talkie batched matmul (PocketAI / upstream notes).
export GGML_METAL_NE11_MM_MIN="${GGML_METAL_NE11_MM_MIN:-1024}"

# Partial Metal offload often yields 0x1f garbage tokens on Talkie; CPU (0) is slower but stable.
# Full GPU (-ngl 99) can OOM on 24 GB unified memory while loading Q5_K_M.
NG_LAYERS="${TALKIE_NGL:-0}"
# Talkie GGUF reports n_ctx_train=2048; llama-server caps slots to that.
CTX="${TALKIE_CTX:-2048}"

exec "$LLAMA_SERVER" \
  -m "$MODEL" \
  -c "$CTX" \
  -ngl "$NG_LAYERS" \
  --parallel 1 \
  --host 127.0.0.1 \
  --port "$PORT" \
  --alias talkie1930
