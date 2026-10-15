# Historical value-bundle experiment (four counterfactual worlds)

## Installation

```sh
python -m pip install requests numpy pandas matplotlib
export OPENROUTER_API_KEY='...'
```

**Talkie:** Use `niklassheth/talkie-1930-13b-it-GGUF` **Q5_K_M** (the *instruction* model, not `-base`). A compatible `llama-server` executable is needed because Talkie has its own `talkie` GGUF architecture; see the model publisher's llama.cpp support instructions. Example with the compatible binary installed:

```sh
export GGML_METAL_NE11_MM_MIN=1024   # Metal fp16 workaround for Talkie
papers/on-hanson-values/run.sh       # defaults: TALKIE_NGL=0, TALKIE_CTX=2048
```

On Apple Silicon, **partial Metal offload often returns garbage tokens**; `run.sh` defaults to **CPU (`TALKIE_NGL=0`)**. Full `-ngl 99` can OOM on 24 GB unified memory. Keep context at **2048** (Talkie train limit). Run interviews in the **same session** as `llama-server` (or use `run_talkie_core_both.sh`). Check `curl http://127.0.0.1:8080/v1/models` before proceeding.

**Both-arm core run:** `papers/on-hanson-values/core_both/` — 128 interviews (64 Talkie + 64 Qwen), judged and analyzed.

**Model arms:** local `talkie` via OpenAI-compatible `/v1/chat/completions` and `qwen1930` via OpenRouter `qwen/qwen3-14b`; both use the same historically situated instruction and wording. Judge model defaults to `openai/gpt-4.1-mini` via OpenRouter. No contemporary-model-only reference is included in the four-world causal contrasts.

## Design

The **four pre-specified core conditions** are:

| ID | Society | Information / exposure |
|---|---|---|
| `A_1930` | historical material setting | familiar historical society |
| `B_transplant` | changed material setting | learned facts only; no exposure to modern people or norms |
| `C_elites` | same changed setting | admired modern elites shown adopting B, without explanations |
| `D_full` | same changed setting | elite practices, peers, media praise, and institutional arrangements |

The **five additional contrast conditions** (extended default) are:

- `C_explained`: as C, but elite reasons are given explicitly; Hanson attribution contrast `C_elites - C_explained`.
- `C_opposite`: as C, but elite behavior favors A; controls for prestige direction.
- `C_moral`: as C, but elites explicitly endorse B as virtuous; isolates overt persuasion.
- `D_peers`: C + peer behavior supporting B.
- `D_media`: `D_peers` + favorable cultural narratives; `D_full - D_media` is institutional exposure.

Eight fixed toy scenarios: marriage, parenting, career loyalty, dual-earner households, kinship/mobility, honor/retaliation, smoking in shared rooms, savings/sharing. These descriptions are **synthetic historical analogies**, not measured historical statistics.

For each scenario/condition/model/replicate, the same three natural questions are asked sequentially. The first asks for advice, the second about character judgments, the third for preference when payoffs and social approval are equal. The first answer is recorded before later questions might elicit moral language. A/B order is counterbalanced across replicates, and constant for each case/replicate across models and conditions.

The independent judge sees ONLY A/B practice definitions, presented order, and the three answers; it is blind to condition, model identity, and ecology description. Strictly verified quotations are kept for auditing. Its output is judgment data, not historical ground truth.

## Run

```sh
python historical_value_experiment.py selftest
python historical_value_experiment.py plan --design core --replicates 2
python historical_value_experiment.py run --design core --replicates 2 --out pilot --limit 8
python historical_value_experiment.py run --design core --replicates 2 --out pilot
python historical_value_experiment.py judge --out pilot
python historical_value_experiment.py analyze --out pilot
```

For the extended experiment (nine conditions, two models, three replications), use a **new** output directory:

```sh
python historical_value_experiment.py plan --design extended --replicates 3
python historical_value_experiment.py run --design extended --replicates 3 --out extended
python historical_value_experiment.py judge --out extended
python historical_value_experiment.py analyze --out extended
```

Default expanded experiment: `8 × 9 × 3 × 2 = 432 interviews`, up to **1296 respondent completions** plus **432 judging calls**, excluding retries. The four-condition run with two replicates is 128 interviews.

**Batteries A+B+C (24 domains), extended, 3 replicates, both arms:** `24 × 9 × 3 × 2 = 1296` interviews. Use `--batteries abc` (default in v3). Battery **A** = original eight (new facts favor B); **B** = eight more eco→B domains; **C** = eight eco→**A** drift-test domains. Long run: `papers/on-hanson-values/run_abc_extended.sh` (Qwen phase first, then Talkie; log: `/tmp/abc-extended-r3.log`).

**Long runs on a Mac:** Jobs started from a Cursor agent shell often die when the agent turn ends (no traceback; `interviews.jsonl` stops updating). Prefer **`papers/on-hanson-values/install_abc_launchagent.sh load`** (LaunchAgent + `caffeinate`), or run in **Terminal.app** / `tmux`. Progress heartbeat: `papers/on-hanson-values/abc_extended_r3/run_state.json`. Do not `pkill` the experiment if one is already running. To reduce inference cost use `--turns 1`, which retains spontaneous behavior and moralization but not elicited moral or equalized-preference measures.

A run is resumable. Successful interview trial IDs and successful judgments are skipped. Retries after partial interview failures re-run that interview from Turn 1, not just the failed turn. Settings are frozen in `manifest.json`; using changed settings requires a new output directory. API credentials are not written to the data logs.

## Outputs

- `interviews.jsonl`: raw, dated full dialogues, including errors and model IDs.
- `judgments.jsonl`: separate blind-coded data with original judge response.
- `coded_observations.csv`: numeric coding and raw Turn 1 answer.
- `case_effects.csv`: paired differences averaged within domain and model.
- `effects.csv`: means + percentile 95% bootstrap intervals, resampling **eight domains** (not responses).
- `model_delta_comparison.csv`: Talkie vs Qwen *differences in within-model treatment effects*.
- `residual_ratios.csv`: conditional fraction `(D-B)/(D-A)` computed on population-mean differences; suppressed for |D-A| < .10.
- `effects_*.png`: separate plots for each outcome.

Coding: A=0, B=1, equal/both/mixed=.5, not stated=missing for preference-direction measures. For `spontaneous_B_moralization`, no expressed moral support of B counts as 0 (including silence); a separate `spontaneous_any_moral` output measures whether any moral position was expressed. This rough ordinal coding is deliberately a *sign/direction* metric; it should not be interpreted as a cardinal psychological strength scale. Moralization's missingness rate is itself an outcome; in particular do NOT recode `not_stated` to neutral `.5`.

## Analysis caveats

1. LLM responses are **not** equivalent to real 1930 or 2026 people. A genuine historical-validation phase needs source-dated human judgments and matched ecological evidence. Talkie has instruction tuning that may contaminate historical priors.
2. D condition *asserts* norms favorable to B; it is a causal cue-sensitivity experiment, not independently derived actual modern cultural exposure.
3. `D_full - C_elites` measures a bundle of social exposure effects; the nested arms attempt to distinguish peers, media, and institutions, but their incremental effects depend on introduction order.
4. The C opacity-vs-explanation contrast changes informational persuasion as well as attribution: signs do not uniquely identify the cognitive mechanism.
5. Pairing is within model, domain, and replicate. The case-bootstrap CI with just eight domains is exploratory; repetitions improve precision within case but do not add independent cultural domains.
6. Model outputs need blinded human cross-coding of a sample. Code should be read as hypotheses, not statistical ground truth.
7. Qwen3 may output reasoning tokens depending on its inference routing; verify visible text in the pilot and configure reasoning separately if required.
