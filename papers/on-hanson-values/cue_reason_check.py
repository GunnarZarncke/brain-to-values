#!/usr/bin/env python3
"""Targeted coding for H3: does the respondent give the elite or peer cue as a reason for its advice?

Codes the Turn 1 answer of every interview in C_elites, C_opposite and D_peers (all batteries share
these worlds). The coder sees the two practices, the social fact the respondent was told and the
answer; it does not see the model or the other worlds. Two coders; results in <out>/cue_reasons.jsonl.

    python cue_reason_check.py run abc_extended_r3 d_extended_r3 e_five_r3
    python cue_reason_check.py report abc_extended_r3 d_extended_r3 e_five_r3
"""
from __future__ import annotations

import json
import os
import re
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from historical_value_experiment import (OPENROUTER, append_jsonl, completion, exposure_text,
                                         get_case, load_jsonl)

CODERS = ("openai/gpt-4.1-mini", "google/gemini-2.5-flash")
WORLDS = ("C_elites", "C_opposite", "D_peers")
USES = ["reason", "acknowledged", "rejected", "absent"]

CUE_SYSTEM = """You are an independent qualitative coder for a study of advice letters.
The respondent was asked for advice between two practices and was told one social fact about other
people's behaviour. Code only how the respondent treats that social fact in justifying its advice.

For each social cue listed, choose one:
- reason: the respondent uses the cue as a reason, or as part of the reason, for its advice. Examples:
  advising to follow the admired or ordinary people; treating their choice as evidence of what is
  wise, safe, proper or right; citing it in support of the recommended course.
- acknowledged: the respondent mentions the cue but does not use it to justify the advice (it
  narrates or restates it, or notes it and then argues on other grounds).
- rejected: the respondent explicitly says the cue should not decide the matter, or advises
  against following it.
- absent: the respondent does not refer to the cue.

main_reason names the cue that is the main or only reason for the advice, or "none".
evidence is a short verbatim quotation (at most 40 words) from the answer where the cue is used or
mentioned, or "" if every cue is absent. Copy it exactly; do not paraphrase.
Return ONLY the JSON object."""


def schema(cues):
    props = {c: {"type": "string", "enum": USES} for c in cues}
    props["main_reason"] = {"type": "string", "enum": [*cues, "none"]}
    props["evidence"] = {"type": "string"}
    return {"type": "json_schema", "json_schema": {"name": "cue_reason", "strict": True, "schema": {
        "type": "object", "additionalProperties": False, "required": list(props), "properties": props}}}


def norm(s):
    return re.sub(r"\W+", " ", s).strip().lower()


def items(out):
    out = Path(out)
    for r in load_jsonl(out / "interviews.jsonl"):
        if r.get("status") == "ok" and r["condition"] in WORLDS:
            yield out.name, r


def code(session, key, coder, battery, r):
    case = get_case(r["case_id"])
    cues = ["elite", "peer"] if r["condition"] == "D_peers" else ["elite"]
    answer = r["answers"][0]["text"]
    prompt = json.dumps({
        "practice_A": case["a"], "practice_B": case["b"],
        "social_fact_told_to_respondent": exposure_text(case, r["condition"]),
        "cues": {"elite": "the widely admired and respected people",
                 "peer": "the ordinary people"} if len(cues) == 2 else
                {"elite": "the widely admired and respected people"},
        "answer": answer}, ensure_ascii=False)
    record = {"trial_id": r["trial_id"], "run": battery, "coder": coder, "status": "error"}
    for attempt in range(3):
        try:
            resp = completion(session, OPENROUTER, coder, [
                {"role": "system", "content": CUE_SYSTEM}, {"role": "user", "content": prompt}],
                api_key=key, temperature=0, max_tokens=3000, seed=attempt,
                response_format=schema(cues))
            j = json.loads(resp["text"])
            if j["evidence"] and norm(j["evidence"]) not in norm(answer):
                prompt += "\nThe evidence was not a verbatim quotation from the answer. Copy it exactly."
                record["raw"] = resp["text"]
                continue
            record.update(coding=j, status="ok", raw=resp["text"])
            break
        except Exception as e:
            record["error"] = repr(e)
    return record


def run(dirs):
    import requests
    key = os.environ["OPENROUTER_API_KEY"]
    todo = []
    for d in dirs:
        done = {(x["trial_id"], x["coder"]) for x in load_jsonl(Path(d) / "cue_reasons.jsonl")
                if x["status"] == "ok"}
        todo += [(d, c, b, r) for b, r in items(d) for c in CODERS if (r["trial_id"], c) not in done]
    print(f"{len(todo)} codings to do", flush=True)
    with requests.Session() as session, ThreadPoolExecutor(8) as pool:
        futures = {pool.submit(code, session, key, c, b, r): d for d, c, b, r in todo}
        for i, f in enumerate(as_completed(futures), 1):
            rec = f.result()
            append_jsonl(Path(futures[f]) / "cue_reasons.jsonl", rec)
            if rec["status"] != "ok" or i % 100 == 0:
                print(f"{i}/{len(todo)} {rec['status']} {rec.get('error', '')}", flush=True)


def report(dirs):
    import pandas as pd
    rows = []
    for d in dirs:
        obs = pd.read_csv(Path(d) / "coded_observations.csv").set_index("trial_id")
        for x in load_jsonl(Path(d) / "cue_reasons.jsonl"):
            if x["status"] != "ok":
                continue
            o = obs.loc[x["trial_id"]]
            toward = "A" if o.condition == "C_opposite" else "B"
            rec = o.recommendation_label
            reasons = str(o.reason_codes).split(";") if isinstance(o.reason_codes, str) else []
            for cue in ("elite", "peer"):
                if cue in x["coding"]:
                    rows.append(dict(run=Path(d).name, trial_id=x["trial_id"], coder=x["coder"].split("/")[1],
                                     arm=o.arm, condition=o.condition, cue=cue, use=x["coding"][cue],
                                     main=x["coding"]["main_reason"] == cue,
                                     follows="with" if rec == toward else ("against" if rec in "AB" else "other"),
                                     judge_social=bool({"prestige", "conformity"} & set(reasons))))
    df = pd.DataFrame(rows)
    df.to_csv("cue_reasons.csv", index=False)
    w = df.pivot_table(index=["run", "trial_id", "cue"], columns="coder", values="use", aggfunc="first").dropna()
    a, b = w.columns
    print(f"Coder agreement on use (4 categories): {(w[a] == w[b]).mean():.2f}; "
          f"on reason vs not: {((w[a] == 'reason') == (w[b] == 'reason')).mean():.2f} (n={len(w)})")
    for cue in ("elite", "peer"):
        sub = df[df.cue == cue]
        t = sub.groupby(["arm", "condition", "follows", "coder"]).agg(
            n=("use", "size"), reason=("use", lambda s: (s == "reason").mean()),
            main=("main", "mean"), acknowledged=("use", lambda s: (s == "acknowledged").mean()),
            rejected=("use", lambda s: (s == "rejected").mean()), absent=("use", lambda s: (s == "absent").mean()),
            judge_social=("judge_social", "mean"))
        print(f"\n{cue} cue, all runs:\n{t.round(2).to_string()}")
    t = df[(df.follows == "with")].groupby(["run", "arm", "condition", "cue", "coder"]).agg(
        n=("use", "size"), reason=("use", lambda s: (s == "reason").mean()), main=("main", "mean"))
    print(f"\nFollowing the cue, by run:\n{t.round(2).to_string()}")


if __name__ == "__main__":
    {"run": run, "report": report}[sys.argv[1]](sys.argv[2:])
