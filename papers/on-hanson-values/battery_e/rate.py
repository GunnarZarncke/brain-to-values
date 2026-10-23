"""Blind moral/practical ratings for battery E candidates and the 32 domains of batteries A-D (see prereg.txt)."""
import csv
import json
import os
import random
import re
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import requests

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from historical_value_experiment import CASES  # noqa: E402

RATERS = ["google/gemini-2.5-flash", "meta-llama/llama-3.3-70b-instruct", "mistralai/mistral-large"]
PROMPT = """Consider an ordinary American adult in 1930 choosing between two practices.
Topic: {topic}
Practice 1: {p1}
Practice 2: {p2}
For such a person, is the choice between these practices mainly a question of right and wrong (morality, decency, duty, religion, obligations to others) or mainly a question of what works best for the person (cost, convenience, health, efficiency, skill)?
Answer on a scale from 1 to 7: 1 = entirely a question of what works best, 4 = equally both, 7 = entirely a question of right and wrong.
Reply with JSON only: {{"rating": <integer 1-7>, "reason": "<one sentence>"}}"""


def items():
    for r in (row for f in ["candidates.csv", "candidates_round2.csv"] for row in csv.DictReader(open(HERE / f))):
        yield {"set": "candidate", "id": r["id"], "cell": r["cell"], "topic": r["topic"], "a": r["a"], "b": r["b"]}
    for c in CASES:
        yield {"set": "abcd", "id": c["id"], "cell": "", "topic": c["topic"], "a": c["a"], "b": c["b"]}


def rate(item, model, key):
    swap = random.Random(f"{item['set']}:{item['id']}").random() < 0.5
    p1, p2 = (item["b"], item["a"]) if swap else (item["a"], item["b"])
    prompt = PROMPT.format(topic=item["topic"], p1=p1, p2=p2)
    for attempt in range(4):
        for wait in [0, 5, 15, 30, 60, 120, 240]:
            # HTTP errors (e.g. upstream 429) are retried with backoff; they are not rating attempts.
            time.sleep(wait)
            r = requests.post("https://openrouter.ai/api/v1/chat/completions", timeout=120,
                              headers={"Authorization": f"Bearer {key}"},
                              json={"model": model, "temperature": 0, "max_tokens": 1500,
                                    "messages": [{"role": "user", "content": prompt}]})
            if r.ok and "choices" in r.json():
                break
        text = r.json().get("choices", [{}])[0].get("message", {}).get("content") or ""
        m = re.search(r"\{.*\}", text, re.S)
        try:
            out = json.loads(m.group(0))
            rating = int(out["rating"])
            if 1 <= rating <= 7:
                return {**item, "rater": model, "swap": swap, "rating": rating,
                        "reason": out.get("reason", ""), "attempts": attempt + 1}
        except Exception:
            pass
    return {**item, "rater": model, "swap": swap, "rating": None, "reason": text[:300], "attempts": 4}


def main():
    key = os.environ["OPENROUTER_API_KEY"]
    path = HERE / "ratings.jsonl"
    done = {}
    if path.exists():
        for r in map(json.loads, open(path)):
            if r["rating"] is not None:
                done[(r["set"], r["id"], r["rater"])] = r
    jobs = [(i, m) for i in items() for m in RATERS if (i["set"], i["id"], m) not in done]
    with ThreadPoolExecutor(int(os.environ.get("RATE_WORKERS", 8))) as ex:
        new = list(ex.map(lambda j: rate(j[0], j[1], key), jobs))
    rows = list(done.values()) + new
    with open(path, "w") as f:
        for r in rows:
            f.write(json.dumps(r) + "\n")
    print(f"{len(rows)} ratings, {sum(r['rating'] is None for r in rows)} failed")


if __name__ == "__main__":
    main()
