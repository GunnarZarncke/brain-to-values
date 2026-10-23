"""Apply the frozen classification rule (prereg.txt) to ratings.jsonl and select domains."""
import csv
import json
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
RATERS = ["google/gemini-2.5-flash", "meta-llama/llama-3.3-70b-instruct", "mistralai/mistral-large"]
# Post-hoc split of batteries A-D used in the discussion that motivated battery E.
POSTHOC_MORAL = {"honor", "charity", "alcohol", "speech", "defer", "parenting", "smoking",
                 "religion", "premarital", "samesex", "intermarriage", "divorce", "spanking", "war", "medicine"}


def classify(ratings):
    if len(ratings) < len(RATERS) or None in ratings:
        return "failed"
    score = sum(ratings) / len(ratings)
    if score >= 5.5 and min(ratings) >= 5:
        return "moral"
    if score <= 2.5 and max(ratings) <= 3:
        return "practical"
    return "ambiguous"


def main():
    by = defaultdict(dict)
    for r in map(json.loads, open(HERE / "ratings.jsonl")):
        by[(r["set"], r["id"])][r["rater"]] = r["rating"]
    sources = {}
    if (HERE / "sources.csv").exists():
        sources = {r["id"]: r["verdict"] for r in csv.DictReader(open(HERE / "sources.csv"))}

    out = []
    for c in csv.DictReader(open(HERE / "candidates.csv")):
        rs = [by[("candidate", c["id"])].get(m) for m in RATERS]
        cls = classify(rs)
        want = c["cell"].split("_")[0]
        out.append({**c, "ratings": rs, "class": cls, "rule_pass": cls == want,
                    "source": sources.get(c["id"], "pending")})
    selected = defaultdict(list)
    for r in sorted(out, key=lambda r: (r["cell"], int(r["priority"]))):
        if r["rule_pass"] and r["source"] == "PASS" and len(selected[r["cell"]]) < 10:
            selected[r["cell"]].append(r["id"])
            r["selected"] = True
    with open(HERE / "screening.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["id", "cell", "priority", *[m.split("/")[1] for m in RATERS], "mean", "class",
                    "rule_pass", "source", "selected"])
        for r in out:
            rs = r["ratings"]
            mean = sum(rs) / len(rs) if None not in rs else ""
            w.writerow([r["id"], r["cell"], r["priority"], *rs, f"{mean:.2f}" if mean != "" else "",
                        r["class"], r["rule_pass"], r["source"], r.get("selected", False)])
    for cell in sorted({r["cell"] for r in out}):
        rows = [r for r in out if r["cell"] == cell]
        print(f"\n{cell}: rule pass {sum(r['rule_pass'] for r in rows)}/{len(rows)}; selected {selected[cell]}")
        for r in rows:
            print(f"  {r['priority']:>2} {r['id']:13s} {str(r['ratings']):14s} {r['class']:9s} source={r['source']}")

    print("\nRule check on batteries A-D (post-hoc split vs frozen rule):")
    agree = defaultdict(int)
    for (s, i), d in sorted(by.items()):
        if s != "abcd":
            continue
        rs = [d.get(m) for m in RATERS]
        cls = classify(rs)
        post = "moral" if i in POSTHOC_MORAL else "practical"
        agree[(post, cls)] += 1
        print(f"  {i:13s} {str(rs):14s} rule={cls:9s} posthoc={post}")
    print("  counts (posthoc, rule):", dict(agree))



def select_round2():
    """Deviation 1: ordered pools in pool_round2.csv, sources from both rounds (round 2 overrides)."""
    by = defaultdict(dict)
    for r in map(json.loads, open(HERE / "ratings.jsonl")):
        if r["set"] == "candidate":
            by[r["id"]][r["rater"]] = r["rating"]
    cells = {r["id"]: r["cell"] for f in ["candidates.csv", "candidates_round2.csv"]
             for r in csv.DictReader(open(HERE / f))}
    sources = {}
    for f in ["sources.csv", "sources_round2.csv"]:
        sources.update({r["id"]: r for r in csv.DictReader(open(HERE / f))})
    rows, selected = [], defaultdict(list)
    for p in csv.DictReader(open(HERE / "pool_round2.csv")):
        i = p["id"]
        rs = [by[i].get(m) for m in RATERS]
        cls = classify(rs)
        src = sources.get(i, {}).get("verdict", "pending")
        ok = cls == p["cell"].split("_")[0] and src.startswith("PASS") and cells[i] == p["cell"]
        sel = ok and len(selected[p["cell"]]) < 10
        if sel:
            selected[p["cell"]].append(i)
        rows.append([p["cell"], p["order"], i, *rs, cls, src, sources.get(i, {}).get("known_failed_by_1930", ""), sel])
    with open(HERE / "selection_round2.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["cell", "order", "id", *[m.split("/")[1] for m in RATERS], "class", "source",
                    "known_failed_by_1930", "selected"])
        w.writerows(rows)
    print("\nRound 2 selection (deviation 1):")
    for r in rows:
        print(f"  {r[0]:20s} {r[1]:>2} {r[2]:13s} {str(r[3:6]):12s} {r[6]:9s} {r[7]:16s} kf={r[8] or '-'} {'SELECTED' if r[9] else ''}")
    for cell, ids in selected.items():
        print(f"  {cell}: {len(ids)}")


if __name__ == "__main__":
    main()
    if (HERE / "pool_round2.csv").exists():
        select_round2()
