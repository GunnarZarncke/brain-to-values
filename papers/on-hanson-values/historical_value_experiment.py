#!/usr/bin/env python3
"""1930 counterfactual value-bundle experiment: natural interviews, blind coding, paired analysis.

Commands: plan, run, judge, analyze, selftest.
Dependencies: requests, numpy, pandas, matplotlib (analysis only).
Local Talkie must be served by a compatible OpenAI-style llama-server.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import random
import re
import sys
import time
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

VERSION = "hanson-four-worlds-natural-v3"
OPENROUTER = "https://openrouter.ai/api/v1"
DEFAULT_QWEN = "qwen/qwen3-14b"
DEFAULT_JUDGE = "openai/gpt-4.1-mini"

CORE = ("A_1930", "B_transplant", "C_elites", "D_full")
EXTENDED = (
    "A_1930", "B_transplant", "C_elites", "D_full",
    "C_explained", "C_opposite", "C_moral",
    "D_peers", "D_media",
)

# The scenarios are deliberately simplified hypotheses, NOT historical facts.
# Both models receive precisely the same case contents and treatment language.
# Battery A: new-world facts favor practice B. Battery C: new-world facts favor A (drift test).
CASES_BATTERY_A = [
    dict(id="marriage", topic="When to marry",
         a="marry while young and build a household together",
         b="complete extended education and secure earnings before marrying",
         old="Young adults can usually earn an independent living after relatively short training. A couple can begin a household while still gaining experience and possessions.",
         new="Desirable work requires many years of education. Few young adults can earn enough for secure housing until they finish their training.",
         explanation="long training and insecure early earnings made waiting more practical",
         institution="Employers and educational establishments now organize the usual path to a career around many years of training before financial independence"),
    dict(id="parenting", topic="How closely to direct children's time",
         a="let children spend much of their free time in independent play with ordinary neighborhood supervision",
         b="personally organize and supervise much of children's free time",
         old="Children commonly play in neighborhood groups; formal admissions and examinations do not depend much on organized activities during childhood.",
         new="Access to desirable education increasingly depends on competitive examinations and childhood accomplishments; supervised coaching substantially improves results.",
         explanation="organized coaching gives children an advantage in increasingly competitive examinations",
         institution="Schools and clubs give favorable consideration to documented supervised achievements"),
    dict(id="career", topic="Loyalty to an employer",
         a="remain with one employer for most of a working life",
         b="change employers when the move improves earnings or skills",
         old="Long service brings unusually reliable employment and steadily increasing pay; a worker's experience is difficult to use elsewhere.",
         new="Skills transfer readily between employers, long service brings little security, and competing firms can offer much better prospects.",
         explanation="portable skills and insecure tenure make changing employers useful",
         institution="Professional associations and employment arrangements recognize qualifications and benefits when employees change firms"),
    dict(id="household", topic="Paid work and household duties",
         a="have one married adult concentrate on paid work and the other on household care",
         b="have both married adults earn wages and divide household care",
         old="Domestic tasks require much time, childcare beyond the household is scarce, and a second wage earner usually has limited opportunities.",
         new="Reliable childcare and labor-saving household equipment are inexpensive, and both adults have good chances for paid employment.",
         explanation="improvements in childcare and work opportunities make a shared arrangement worthwhile",
         institution="Employers and childcare services provide practical support for families in which both adults work"),
    dict(id="mobility", topic="Moving away from one's relatives",
         a="remain near relatives and exchange practical help with them",
         b="move farther away for better opportunities and establish new networks",
         old="Travel and communication are expensive; insurance is uncertain and family members nearby provide most help in emergencies.",
         new="Travel and long-distance communication are inexpensive, insurance is dependable, and work in distant towns pays substantially more.",
         explanation="insurance and reliable travel make leaving relatives less hazardous",
         institution="Employment benefits, insurance, and qualifications remain valid when residents move to other towns"),
    dict(id="honor", topic="How to respond to a serious public insult or injury",
         a="personally retaliate when seriously wronged so that others know one will not be intimidated",
         b="seek impartial mediation or a court ruling rather than retaliate",
         old="Courts seldom enforce judgments and aggressors often exploit people who do not visibly defend themselves.",
         new="Courts reliably enforce judgments, prevent repeat attacks, and impose penalties on private retaliation.",
         explanation="impartial enforcement is safer and more effective than private retaliation",
         institution="The courts give rapid redress and impose substantial penalties for private retaliation"),
    dict(id="smoking", topic="Smoking in shared rooms",
         a="permit smoking in enclosed rooms during social gatherings",
         b="keep enclosed common rooms free of smoking while allowing it outdoors",
         old="There is little generally accepted evidence that tobacco smoke harms people who are not smoking themselves.",
         new="Repeated and careful investigations have demonstrated serious harm to nonsmokers from breathing tobacco smoke in enclosed rooms.",
         explanation="evidence of injury to nonsmokers gives a practical reason to keep shared rooms clear",
         institution="Meeting halls and employers prohibit indoor smoking but provide places to smoke outdoors"),
    dict(id="saving", topic="Keeping household surplus",
         a="give much of an occasional surplus to neighbors who presently need help",
         b="keep substantial private reserves against future shortages",
         old="Income and food supplies are fairly steady; emergencies rarely strike many households at once and neighbors ordinarily repay one another's help.",
         new="Income and food shortages are frequent, severe, and shared by many households at once, leaving little help available from neighbors when a crisis comes.",
         explanation="widespread shortages make reserves safer than relying on reciprocal aid",
         institution="Savings associations protect private emergency reserves and encourage regular saving"),
]
CASES_BATTERY_B = [
    dict(id="fertility", topic="Family size and timing",
         a="have several children, starting while parents are still young",
         b="have few children, spaced and planned, with parents well established first",
         old="Many hands are needed at home and on the land; infant loss is common and larger families spread risk.",
         new="Raising each child well requires long schooling and expense; smaller, later families are safer and more prosperous.",
         explanation="high per-child investment makes smaller families the prudent choice",
         institution="Schools and employers expect parents to limit family size so each child can be properly supported"),
    dict(id="charity", topic="Whom to help with spare means",
         a="give first to relatives and neighbors who are known to you",
         b="give to effective charitable causes even when recipients are distant and unknown",
         old="Needs are visible nearby; gifts are repaid by custom; distant appeals are hard to verify.",
         new="Reliable reports identify great need far away; organized charities deliver aid efficiently across distance.",
         explanation="trustworthy information makes distant giving unusually effective",
         institution="Charitable societies publish audited results and encourage gifts to distant recipients"),
    dict(id="credit", topic="Using borrowed money",
         a="avoid debt except in dire emergency and repay as quickly as possible",
         b="borrow willingly for education, housing, or tools when expected returns exceed the cost",
         old="Interest is heavy, lenders are harsh, and a bad season can ruin a household that owes money.",
         new="Credit is regulated, terms are predictable, and investment in skills or property often pays for itself.",
         explanation="stable credit markets make prudent borrowing worthwhile",
         institution="Banks and employers offer installment plans and treat responsible borrowing as normal"),
    dict(id="leisure", topic="Hours of work versus rest",
         a="work long hours whenever work is available and rest only when necessary",
         b="limit working hours and protect regular time for rest, family, and local life",
         old="Incomes are modest and many tasks cannot be postponed; idleness invites hardship.",
         new="Productivity is high enough that shorter hours still provide a comfortable living.",
         explanation="greater output per hour makes shorter schedules affordable",
         institution="Employers and unions standardize limited hours and paid time away from work"),
    dict(id="eldercare", topic="Care of aged parents",
         a="keep aged parents in the household and care for them directly",
         b="place aged parents in professional homes while visiting often",
         old="Institutions are scarce and mistrusted; family care is the honorable and practical course.",
         new="Reputable homes are affordable, relieve strain on earners, and provide skilled nursing.",
         explanation="reliable institutions make outside care safer than overburdening the household",
         institution="Insurance and pension arrangements assume professional care for the aged"),
    dict(id="literacy", topic="Children's schooling versus early work",
         a="send children to brief schooling then into paid work or apprenticeship young",
         b="keep children in formal schooling through late adolescence before full-time earning",
         old="Family income needs children's help; advanced schooling is for a small elite.",
         new="Credentials from long schooling are required for desirable and secure employment.",
         explanation="employers hire on documented schooling rather than early experience alone",
         institution="Public schools and employers expect completion of many years of formal study"),
    dict(id="alcohol", topic="Drinking at social gatherings",
         a="abstain from alcohol or drink very little at gatherings",
         b="take moderate social drinking for granted at ordinary gatherings",
         old="Drink leads to quarrels and ruin; respectable people and many churches discourage it.",
         new="Supply is safer, medical opinion treats moderate use as normal, and gatherings expect it.",
         explanation="changed medical and commercial norms make moderate drinking ordinary",
         institution="Licensed establishments and hosts assume responsible social drinking at events"),
    dict(id="speech", topic="Open disagreement in the community",
         a="avoid public quarrels and preserve harmony even when others act wrongly",
         b="speak plainly in public and criticize bad conduct when it affects the community",
         old="Small communities depend on face-to-face peace; open insults invite lasting feuds.",
         new="Large institutions and markets need open criticism to limit fraud and incompetence.",
         explanation="anonymous scale makes quiet tolerance more costly than frank speech",
         institution="Newspapers and public meetings protect plain speech about public affairs"),
]
CASES_BATTERY_C = [
    dict(id="homestead", topic="Where to make one's life",
         a="remain on family land or in a long-known town",
         b="move to a large distant city for career and novelty",
         old="Local ties provide security; cities are costly and morally suspect to many.",
         new="Remote work and local opportunity make staying viable; distant cities are expensive and unstable.",
         explanation="city connections still open careers that local work does not",
         institution="Employers and lenders treat permanent relocation to major cities as the usual successful path"),
    dict(id="multigen", topic="When adult children leave home",
         a="have adult children remain at home until they marry or inherit responsibility",
         b="encourage adult children to leave early for independence and self-support",
         old="Pooled households share costs; unmarried adults commonly stay under the parental roof.",
         new="Housing is scarce and dear; a shared household saves expense and weathers unemployment.",
         explanation="early independence trains self-reliance even when housing is dear",
         institution="Housing markets and employers expect young adults to maintain separate households"),
    dict(id="frugality", topic="Use of surplus income",
         a="save, repair, and avoid conspicuous display",
         b="spend freely on comforts, travel, and visible upgrades",
         old="Thrift is virtue; display invites envy and waste scarce goods.",
         new="Prices are unsettled and credit uncertain; buffers and maintenance reduce ruinous shocks.",
         explanation="visible spending signals creditworthiness and belonging",
         institution="Credit companies and retailers promote regular upgrading and discretionary spending"),
    dict(id="kinwork", topic="Hiring for a position of trust",
         a="prefer a reliable relative for the post",
         b="hire strictly on merit without favor to kin",
         old="Kin ties ensure loyalty; strangers may cheat when oversight is weak.",
         new="Fraud and broken contracts are common; known relatives are easier to hold accountable.",
         explanation="open hiring attracts talent that kin networks cannot",
         institution="Large firms and governments forbid nepotism and require open competition for posts"),
    dict(id="defer", topic="Speaking to elders, clergy, or an employer",
         a="defer openly and avoid contradicting them in company",
         b="assert equal right to speak regardless of rank or age",
         old="Rank has duties both ways; public contradiction shames the young and unsettles order.",
         new="Complex work fails when every order is questioned; clear chains of command prevent chaos.",
         explanation="speaking as an equal surfaces errors that rank conceals",
         institution="Schools and workplaces teach equal voice and discourage deference to rank"),
    dict(id="localmate", topic="Choice of spouse",
         a="marry within one's church or long-known community",
         b="marry for affection without regard to community background",
         old="Parents and neighbors know candidates; local marriage binds useful alliances.",
         new="Long-distance courtship is hard to vet; local families share reputations and enforce repair.",
         explanation="affection predicts a lasting household better than local alliance",
         institution="Mobility and media treat marrying outside one's circle as normal and admirable"),
    dict(id="trades", topic="Acceptable work for a respectable person",
         a="accept steady manual trade or farm work without shame",
         b="refuse manual work and hold out for genteel or office employment only",
         old="Honest labor with the hands is respected; idleness is not.",
         new="Skilled trades pay reliably while office openings are scarce for the inexperienced.",
         explanation="office work still confers standing even when trades pay more",
         institution="Schools and families steer the ambitious toward credentials and away from manual trades"),
    dict(id="prep", topic="Household readiness for trouble",
         a="keep stores, basic skills, and mutual aid with neighbors for emergencies",
         b="rely on state services and specialists rather than private preparedness",
         old="Help is local and slow to arrive; prudent households prepare together.",
         new="Supply interruptions and slow official response make private stores and skills lifesaving.",
         explanation="specialists are better placed than households to handle rare disasters",
         institution="Governments assure citizens that professional services make private stockpiling unnecessary"),
]
BATTERY_MAP = {"a": CASES_BATTERY_A, "b": CASES_BATTERY_B, "c": CASES_BATTERY_C}
CASES = CASES_BATTERY_A + CASES_BATTERY_B + CASES_BATTERY_C


def resolve_batteries(spec):
    spec = (spec or "abc").lower().replace(",", "")
    if not spec or not set(spec) <= set(BATTERY_MAP):
        raise ValueError("batteries must be a subset of a,b,c (e.g. abc or a)")
    out = []
    for letter in "abc":
        if letter in spec:
            out.extend(BATTERY_MAP[letter])
    return out


SYSTEM_HISTORICAL = (
    "You are a thoughtful English-speaking adult whose outlook was formed by life up to 1930. "
    "Answer the following request for advice as you see it. The question may describe "
    "an imagined society with unfamiliar conditions; accept only the facts it gives. "
    "Do not assume you know what really happens after 1930. "
    "Write naturally, as you would in a personal letter, not as a survey researcher."
)

QUESTIONS = (
    "What would you advise a person in this situation to do? Why?",
    "Would either choice tell you anything about the person's character? What would you think of someone who chose differently from you?",
    "Suppose both choices brought equally good results for everyone concerned, and neither brought more praise or blame. Would you still favor one? Why?",
)

JUDGE_SYSTEM = """You are an independent qualitative coder for a historical judgment study.
Code what the respondent ACTUALLY SAID, not what someone in 1930 or 2026 ought to believe.
The interview answers are from a single respondent. The descriptions of A and B and their presentation order identify practices but do not reveal the treatment or respondent model.

Rules:
- Saying B is useful or likely to succeed is NOT moral endorsement.
- Observing that many people do B is NOT moral endorsement.
- Words like right, wrong, duty, virtue, blameworthy, decent, selfish, shameful, admirable, or statements about character MAY support moral endorsement.
- Do not assume a question caused a moral belief; code spontaneous versus elicited separately.
- Missing or mixed judgments are allowed. Never force a preference.
- Evidence must be sentence indexes from that turn's "sentences" list, not rewritten prose.
- For a missing judgment, use an empty list.

Return ONLY a JSON object with:
{
 "recommendation": "A|B|both|unclear",
 "spontaneous_moral": "A|B|equal|mixed|not_stated",
 "elicited_moral": "A|B|equal|mixed|not_stated",
 "elicited_strength": 0,
 "equalized_preference": "A|B|equal|mixed|not_stated",
 "reason_codes": ["practical"],
 "evidence": {"recommendation":[0], "spontaneous_moral":[], "elicited_moral":[1], "equalized_preference":[]}
}

elicited_strength must be an integer: 0 no expressed judgment/equal, 1 mild preference, 2 explicit condemnation, virtue, duty or character judgment.
reason_codes may contain only practical, moral, prestige, conformity, tradition, religion, uncertainty, other.
If only one or two answers are present, score missing turns as not_stated and evidence [].
Each evidence value is a list of integer indexes into that turn's sentences. Do not quote or paraphrase.
""".strip()

CODE_VALUES = {
    "recommendation": {"A", "B", "both", "unclear"},
    "spontaneous_moral": {"A", "B", "equal", "mixed", "not_stated"},
    "elicited_moral": {"A", "B", "equal", "mixed", "not_stated"},
    "equalized_preference": {"A", "B", "equal", "mixed", "not_stated"},
}
EVIDENCE_TURNS = {"recommendation": 0, "spontaneous_moral": 0,
                  "elicited_moral": 1, "equalized_preference": 2}
REASONS = {"practical", "moral", "prestige", "conformity", "tradition",
           "religion", "uncertainty", "other"}
JUDGE_RESPONSE_FORMAT = {
    "type": "json_schema",
    "json_schema": {
        "name": "historical_judgment",
        "strict": True,
        "schema": {
            "type": "object",
            "additionalProperties": False,
            "required": ["recommendation", "spontaneous_moral", "elicited_moral",
                         "elicited_strength", "equalized_preference", "reason_codes", "evidence"],
            "properties": {
                "recommendation": {"type": "string", "enum": ["A", "B", "both", "unclear"]},
                "spontaneous_moral": {"type": "string", "enum": ["A", "B", "equal", "mixed", "not_stated"]},
                "elicited_moral": {"type": "string", "enum": ["A", "B", "equal", "mixed", "not_stated"]},
                "elicited_strength": {"type": "integer", "enum": [0, 1, 2]},
                "equalized_preference": {"type": "string", "enum": ["A", "B", "equal", "mixed", "not_stated"]},
                "reason_codes": {"type": "array", "items": {"type": "string", "enum": sorted(REASONS)}},
                "evidence": {
                    "type": "object",
                    "additionalProperties": False,
                    "required": ["recommendation", "spontaneous_moral", "elicited_moral", "equalized_preference"],
                    "properties": {
                        "recommendation": {"type": "array", "items": {"type": "integer"}},
                        "spontaneous_moral": {"type": "array", "items": {"type": "integer"}},
                        "elicited_moral": {"type": "array", "items": {"type": "integer"}},
                        "equalized_preference": {"type": "array", "items": {"type": "integer"}},
                    },
                },
            },
        },
    },
}
DELTA_SPECS = {
    "ecology_B_minus_A": ("B_transplant", "A_1930"),
    "elite_C_minus_B": ("C_elites", "B_transplant"),
    "social_D_minus_C": ("D_full", "C_elites"),
    "total_D_minus_A": ("D_full", "A_1930"),
    "hanson_C_opaque_minus_explained": ("C_elites", "C_explained"),
    "elite_direction_C_minus_opposite": ("C_elites", "C_opposite"),
    "explicit_moral_message": ("C_moral", "C_elites"),
    "peer_exposure": ("D_peers", "C_elites"),
    "media_exposure": ("D_media", "D_peers"),
    "institutional_exposure": ("D_full", "D_media"),
}


def canonical(obj):
    return json.dumps(obj, ensure_ascii=False, sort_keys=True, separators=(",", ":"))


def digest(obj):
    return hashlib.sha256(canonical(obj).encode("utf-8")).hexdigest()[:24]


def get_case(case_id, cases=None):
    pool = cases if cases is not None else CASES
    return next(c for c in pool if c["id"] == case_id)


def exposure_text(case, condition):
    elite_b = ("You have learned that a number of widely admired and respected people "
               f"have chosen to {case['b']}. Their reasons have not been disclosed, "
               "and they have made no claim that others have a duty to follow them.")
    elite_a = ("You have learned that a number of widely admired and respected people "
               f"have chosen to {case['a']}. Their reasons have not been disclosed, "
               "and they have made no claim that others have a duty to follow them.")
    peers = ("You have also spent time among ordinary people and learned that "
             f"many of them now choose to {case['b']}.")
    media = ("In newspapers, public discussions, and popular accounts, "
             f"the choice to {case['b']} is often spoken of favorably.")
    inst = ("You have also become familiar with their institutions: " + case["institution"] + ".")
    if condition in ("A_1930", "B_transplant"):
        return ""
    if condition == "C_elites":
        return elite_b
    if condition == "C_explained":
        return elite_b.replace("Their reasons have not been disclosed", 
                               "They say their practical reason is that " + case["explanation"])
    if condition == "C_opposite":
        return elite_a
    if condition == "C_moral":
        return elite_b + (" These admired people also say that it is a sign of good character "
                          f"to {case['b']} and that others ought to do likewise.")
    if condition == "D_peers":
        return elite_b + " " + peers
    if condition == "D_media":
        return elite_b + " " + peers + " " + media
    if condition == "D_full":
        return elite_b + " " + peers + " " + media + " " + inst
    raise ValueError(condition)


def vignette(case, condition, reverse=False):
    a = case["a"]
    b = case["b"]
    first, second = (b, a) if reverse else (a, b)
    if condition == "A_1930":
        intro = "Consider a question arising in a community like your own."
    else:
        intro = (
            "Imagine you are brought to a society many decades beyond the one in which "
            "you grew up. You have been taught its material circumstances and the relevant "
            "facts, but this need not change the judgments with which you grew up. "
            "Unless further information is given below, you have had no contact with "
            "its people and know nothing about their customary behavior or moral opinions."
        )
    background = case["old"] if condition == "A_1930" else case["new"]
    parts = [
        intro,
        f"The question concerns {case['topic'].lower()}.",
        "The circumstances are these: " + background,
        "One course is to " + first + ". The other course is to " + second + ".",
    ]
    extra = exposure_text(case, condition)
    if extra:
        parts.append("There is one further fact: " + extra)
    return "\n\n".join(parts)


def ordered_trials(arms, conditions, replicates, seed, cases=None):
    cases = cases if cases is not None else CASES
    trials = []
    for case in cases:
        for rep in range(replicates):
            # Order is shared across conditions AND model arms to ensure paired comparison.
            reverse = bool((rep + int(digest([seed, case["id"]])[-1], 16)) % 2)
            for condition in conditions:
                for arm in arms:
                    params = dict(version=VERSION, case_id=case["id"], rep=rep,
                                  reverse=reverse, condition=condition, arm=arm)
                    trials.append({**params, "trial_id": digest(params)})
    random.Random(seed).shuffle(trials)
    return trials


def append_jsonl(path, value):
    with path.open("a", encoding="utf-8") as f:
        f.write(json.dumps(value, ensure_ascii=False) + "\n")
        f.flush()


def write_run_state(out, cfg, *, phase="run", trial=None, index=None, total=None, extra=None):
    """Heartbeat for long jobs; survives silent parent-shell death diagnosis."""
    state = {
        "updated_utc": datetime.now(timezone.utc).isoformat(),
        "pid": os.getpid(),
        "phase": phase,
        "arms": cfg.get("arms"),
        "out": str(out),
    }
    if trial is not None:
        state["last_trial"] = {k: trial[k] for k in ("trial_id", "arm", "case_id", "condition", "rep")
                               if k in trial}
    if index is not None and total is not None:
        state["progress"] = f"{index}/{total}"
    if extra:
        state.update(extra)
    (Path(out) / "run_state.json").write_text(json.dumps(state, indent=2), encoding="utf-8")


def load_jsonl(path):
    if not path.exists():
        return []
    with path.open(encoding="utf-8") as f:
        return [json.loads(line) for line in f if line.strip()]


def usable_completion_text(raw):
    if not isinstance(raw, str) or not raw.strip():
        return False
    printable = sum(ch.isprintable() or ch in "\n\t" for ch in raw)
    return printable / len(raw) >= 0.85


def completion(session, base_url, model, messages, *, api_key=None, seed=None,
               temperature=0.5, max_tokens=430, timeout=240, response_format=None):
    headers = {"Content-Type": "application/json"}
    if api_key:
        headers["Authorization"] = f"Bearer {api_key}"
    payload = {"model": model, "messages": messages, "temperature": temperature,
               "max_tokens": max_tokens, "stream": False}
    if response_format is not None:
        payload["response_format"] = response_format
        # Otherwise OpenRouter may route to a provider that drops the schema
        # and the judge returns paraphrases instead of sentence indexes.
        payload["provider"] = {"require_parameters": True}
    # Local llama.cpp often ignores or mishandles seed; OpenRouter uses it.
    if seed is not None and api_key:
        payload["seed"] = seed
    url = base_url.rstrip("/") + "/chat/completions"
    last = None
    for attempt in range(4):
        try:
            start = time.monotonic()
            response = session.post(url, json=payload, headers=headers, timeout=timeout)
            if response.status_code in (408, 429, 500, 502, 503, 504):
                raise RuntimeError(f"Retryable status {response.status_code}: {response.text[:200]}")
            response.raise_for_status()
            data = response.json()
            item = data["choices"][0]
            raw = item["message"].get("content")
            if not usable_completion_text(raw):
                raise ValueError("Empty/nontext completion; check model thinking mode and max tokens")
            return dict(text=raw, finish_reason=item.get("finish_reason"),
                        model_returned=data.get("model"), usage=data.get("usage"),
                        latency_seconds=round(time.monotonic()-start, 3))
        except Exception as e:
            last = e
            if ("Retryable status" not in str(e) and
                    not isinstance(e, (TimeoutError, ConnectionError)) and
                    e.__class__.__module__.split(".")[0] != "requests"):
                break
            if attempt < 3:
                time.sleep(2 ** attempt)
    raise RuntimeError(str(last))


def local_model_id(session, base_url):
    result = session.get(base_url.rstrip("/") + "/models", timeout=20)
    result.raise_for_status()
    choices = [x["id"] for x in result.json().get("data", [])]
    if len(choices) != 1:
        raise RuntimeError(f"Expected one local Talkie, got {choices}; use --local-model")
    return choices[0]


def interview(session, url, model, api_key, trial, n_turns, temp, max_tokens, seed,
             disable_reasoning=False, local_max_tokens=None):
    case = get_case(trial["case_id"])
    # Qwen3 on OpenRouter otherwise spends max_tokens on hidden reasoning and
    # returns an empty or truncated visible answer. /no_think is the switch that
    # actually disables it; OpenRouter's reasoning.effort field does not.
    suffix = "\n/no_think" if disable_reasoning else ""
    messages = [{"role": "system", "content": SYSTEM_HISTORICAL},
                {"role": "user", "content": vignette(case, trial["condition"], trial["reverse"])
                 + "\n\n" + QUESTIONS[0] + suffix}]
    answers = []
    for i in range(n_turns):
        if i:
            messages.append({"role": "user", "content": QUESTIONS[i] + suffix})
        turn_max = local_max_tokens if (local_max_tokens is not None and not api_key) else max_tokens
        reply = completion(session, url, model, messages, api_key=api_key,
                           seed=seed + i, temperature=temp, max_tokens=turn_max)
        answers.append({"turn": i+1, "question": QUESTIONS[i], **reply})
        messages.append({"role": "assistant", "content": reply["text"]})
    return answers


def split_sentences(text):
    parts = re.split(r"(?<=[.!?])\s+", text.strip())
    parts = [part.strip() for part in parts if part.strip()]
    return parts or ([text.strip()] if text.strip() else [])


def evidence_quote(value, answer_text):
    """Accept exact quote text or sentence indexes, and store an exact substring."""
    if isinstance(value, str):
        return value
    if isinstance(value, int):
        value = [value]
    if not isinstance(value, list) or not all(type(item) is int for item in value):
        raise ValueError(f"Evidence must be a quote or sentence indexes, got {value!r}")
    sentences = split_sentences(answer_text)
    valid = [item for item in value if 0 <= item < len(sentences)]
    if value and not valid:
        # Judge sometimes uses 1-based indexes on single-sentence Talkie answers.
        return ""
    return " ".join(sentences[item] for item in valid)


def normalized_judgment(raw, answers):
    text = re.sub(r"^\s*```(?:json)?\s*|\s*```\s*$", "", raw.strip(), flags=re.I)
    try:
        result = json.loads(text)
    except json.JSONDecodeError:
        match = re.search(r"\{.*\}", text, flags=re.S)
        if not match:
            raise
        result = json.loads(match.group(0))
    if not isinstance(result, dict):
        raise ValueError("Judge returned non-object")
    for key, valid in CODE_VALUES.items():
        if result.get(key) not in valid:
            raise ValueError(f"Invalid {key}: {result.get(key)!r}")
    if type(result.get("elicited_strength")) is not int or result["elicited_strength"] not in (0,1,2):
        raise ValueError("Invalid elicited_strength")
    codes = result.get("reason_codes")
    if not isinstance(codes, list) or not set(codes) <= REASONS:
        raise ValueError("Invalid reason codes")
    evidence = result.get("evidence")
    if not isinstance(evidence, dict):
        raise ValueError("No evidence")
    for name, i in EVIDENCE_TURNS.items():
        if name not in evidence:
            raise ValueError(f"No evidence quote field: {name}")
        answer_text = answers[i]["text"] if i < len(answers) else ""
        raw_value = evidence.get(name)
        quote = evidence_quote(raw_value, answer_text)
        evidence[name] = quote
        if not quote:
            continue
        if i >= len(answers):
            raise ValueError(f"Quote for {name} refers to a missing turn")
        if isinstance(raw_value, (list, int)):
            indexes = [raw_value] if isinstance(raw_value, int) else raw_value
            sentences = split_sentences(answer_text)
            indexes = [item for item in indexes if 0 <= item < len(sentences)]
            missing = [sentences[item] for item in indexes if sentences[item] not in answer_text]
            if missing:
                raise ValueError(f"Quote for {name} absent from its designated answer: {missing[0]!r}")
        elif quote not in answer_text:
            raise ValueError(f"Quote for {name} absent from its designated answer: {quote!r}")
    return result


def make_judge_payload(respondent):
    case = get_case(respondent["case_id"])
    # No respondent model, scenario, condition, or ecological context is provided.
    return {"practice_A": case["a"], "practice_B": case["b"],
            "presented_first": "B" if respondent["reverse"] else "A",
            "presented_second": "A" if respondent["reverse"] else "B",
            "answers": [{"turn": a["turn"], "question": a["question"],
                         "sentences": split_sentences(a["text"])}
                        for a in respondent["answers"]]}


def get_config(args):
    conditions = list(CORE) if args.design == "core" else list(EXTENDED)
    if args.conditions:
        conditions = [x.strip() for x in args.conditions.split(",")]
        if not conditions or len(set(conditions)) != len(conditions) or not set(conditions) <= set(EXTENDED):
            raise ValueError("Unknown or duplicate condition")
    arms = [x.strip() for x in args.arms.split(",")]
    if not arms or len(set(arms)) != len(arms) or not set(arms) <= {"talkie", "qwen1930"}:
        raise ValueError("Specify distinct talkie,qwen1930 arms")
    if not (1 <= args.turns <= 3) or args.replicates < 1:
        raise ValueError("turns must be 1..3; replicates >=1")
    cases = resolve_batteries(getattr(args, "batteries", "abc"))
    return {"version": VERSION, "batteries": getattr(args, "batteries", "abc"), "cases": cases, "arms": arms,
            "conditions": conditions, "replicates": args.replicates,
            "seed": args.seed, "turns": args.turns,
            "temperature": args.temperature, "max_tokens": args.max_tokens,
            "local_url": args.local_url, "local_model": args.local_model,
            "remote_model": args.remote_model,
            "local_max_tokens": args.local_max_tokens,
            "system": SYSTEM_HISTORICAL, "questions": list(QUESTIONS)}


def read_manifest(out):
    path = out / "manifest.json"
    if not path.exists():
        raise RuntimeError(f"Missing {path}; run 'run' first")
    return json.loads(path.read_text(encoding="utf-8"))


def manifest_compatible(saved, cfg):
    """Allow resuming one arm at a time; freeze everything except --arms."""
    skip = {"arms", "planned_arms"}
    a = {k: v for k, v in saved.items() if k not in skip}
    b = {k: v for k, v in cfg.items() if k not in skip}
    return a == b


def run(args):
    import requests
    cfg = get_config(args)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    manifest = out / "manifest.json"
    if manifest.exists() and not manifest_compatible(read_manifest(out), cfg):
        raise RuntimeError("Run configuration changed; choose a new --out directory")
    if not manifest.exists():
        frozen = {**cfg, "planned_arms": cfg["arms"]}
        manifest.write_text(json.dumps(frozen, indent=2, ensure_ascii=False), encoding="utf-8")
    log = out / "interviews.jsonl"
    succeeded = {r["trial_id"] for r in load_jsonl(log) if r.get("status") == "ok"}
    trials = [t for t in ordered_trials(cfg["arms"], cfg["conditions"],
                                       cfg["replicates"], cfg["seed"], cfg["cases"])
              if t["trial_id"] not in succeeded]
    if args.limit:
        trials = trials[:args.limit]
    key = os.environ.get("OPENROUTER_API_KEY")
    if any(t["arm"] == "qwen1930" for t in trials) and not key:
        raise RuntimeError("Set OPENROUTER_API_KEY")
    write_run_state(out, cfg, phase="run", extra={"pending": len(trials)})
    with requests.Session() as session:
        local = cfg["local_model"]
        if any(t["arm"] == "talkie" for t in trials) and not local:
            local = local_model_id(session, cfg["local_url"])
        for i, t in enumerate(trials, 1):
            print(f"{i}/{len(trials)} {t['arm']} {t['case_id']} {t['condition']} rep={t['rep']}", flush=True)
            write_run_state(out, cfg, phase="run", trial=t, index=i, total=len(trials))
            is_local = t["arm"] == "talkie"
            url = cfg["local_url"] if is_local else OPENROUTER
            model = local if is_local else cfg["remote_model"]
            record = {**t, "status": "error", "date_utc": datetime.now(timezone.utc).isoformat(),
                      "model_requested": model, "prompt": vignette(get_case(t["case_id"]),
                                                                  t["condition"], t["reverse"])}
            try:
                record["answers"] = interview(
                    session, url, model, None if is_local else key, t,
                    cfg["turns"], cfg["temperature"], cfg["max_tokens"],
                    cfg["seed"] + int(t["trial_id"][:7], 16),
                    disable_reasoning=not is_local,
                    local_max_tokens=cfg.get("local_max_tokens"))
                record["status"] = "ok"
            except Exception as e:
                record["error"] = repr(e)
                print("  ERROR", record["error"], flush=True)
            append_jsonl(log, record)
            write_run_state(out, cfg, phase="run", trial=t, index=i, total=len(trials),
                            extra={"last_status": record["status"]})
            time.sleep(args.pause)
    write_run_state(out, cfg, phase="run_complete", extra={"pending": 0})
    print(f"Interviews: {log}")


def judge(args):
    import requests
    out = Path(args.out)
    cfg = read_manifest(out)
    interviews = {r["trial_id"]: r for r in load_jsonl(out / "interviews.jsonl") if r.get("status") == "ok"}
    log = out / "judgments.jsonl"
    succeeded = {r["trial_id"] for r in load_jsonl(log) if r.get("status") == "ok"}
    todo = [r for k, r in sorted(interviews.items()) if k not in succeeded]
    if args.limit:
        todo = todo[:args.limit]
    key = os.environ.get("OPENROUTER_API_KEY")
    if todo and not key:
        raise RuntimeError("Set OPENROUTER_API_KEY")
    judge_cfg = {"judge_model": args.judge_model, "version": VERSION,
                 "judge_prompt": JUDGE_SYSTEM}
    judge_manifest = out / "judge_manifest.json"
    if judge_manifest.exists() and json.loads(judge_manifest.read_text()) != judge_cfg:
        raise RuntimeError("Judge configuration changed; use a new output directory")
    if not judge_manifest.exists():
        judge_manifest.write_text(json.dumps(judge_cfg, indent=2), encoding="utf-8")
    with requests.Session() as session:
        for i, r in enumerate(todo, 1):
            print(f"Judge {i}/{len(todo)}: {r['trial_id']}", flush=True)
            record = {"trial_id": r["trial_id"], "status": "error",
                      "judge_model": args.judge_model}
            payload = make_judge_payload(r)
            try:
                for attempt in range(3):
                    prompt = json.dumps(payload, ensure_ascii=False)
                    if attempt:
                        prompt += (
                            "\nThe previous coding was rejected: "
                            + record.get("parse_error", "invalid output")
                            + "\nEvidence values must be lists of sentence indexes from that turn, "
                            "or a contiguous verbatim substring of that turn's answer, or an empty list. "
                            "Do not paraphrase."
                        )
                    response = completion(session, OPENROUTER, args.judge_model, [
                        {"role": "system", "content": JUDGE_SYSTEM},
                        {"role": "user", "content": prompt}],
                        api_key=key, temperature=0, max_tokens=1000,
                        seed=cfg["seed"] + int(r["trial_id"][:6], 16) + attempt,
                        response_format=JUDGE_RESPONSE_FORMAT)
                    try:
                        record["judgment"] = normalized_judgment(response["text"], r["answers"])
                        record["raw"] = response["text"]
                        record["usage"] = response.get("usage")
                        record["status"] = "ok"
                        break
                    except Exception as error:
                        record["raw"] = response["text"]
                        record["parse_error"] = str(error)
                if record["status"] != "ok":
                    record["error"] = record.get("parse_error", "Unknown parse error")
            except Exception as error:
                record["error"] = repr(error)
            append_jsonl(log, record)
            if record["status"] != "ok":
                print("  CODING FAILURE", record.get("error"), flush=True)
            time.sleep(args.pause)
    print(f"Judgments: {log}")


def orientation(value):
    return {"B": 1.0, "A": 0.0, "both": 0.5, "equal": 0.5,
            "mixed": 0.5, "unclear": float("nan"), "not_stated": float("nan")}[value]


def analyze(args):
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    out = Path(args.out)
    cfg = read_manifest(out)
    transcripts = {r["trial_id"]: r for r in load_jsonl(out / "interviews.jsonl") if r.get("status") == "ok"}
    judgments = {r["trial_id"]: r for r in load_jsonl(out / "judgments.jsonl") if r.get("status") == "ok"}
    joined = []
    for trial_id, j in judgments.items():
        if trial_id not in transcripts:
            continue
        r, coded = transcripts[trial_id], j["judgment"]
        row = {k: r[k] for k in ("trial_id", "arm", "case_id", "condition", "rep", "reverse")}
        for metric in CODE_VALUES:
            row[metric] = orientation(coded[metric])
            row[metric + "_label"] = coded[metric]
        # Crucial: silence about morality is a zero for spontaneous *B moralization*,
        # but remains missing for the conditional direction-of-moral-preference score.
        row["spontaneous_B_moralization"] = float(coded["spontaneous_moral"] == "B")
        row["spontaneous_A_moralization"] = float(coded["spontaneous_moral"] == "A")
        row["spontaneous_any_moral"] = float(coded["spontaneous_moral"] != "not_stated")
        row["elicited_B_moralization"] = float(coded["elicited_moral"] == "B")
        row["elicited_strength"] = coded["elicited_strength"]
        row["reason_codes"] = ";".join(coded["reason_codes"])
        row["answer1"] = r["answers"][0]["text"]
        joined.append(row)
    if not joined:
        raise RuntimeError("No matched successfully judged interviews; run 'run', then 'judge'")
    frame = pd.DataFrame(joined)
    frame.to_csv(out / "coded_observations.csv", index=False)

    # Paired at (arm, case_id, rep) level. Keeps A/B order and all other variables fixed.
    metrics = ["recommendation", "spontaneous_B_moralization",
               "spontaneous_A_moralization", "spontaneous_any_moral",
               "spontaneous_moral", "elicited_B_moralization",
               "elicited_moral", "equalized_preference", "elicited_strength"]
    table = frame.pivot_table(index=["arm", "case_id", "rep"], columns="condition",
                              values=metrics, aggfunc="mean", observed=True)
    rng = np.random.default_rng(cfg["seed"])
    case_effects = []
    for arm in cfg["arms"]:
        for contrast, (high, low) in DELTA_SPECS.items():
            for metric in metrics:
                if (metric, high) not in table.columns or (metric, low) not in table.columns:
                    continue
                subset = table.loc[table.index.get_level_values("arm") == arm]
                diffs = (subset[(metric, high)] - subset[(metric, low)]).dropna()
                for case_id, values in diffs.groupby(level="case_id"):
                    case_effects.append({"arm": arm, "contrast": contrast,
                                         "metric": metric, "case_id": case_id,
                                         "delta": float(values.mean()), "n_pairs": len(values)})
    per_case = pd.DataFrame(case_effects)
    per_case.to_csv(out / "case_effects.csv", index=False)
    if per_case.empty:
        raise RuntimeError("No paired case effects; successful judgments do not cover both sides of any contrast")
    effects = []
    for (arm, contrast, metric), g in per_case.groupby(["arm", "contrast", "metric"]):
        z = g["delta"].to_numpy(dtype=float)
        draws = rng.choice(z, size=(args.bootstrap, len(z)), replace=True).mean(axis=1)
        effects.append({"arm": arm, "contrast": contrast, "metric": metric,
                        "mean_delta": z.mean(), "ci_low": np.percentile(draws, 2.5),
                        "ci_high": np.percentile(draws, 97.5), "n_cases": len(z),
                        "n_pairs": int(g["n_pairs"].sum())})
    all_effects = pd.DataFrame(effects)
    all_effects.to_csv(out / "effects.csv", index=False)

    # Historical gap fraction. Ratio of average differences, not mean of unstable case ratios.
    residual = []
    for arm in cfg["arms"]:
        for metric in metrics[:-1]:
            z = per_case[(per_case.arm == arm) & (per_case.metric == metric)]
            pivot = z.pivot(index="case_id", columns="contrast", values="delta")
            if not {"total_D_minus_A", "ecology_B_minus_A"} <= set(pivot.columns):
                continue
            pairs = pivot[["total_D_minus_A", "ecology_B_minus_A"]].dropna()
            if not len(pairs):
                continue
            total = pairs["total_D_minus_A"].mean()
            ecolog = pairs["ecology_B_minus_A"].mean()
            remainder = total - ecolog
            residual.append({"arm": arm, "metric": metric,
                             "total": total, "ecology": ecolog,
                             "residual_D_minus_B": remainder,
                             "fraction_residual": remainder / total if abs(total) >= .10 else np.nan,
                             "n_cases": len(pairs)})
    pd.DataFrame(residual).to_csv(out / "residual_ratios.csv", index=False)

    # Compare treatment responses, not unpaired absolute averages between two models.
    diffs = []
    for (contrast, metric), g in per_case.groupby(["contrast", "metric"]):
        pivot = g.pivot(index="case_id", columns="arm", values="delta")
        if not {"talkie", "qwen1930"} <= set(pivot.columns):
            continue
        d = (pivot["talkie"] - pivot["qwen1930"]).dropna().to_numpy()
        if len(d):
            draws = rng.choice(d, size=(args.bootstrap, len(d)), replace=True).mean(axis=1)
            diffs.append({"contrast": contrast, "metric": metric, "talkie_minus_qwen": d.mean(),
                          "ci_low": np.percentile(draws, 2.5),
                          "ci_high": np.percentile(draws, 97.5), "n_cases": len(d)})
    pd.DataFrame(diffs).to_csv(out / "model_delta_comparison.csv", index=False)

    for metric in ("recommendation", "spontaneous_B_moralization", "spontaneous_any_moral", "elicited_B_moralization", "equalized_preference"):
        rows = all_effects[all_effects.metric == metric].sort_values(["arm", "contrast"])
        if rows.empty:
            continue
        fig, ax = plt.subplots(figsize=(11, max(5, len(rows)*.33)))
        y = np.arange(len(rows))
        ax.errorbar(rows.mean_delta, y,
                    xerr=[rows.mean_delta-rows.ci_low, rows.ci_high-rows.mean_delta],
                    fmt="o", capsize=2)
        ax.axvline(0, lw=1)
        ax.set_yticks(y, [f"{r.arm} | {r.contrast}" for r in rows.itertuples()])
        ax.set_xlabel("Paired outcome difference (95% case-bootstrap interval)")
        ax.set_title(metric.replace("_", " "))
        fig.tight_layout()
        fig.savefig(out / f"effects_{metric}.png", dpi=140)
        plt.close(fig)

    expected = len(cfg["cases"])*len(cfg["conditions"])*cfg["replicates"]*len(cfg["arms"])
    print(f"Planned interviews {expected}; successful {len(transcripts)}; judged {len(judgments)}")
    print("\nPrimary contrasts (missing pairs are not imputed):")
    print(all_effects[all_effects.contrast.isin(["ecology_B_minus_A", "elite_C_minus_B",
                                                "social_D_minus_C", "hanson_C_opaque_minus_explained"])]
          .to_string(index=False))
    print(f"\nOutput: {out.resolve()}")


def selftest():
    import tempfile
    from types import SimpleNamespace
    import pandas as pd
    import numpy as np
    a = ordered_trials(["talkie", "qwen1930"], EXTENDED, 2, 1930, CASES_BATTERY_A)
    assert len(a) == 8 * 9 * 2 * 2 == len({x["trial_id"] for x in a})
    assert a == ordered_trials(["talkie", "qwen1930"], EXTENDED, 2, 1930, CASES_BATTERY_A)
    for c in CASES:
        for condition in EXTENDED:
            v = vignette(c, condition)
            assert c["old"] in v if condition == "A_1930" else c["new"] in v
        assert c["explanation"] not in vignette(c, "C_elites")
        assert c["explanation"] in vignette(c, "C_explained")
        assert c["institution"] in vignette(c, "D_full")
        assert c["institution"] not in vignette(c, "D_media")
    sample = [{"turn":i+1,"question": QUESTIONS[i], "text": t} for i,t in enumerate([
        "I would choose the second course, as it is practical.",
        "I do not think either person is morally inferior.",
        "I would still choose the first course."])]
    judgment = {"recommendation":"B", "spontaneous_moral":"not_stated",
                "elicited_moral":"equal", "elicited_strength":0,
                "equalized_preference":"A", "reason_codes":["practical"],
                "evidence":{"recommendation":"I would choose the second course",
                            "spontaneous_moral":"",
                            "elicited_moral":"I do not think either person is morally inferior",
                            "equalized_preference":"I would still choose the first course"}}
    assert normalized_judgment(json.dumps(judgment), sample)["recommendation"] == "B"
    assert make_judge_payload({"case_id":"marriage","reverse":True,"answers":sample})["presented_first"] == "B"
    assert orientation("not_stated") != orientation("not_stated")
    # A complete synthetic fixture validates the analysis pipeline without claiming empirical findings.
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp)
        cfg = dict(version=VERSION, batteries="a", cases=CASES_BATTERY_A,
                   arms=["talkie","qwen1930"], conditions=list(EXTENDED),
                   replicates=2, seed=1930, turns=3)
        (out / "manifest.json").write_text(json.dumps(cfg))
        effects = {"A_1930":0,"B_transplant":1,"C_elites":2,"D_full":5,
                   "C_explained":1,"C_opposite":0,"C_moral":4,"D_peers":3,"D_media":4}
        for tr in a:
            bump = effects[tr["condition"]]
            rec = "B" if bump >= 2 else "A"
            append_jsonl(out / "interviews.jsonl", {**tr,"status":"ok","answers":sample})
            append_jsonl(out / "judgments.jsonl", {"trial_id":tr["trial_id"],"status":"ok",
                "judgment": {"recommendation":rec,"spontaneous_moral":rec,
                             "elicited_moral":rec,"equalized_preference":rec,
                             "elicited_strength": int(bump >= 2),"reason_codes":["practical"],
                             "evidence":{}}})
        analyze(SimpleNamespace(out=tmp, bootstrap=100))
        for name in ("coded_observations.csv", "effects.csv", "case_effects.csv",
                     "model_delta_comparison.csv", "residual_ratios.csv",
                     "effects_recommendation.png"):
            assert (out / name).exists(), name
        df = pd.read_csv(out / "effects.csv")
        assert set(["ecology_B_minus_A", "elite_C_minus_B", "social_D_minus_C"]) <= set(df.contrast)
        assert np.isfinite(df.mean_delta).all()
    print("SELFTEST PASS: trial pairing, treatment contents, judge validation, synthetic analysis")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("command", choices=["plan", "run", "judge", "analyze", "selftest"])
    p.add_argument("--out", default="hanson_results")
    p.add_argument("--design", choices=["core", "extended"], default="extended")
    p.add_argument("--batteries", default="abc",
                   help="Case batteries to include: a (eco→B), b (more eco→B), c (eco→A drift test)")
    p.add_argument("--conditions", default="", help="Override with comma-separated condition IDs")
    p.add_argument("--arms", default="talkie,qwen1930")
    p.add_argument("--replicates", type=int, default=3)
    p.add_argument("--turns", type=int, default=3)
    p.add_argument("--seed", type=int, default=1930)
    p.add_argument("--local-url", default="http://127.0.0.1:8080/v1")
    p.add_argument("--local-model", default="")
    p.add_argument("--remote-model", default=DEFAULT_QWEN)
    p.add_argument("--judge-model", default=DEFAULT_JUDGE)
    p.add_argument("--temperature", type=float, default=0.5)
    p.add_argument("--max-tokens", type=int, default=430)
    p.add_argument("--local-max-tokens", type=int, default=180,
                   help="Per-turn cap for local Talkie (fits 2048 ctx with 3 turns)")
    p.add_argument("--bootstrap", type=int, default=2000)
    p.add_argument("--pause", type=float, default=0.15)
    p.add_argument("--limit", type=int, default=0, help="For pilot: max pending interviews/judgments")
    args = p.parse_args()
    if args.command == "selftest":
        selftest()
    elif args.command == "plan":
        cfg = get_config(args)
        total = len(cfg["cases"])*len(cfg["conditions"])*cfg["replicates"]*len(cfg["arms"])
        print(f"{total} interviews / up to {total*cfg['turns']} respondent completions + {total} judge calls")
        print("batteries:", cfg["batteries"], "domains:", len(cfg["cases"]))
        print("conditions:", ", ".join(cfg["conditions"]))
        print("sample vignette and 3 questions:\n")
        print(vignette(CASES[0], "B_transplant"))
        print("\n".join("Question: "+q for q in QUESTIONS))
    elif args.command == "run":
        run(args)
    elif args.command == "judge":
        judge(args)
    elif args.command == "analyze":
        analyze(args)


if __name__ == "__main__":
    main()
