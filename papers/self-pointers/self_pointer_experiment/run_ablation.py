import argparse
import json
from pathlib import Path
import torch
from torch.nn import functional as F

from self_pointer.env import CueConfig, generate_paired_batch
from self_pointer.model import ModelConfig, SelfPointerNet
from self_pointer.analysis import task_accuracy


CONDS = {
    "vision_only": CueConfig(),
    "orient_only": CueConfig(orient=True),
    "efference_only": CueConfig(efference=True),
    "proprio_only": CueConfig(proprio=True),
    "orient_efference": CueConfig(orient=True, efference=True),
    "all_internal": CueConfig(orient=True, efference=True, proprio=True),
    "explicit_self": CueConfig(explicit_self=True),
}


def train(cond, recirc, policy_reads_state, steps, batch_worlds, horizon, device, seed):
    torch.manual_seed(seed)
    model = SelfPointerNet(ModelConfig(recirculate=recirc, policy_reads_state=policy_reads_state)).to(device)
    opt = torch.optim.AdamW(model.parameters(), lr=3e-4)
    for s in range(steps):
        b = generate_paired_batch(batch_worlds, horizon, CONDS[cond], device=device, seed=seed+s)
        loss = F.cross_entropy(model(b.visual, b.cues), b.label)
        opt.zero_grad(); loss.backward(); opt.step()
    a = task_accuracy(model, CONDS[cond], 1024, horizon, device, seed+50000)
    ap = task_accuracy(model, CONDS[cond], 1024, horizon, device, seed+50000, swapped=True)
    return {"condition": cond, "recirculate": recirc, "policy_reads_state": policy_reads_state, "accuracy": a, "swap_accuracy": ap}


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--steps", type=int, default=1200)
    p.add_argument("--batch-worlds", type=int, default=96)
    p.add_argument("--horizon", type=int, default=8)
    p.add_argument("--device", default="cuda" if torch.cuda.is_available() else "cpu")
    p.add_argument("--out", default="results/ablation.json")
    a = p.parse_args()
    results=[]
    k=0
    for cond in CONDS:
        architectures = [(True, True)]
        if cond in {"orient_only", "all_internal"}:
            architectures = [
                (True, False),   # recirculation is the only route from persistent state to decision
                (False, True),   # standard recurrent state goes directly to policy
                (True, True),    # both routes available
            ]
        for recirc, policy_reads_state in architectures:
            r=train(cond, recirc, policy_reads_state, a.steps, a.batch_worlds, a.horizon, a.device, 10_000+k*1000)
            print(r); results.append(r); k += 1
    out=Path(a.out); out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(results, indent=2))

if __name__ == "__main__":
    main()
