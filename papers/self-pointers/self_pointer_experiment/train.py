import argparse
import json
from pathlib import Path
import torch
from torch.nn import functional as F

from self_pointer.env import CueConfig, generate_paired_batch
from self_pointer.model import ModelConfig, SelfPointerNet
from self_pointer.analysis import task_accuracy, self_probe_curves


EXPERIMENTS = {
    "vision_only": CueConfig(),
    "orient_only": CueConfig(orient=True),
    "efference_only": CueConfig(efference=True),
    "proprio_only": CueConfig(proprio=True),
    "orient_efference": CueConfig(orient=True, efference=True),
    "orient_proprio": CueConfig(orient=True, proprio=True),
    "all_internal": CueConfig(orient=True, efference=True, proprio=True),
    "explicit_self": CueConfig(explicit_self=True),
}


def train_one(name, steps, batch_worlds, horizon, lr, device, seed, recirculate, outdir):
    torch.manual_seed(seed)
    cues = EXPERIMENTS[name]
    model = SelfPointerNet(ModelConfig(recirculate=recirculate)).to(device)
    opt = torch.optim.AdamW(model.parameters(), lr=lr)

    for step in range(steps):
        b = generate_paired_batch(batch_worlds, horizon, cues, device=device, seed=seed + step)
        logits = model(b.visual, b.cues)
        loss = F.cross_entropy(logits, b.label)
        opt.zero_grad(); loss.backward(); opt.step()
        if (step + 1) % max(1, steps // 10) == 0:
            acc = (logits.argmax(-1) == b.label).float().mean().item()
            print(f"{name:20s} step={step+1:5d} loss={loss.item():.4f} train_acc={acc:.3f}")

    eval_acc = task_accuracy(model, cues, n_worlds=1024, horizon=horizon, device=device, seed=seed + 100000)
    patch_acc = task_accuracy(model, cues, n_worlds=1024, horizon=horizon, device=device, seed=seed + 100000, swapped=True)
    z_curve, obj_curve = self_probe_curves(model, cues, n_worlds=1024, horizon=horizon, device=device, seed=seed + 200000)

    result = {
        "name": name,
        "recirculate": recirculate,
        "eval_accuracy": eval_acc,
        "state_swap_counterfactual_accuracy": patch_acc,
        "z_self_probe": z_curve,
        "object_self_probe": obj_curve,
    }
    suffix = "recirc" if recirculate else "no_recirc"
    torch.save(model.state_dict(), outdir / f"{name}_{suffix}.pt")
    (outdir / f"{name}_{suffix}.json").write_text(json.dumps(result, indent=2))
    print(json.dumps(result, indent=2))
    return result


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--experiment", choices=list(EXPERIMENTS) + ["all"], default="all")
    p.add_argument("--steps", type=int, default=2000)
    p.add_argument("--batch-worlds", type=int, default=128)
    p.add_argument("--horizon", type=int, default=8)
    p.add_argument("--lr", type=float, default=3e-4)
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--device", default="cuda" if torch.cuda.is_available() else "cpu")
    p.add_argument("--no-recirculation", action="store_true")
    p.add_argument("--outdir", default="results")
    args = p.parse_args()
    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)
    names = list(EXPERIMENTS) if args.experiment == "all" else [args.experiment]
    all_results = []
    for i, name in enumerate(names):
        all_results.append(train_one(name, args.steps, args.batch_worlds, args.horizon, args.lr,
                                     args.device, args.seed + 10000 * i, not args.no_recirculation, outdir))
    (outdir / "summary.json").write_text(json.dumps(all_results, indent=2))


if __name__ == "__main__":
    main()
