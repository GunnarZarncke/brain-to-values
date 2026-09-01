import torch
from torch import nn
from .env import generate_paired_batch, CueConfig


@torch.no_grad()
def task_accuracy(model, cue_config, n_worlds=2048, horizon=8, device="cpu", seed=1000, swapped=False):
    b = generate_paired_batch(n_worlds, horizon, cue_config, device=device, seed=seed)
    logits = model(b.visual, b.cues, swap_state_before_final=swapped)
    pred = logits.argmax(-1)
    target = b.label
    if swapped:
        target = target.view(-1, 2).flip(1).reshape(-1)
    return (pred == target).float().mean().item()


def _linear_probe(X, y, n_classes, steps=300, lr=0.05):
    n = X.shape[0]
    idx = torch.randperm(n, device=X.device)
    cut = int(0.7 * n)
    tr, te = idx[:cut], idx[cut:]
    probe = nn.Linear(X.shape[-1], n_classes).to(X.device)
    opt = torch.optim.Adam(probe.parameters(), lr=lr)
    for _ in range(steps):
        loss = nn.functional.cross_entropy(probe(X[tr]), y[tr])
        opt.zero_grad(); loss.backward(); opt.step()
    with torch.no_grad():
        return (probe(X[te]).argmax(-1) == y[te]).float().mean().item()


def self_probe_curves(model, cue_config, n_worlds=2048, horizon=8, device="cpu", seed=2000):
    b = generate_paired_batch(n_worlds, horizon, cue_config, device=device, seed=seed)
    with torch.no_grad():
        _, z, agents = model(b.visual, b.cues, return_traces=True)
    z_acc, obj_acc = [], []
    for t in range(z.shape[1]):
        z_acc.append(_linear_probe(z[:, t].detach(), b.self_id, 2))
        X = agents[:, t].reshape(-1, agents.shape[-1]).detach()
        y = torch.zeros(agents.shape[0], 2, dtype=torch.long, device=device)
        y[torch.arange(agents.shape[0], device=device), b.self_id] = 1
        obj_acc.append(_linear_probe(X, y.reshape(-1), 2))
    return z_acc, obj_acc
