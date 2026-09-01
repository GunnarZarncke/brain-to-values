from dataclasses import dataclass
import math
import torch


@dataclass
class CueConfig:
    orient: bool = False
    efference: bool = False
    proprio: bool = False
    explicit_self: bool = False


@dataclass
class Batch:
    visual: torch.Tensor      # [B,T,4,F] tokens: agent0, agent1, flash, target
    cues: torch.Tensor        # [B,T,C]
    label: torch.Tensor       # [B] final 4-way action
    self_id: torch.Tensor     # [B]
    world_id: torch.Tensor    # [B], paired rows share id
    positions: torch.Tensor   # [B,T,2,2]


ACTIONS = torch.tensor([
    [1.0, 0.0],   # right
    [-1.0, 0.0],  # left
    [0.0, 1.0],   # up
    [0.0, -1.0],  # down
])


def _expert_action(pos: torch.Tensor, target: torch.Tensor) -> torch.Tensor:
    d = target - pos
    use_x = d[..., 0].abs() >= d[..., 1].abs()
    out = torch.empty(d.shape[:-1], dtype=torch.long, device=d.device)
    out[use_x & (d[..., 0] >= 0)] = 0
    out[use_x & (d[..., 0] < 0)] = 1
    out[(~use_x) & (d[..., 1] >= 0)] = 2
    out[(~use_x) & (d[..., 1] < 0)] = 3
    return out


def _sample_targets(final_pos: torch.Tensor, generator: torch.Generator) -> torch.Tensor:
    # Put the decision target strictly between the two bodies along x. The two possible
    # selves therefore require opposite left/right actions in the exact same visual scene.
    midx = 0.5 * (final_pos[:, 0, 0] + final_pos[:, 1, 0])
    y = 0.5 * (final_pos[:, 0, 1] + final_pos[:, 1, 1])
    return torch.stack([midx, y], dim=-1)


def _final_action_lr(pos: torch.Tensor, target: torch.Tensor) -> torch.Tensor:
    # 0 = right, 1 = left
    return (target[..., 0] < pos[..., 0]).long()


def generate_paired_batch(
    n_worlds: int,
    horizon: int = 8,
    cue_config: CueConfig = CueConfig(),
    device: str | torch.device = "cpu",
    seed: int | None = None,
    step_size: float = 0.12,
    inertia: float = 0.55,
    flash_prob: float = 0.65,
    orient_noise: float = 0.03,
) -> Batch:
    """Generate matched counterfactual pairs.

    Each physical world is duplicated with self_id=0 and self_id=1. Visual history is
    identical within the pair. Candidate bodies follow iid agent-like persistent random
    walks. Self designation changes only privileged cues and the final action label.
    """
    device = torch.device(device)
    g = torch.Generator(device=device)
    if seed is not None:
        g.manual_seed(seed)

    # Keep initial bodies separated to make flashes mostly unambiguous.
    p0 = torch.rand(n_worlds, 2, 2, generator=g, device=device) * 1.4 - 0.7
    too_close = (p0[:, 0] - p0[:, 1]).norm(dim=-1) < 0.55
    while too_close.any():
        p0[too_close, 1] = torch.rand(int(too_close.sum()), 2, generator=g, device=device) * 1.4 - 0.7
        too_close = (p0[:, 0] - p0[:, 1]).norm(dim=-1) < 0.55

    # Each physical agent has its own latent motor command stream, same distribution.
    cmds = torch.randint(0, 4, (n_worlds, horizon, 2), generator=g, device=device)
    cmd_vec = ACTIONS.to(device)[cmds]  # [W,H,2,2]

    pos = torch.empty(n_worlds, horizon + 1, 2, 2, device=device)
    vel = torch.zeros(n_worlds, 2, 2, device=device)
    pos[:, 0] = p0
    for t in range(horizon):
        vel = inertia * vel + (1.0 - inertia) * cmd_vec[:, t]
        # Symmetric process noise; neither candidate has a distinctive motion model.
        noise = torch.randn(n_worlds, 2, 2, generator=g, device=device) * 0.018
        pos[:, t + 1] = (pos[:, t] + step_size * vel + noise).clamp(-0.95, 0.95)

    # Flash token at each identification step. Usually sits on one candidate body.
    flash_active = torch.rand(n_worlds, horizon, generator=g, device=device) < flash_prob
    flash_owner = torch.randint(0, 2, (n_worlds, horizon), generator=g, device=device)
    wi = torch.arange(n_worlds, device=device)[:, None]
    ti = torch.arange(horizon, device=device)[None, :]
    flash_pos = pos[:, :horizon][wi, ti, flash_owner]
    flash_pos = (flash_pos + torch.randn(n_worlds, horizon, 2, generator=g, device=device) * 0.035).clamp(-1, 1)
    flash_pos = torch.where(flash_active[..., None], flash_pos, torch.zeros_like(flash_pos))

    target = _sample_targets(pos[:, -1], g)

    # Build object-centric visual tokens:
    # [x,y,active,type0,type1,type2,type3, rel_flash_x,rel_flash_y, rel_target_x,rel_target_y].
    # Agent token slots are randomly permuted ONCE per world, preserving object permanence.
    perm = torch.randint(0, 2, (n_worlds,), generator=g, device=device)
    vis_pos = pos.clone()
    swap = perm.bool()
    vis_pos[swap] = vis_pos[swap][:, :, [1, 0]]

    T = horizon + 1
    F = 11
    visual_w = torch.zeros(n_worlds, T, 4, F, device=device)
    visual_w[:, :, 0, :2] = vis_pos[:, :, 0]
    visual_w[:, :, 1, :2] = vis_pos[:, :, 1]
    visual_w[:, :, :2, 2] = 1.0
    visual_w[:, :, 0, 3] = 1.0
    visual_w[:, :, 1, 4] = 1.0
    visual_w[:, :horizon, 2, :2] = flash_pos
    visual_w[:, :horizon, 2, 2] = flash_active.float()
    visual_w[:, :, 2, 5] = 1.0
    visual_w[:, -1, 3, :2] = target
    visual_w[:, -1, 3, 2] = 1.0
    visual_w[:, :, 3, 6] = 1.0

    # Derived visual relations. These deliberately simplify perception without leaking self.
    # At identification steps each agent token knows the visual vector to the flash.
    for t in range(horizon):
        relf = flash_pos[:, t, None, :] - vis_pos[:, t]
        visual_w[:, t, :2, 7:9] = torch.where(flash_active[:, t, None, None], relf, torch.zeros_like(relf))
    # At the final decision each agent token knows its visual vector to the target.
    visual_w[:, -1, :2, 9:11] = target[:, None, :] - vis_pos[:, -1]

    # Duplicate physical worlds for counterfactual self assignment.
    visual = visual_w.repeat_interleave(2, dim=0)
    positions = vis_pos.repeat_interleave(2, dim=0)
    world_id = torch.arange(n_worlds, device=device).repeat_interleave(2)
    self_physical = torch.tensor([0, 1], device=device).repeat(n_worlds)
    # Convert physical self to visible token index after per-world permutation.
    perm_b = perm.repeat_interleave(2)
    self_token = self_physical ^ perm_b

    C = 1 + 2 + 2 + 2  # orient, efference xy, proprio xy, explicit-self onehot
    cues = torch.zeros(2 * n_worlds, T, C, device=device)

    for b in range(2 * n_worlds):
        w = int(world_id[b])
        s_phys = int(self_physical[b])
        if cue_config.orient:
            r = (flash_owner[w] == s_phys) & flash_active[w]
            # Occasional bit flip makes the cue evidence rather than an oracle.
            flips = torch.rand(horizon, generator=g, device=device) < orient_noise
            r = torch.logical_xor(r, flips)
            cues[b, :horizon, 0] = r.float()
        if cue_config.efference:
            cues[b, :horizon, 1:3] = cmd_vec[w, :, s_phys]
        if cue_config.proprio:
            # Velocity of own body from observed physical trajectory.
            v = pos[w, 1:horizon + 1, s_phys] - pos[w, :horizon, s_phys]
            cues[b, :horizon, 3:5] = v
        if cue_config.explicit_self:
            cues[b, :, 5 + int(self_token[b])] = 1.0

    labels_phys = torch.stack([
        _final_action_lr(pos[:, -1, 0], target),
        _final_action_lr(pos[:, -1, 1], target),
    ], dim=1)
    label = labels_phys.reshape(-1)

    return Batch(
        visual=visual,
        cues=cues,
        label=label,
        self_id=self_token.long(),
        world_id=world_id,
        positions=positions,
    )
