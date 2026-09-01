from dataclasses import dataclass
import torch
from torch import nn


@dataclass
class ModelConfig:
    visual_dim: int = 11
    cue_dim: int = 7
    d_model: int = 48
    z_dim: int = 32
    nhead: int = 4
    recirculate: bool = True
    policy_reads_state: bool = True


class SelfPointerNet(nn.Module):
    """Small recurrent object-centric network with an explicit recirculation bottleneck."""

    def __init__(self, cfg: ModelConfig = ModelConfig()):
        super().__init__()
        self.cfg = cfg
        self.visual_in = nn.Linear(cfg.visual_dim, cfg.d_model)
        enc1 = nn.TransformerEncoderLayer(cfg.d_model, cfg.nhead, 2 * cfg.d_model, batch_first=True)
        self.encoder1 = nn.TransformerEncoder(enc1, num_layers=1)
        self.gru = nn.GRUCell(cfg.d_model + cfg.cue_dim, cfg.z_dim)
        self.z_to_token = nn.Linear(cfg.z_dim, cfg.d_model)
        self.token_gate = nn.Linear(cfg.d_model, cfg.d_model)
        enc2 = nn.TransformerEncoderLayer(cfg.d_model, cfg.nhead, 2 * cfg.d_model, batch_first=True)
        self.encoder2 = nn.TransformerEncoder(enc2, num_layers=1)
        policy_in = 2 * cfg.d_model + (cfg.z_dim if cfg.policy_reads_state else 0)
        self.policy = nn.Sequential(
            nn.Linear(policy_in, cfg.d_model),
            nn.GELU(),
            nn.Linear(cfg.d_model, 2),
        )

    def init_state(self, batch_size, device):
        return torch.zeros(batch_size, self.cfg.z_dim, device=device)

    def percept(self, visual_t, z):
        h = self.encoder1(self.visual_in(visual_t))
        if self.cfg.recirculate:
            # Generic multiplicative state-object interaction. The same recurrent state
            # can therefore alter different tokens differently based on their current
            # representation, without any hard-coded self slot.
            h = h + torch.tanh(self.token_gate(h)) * self.z_to_token(z)[:, None, :]
        h = self.encoder2(h)
        return h

    def forward(self, visual, cues, swap_state_before_final=False, return_traces=False):
        B, T, _, _ = visual.shape
        z = self.init_state(B, visual.device)
        z_trace, agent_trace = [], []

        for t in range(T):
            if swap_state_before_final and t == T - 1:
                if B % 2:
                    raise ValueError("state swapping expects paired rows")
                z = z.view(B // 2, 2, -1).flip(1).reshape(B, -1)

            h = self.percept(visual[:, t], z)
            pooled = h.mean(dim=1)
            z = self.gru(torch.cat([pooled, cues[:, t]], dim=-1), z)
            z_trace.append(z)
            agent_trace.append(h[:, :2])

        agents_flat = h[:, :2].reshape(B, -1)
        policy_input = torch.cat([agents_flat, z], dim=-1) if self.cfg.policy_reads_state else agents_flat
        logits = self.policy(policy_input)
        if return_traces:
            return logits, torch.stack(z_trace, 1), torch.stack(agent_trace, 1)
        return logits
