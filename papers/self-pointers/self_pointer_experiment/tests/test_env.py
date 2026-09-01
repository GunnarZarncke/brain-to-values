import torch
from self_pointer.env import CueConfig, generate_paired_batch


def test_counterfactual_pairs_share_visuals_and_need_different_actions():
    b = generate_paired_batch(32, cue_config=CueConfig(orient=True), seed=1)
    v = b.visual.view(32,2,*b.visual.shape[1:])
    assert torch.allclose(v[:,0], v[:,1])
    y = b.label.view(32,2)
    assert torch.all(y[:,0] != y[:,1])


def test_orient_cue_differs_with_self_on_informative_events():
    b = generate_paired_batch(32, cue_config=CueConfig(orient=True), seed=2, orient_noise=0.0)
    c = b.cues.view(32,2,*b.cues.shape[1:])
    assert (c[:,0,:,0] != c[:,1,:,0]).any()
