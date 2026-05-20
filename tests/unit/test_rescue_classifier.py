"""Unit tests for rescue classifier model."""

import pytest
import torch

from refold.models.rescue_classifier.model import RescueClassifier, FocalLoss
from refold.constants import (
    GNN_NODE_DIM, GNN_EDGE_DIM, ESM2_EMBED_DIM,
    THERMO_FEAT_DIM, EVO_FEAT_DIM,
)


def _make_batch(batch_size=2, n_nodes=10):
    node_feats = torch.randn(batch_size, n_nodes, GNN_NODE_DIM)
    edge_index = torch.zeros(2, 1, dtype=torch.long)
    edge_feats = torch.randn(1, GNN_EDGE_DIM)
    esm2 = torch.randn(batch_size, ESM2_EMBED_DIM)
    thermo = torch.randn(batch_size, THERMO_FEAT_DIM)
    evo = torch.randn(batch_size, EVO_FEAT_DIM)
    return node_feats, edge_index, edge_feats, esm2, thermo, evo


def test_rescue_classifier_forward_pass():
    model = RescueClassifier()
    model.eval()
    nf, ei, ef, esm2, thermo, evo = _make_batch(batch_size=2)
    with torch.no_grad():
        logits = model(nf, ei, ef, esm2, thermo, evo)
    assert logits.shape == (2, 1)


def test_rescue_probability_is_in_0_1():
    model = RescueClassifier()
    model.eval()
    nf, ei, ef, esm2, thermo, evo = _make_batch(batch_size=4)
    with torch.no_grad():
        probs = model.predict_proba(nf, ei, ef, esm2, thermo, evo)
    assert probs.shape == (4, 1)
    assert (probs >= 0).all() and (probs <= 1).all()


def test_focal_loss_computation():
    loss_fn = FocalLoss(alpha=0.25, gamma=2.0)
    logits = torch.randn(8, 1)
    targets = torch.randint(0, 2, (8, 1)).float()
    loss = loss_fn(logits, targets)
    assert loss.item() >= 0
    assert not torch.isnan(loss)


def test_rescue_classifier_single_residue():
    model = RescueClassifier()
    model.eval()
    nf, ei, ef, esm2, thermo, evo = _make_batch(batch_size=1, n_nodes=1)
    with torch.no_grad():
        logits = model(nf, ei, ef, esm2, thermo, evo)
    assert logits.shape == (1, 1)
