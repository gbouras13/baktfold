"""
Function-level tests for baktfold.features.predict_3Di write helpers.

``write_predictions`` (3Di FASTA + in-place Bakta-feature update) and
``write_embeddings`` (HDF5) are pure given a pre-built predictions dict — no
torch / GPU needed at call time (importing the module does pull torch via
pholdlib; see conftest note on #29).

Guards:
  - REVIEW_FINDINGS #26 (per-residue masking loop -> vectorise): the masking
    snapshot pins exactly which residues become 'X'.
  - REVIEW_FINDINGS #16 (h5py atomic write): the round-trip test gates the rewrite.
  - REVIEW_FINDINGS #38 (w+ -> w): unaffected, write content pinned here.

Run::

    pytest tests/unit/test_predict_3Di.py -v
"""
from __future__ import annotations

import h5py
import numpy as np
import pytest

from baktfold.features.predict_3Di import write_predictions, write_embeddings
from pholdlib.prostt5.output import SS_MAPPING


def _pred(indices, probs):
    """Build a (pred, mean_prob, all_prob) tuple like the inference engine emits.

    pred     : int8 array of 3Di class indices (0..20)
    all_prob : float32 (1, L) max-softmax confidence per residue
    """
    pred = np.asarray(indices, dtype=np.byte)
    all_prob = np.asarray([probs], dtype=np.float32)  # shape (1, L)
    mean_prob = float(np.mean(probs)) if len(probs) else 0.0
    return (pred, mean_prob, all_prob)


# ─── write_predictions: basic ───────────────────────────────────────────────


def test_write_predictions_basic(tmp_path):
    """No masking (threshold 0): 3Di string is SS_MAPPING of the indices,
    feature dict gets the 3Di string, FASTA has one record per prediction."""
    feats = [{"locus": "cds_1", "id": "cds_1"}, {"locus": "cds_2", "id": "cds_2"}]
    preds = {
        "cds_1": _pred([0, 1, 2, 3], [0.9, 0.9, 0.9, 0.9]),   # ACDE
        "cds_2": _pred([5, 6, 7], [0.9, 0.9, 0.9]),           # GHI
    }
    out = tmp_path / "out_3di.fasta"
    write_predictions(feats, preds, out, mask_threshold=0)

    assert feats[0]["3di"] == "ACDE"
    assert feats[1]["3di"] == "GHI"
    assert out.read_text() == ">cds_1\nACDE\n>cds_2\nGHI\n"


# ─── write_predictions: confidence masking (#26 gate) ───────────────────────


def test_write_predictions_masks_low_confidence(tmp_path):
    """mask_threshold=50 -> mask_prop=0.5. Residue with all_prob 0.3 (<0.5)
    has its index forced to 20 ('X')."""
    feats = [{"locus": "cds_1", "id": "cds_1"}]
    preds = {"cds_1": _pred([0, 1, 2, 3], [0.9, 0.3, 0.9, 0.9])}  # ACDE, idx1 low
    out = tmp_path / "out_3di.fasta"
    write_predictions(feats, preds, out, mask_threshold=50)

    assert feats[0]["3di"] == "AXDE"  # position 1 masked
    assert SS_MAPPING[20] == "X"


def test_write_predictions_no_mask_at_threshold_zero(tmp_path):
    """threshold 0 -> mask_prop 0 -> nothing masked even at low probs."""
    feats = [{"locus": "cds_1", "id": "cds_1"}]
    preds = {"cds_1": _pred([0, 1, 2], [0.0, 0.0, 0.0])}
    write_predictions(feats, preds, tmp_path / "o.fasta", mask_threshold=0)
    assert feats[0]["3di"] == "ACD"  # all_prob 0.0 is not < 0.0


# ─── write_predictions: id vs locus keying ──────────────────────────────────


def test_write_predictions_has_duplicate_locus_uses_id(tmp_path):
    feats = [{"locus": "shared", "id": "uniq_1"}]
    preds = {"uniq_1": _pred([0, 1], [0.9, 0.9])}
    out = tmp_path / "o.fasta"
    write_predictions(feats, preds, out, mask_threshold=0, has_duplicate_locus=True)
    assert feats[0]["3di"] == "AC"
    assert out.read_text() == ">uniq_1\nAC\n"


# ─── write_predictions: missing / zero-length predictions ───────────────────


def test_write_predictions_missing_prediction_sets_none(tmp_path):
    """A feature with no prediction (OOM / corrupt) -> feat['3di']=None, no FASTA line."""
    feats = [{"locus": "cds_1", "id": "cds_1"}, {"locus": "cds_2", "id": "cds_2"}]
    preds = {"cds_1": _pred([0, 1], [0.9, 0.9])}  # cds_2 absent
    out = tmp_path / "o.fasta"
    write_predictions(feats, preds, out, mask_threshold=0)
    assert feats[1]["3di"] is None
    assert out.read_text() == ">cds_1\nAC\n"  # only cds_1 written


def test_write_predictions_zero_length_dropped(tmp_path):
    """Zero-length prediction (issue #47) is filtered -> treated as missing."""
    feats = [{"locus": "cds_1", "id": "cds_1"}]
    preds = {"cds_1": _pred([], [])}
    out = tmp_path / "o.fasta"
    write_predictions(feats, preds, out, mask_threshold=0)
    assert feats[0]["3di"] is None
    assert out.read_text() == ""


# ─── write_embeddings: h5py round-trip (#16 gate) ───────────────────────────


def test_write_embeddings_roundtrip(tmp_path):
    emb = {
        "cds_1": np.arange(4, dtype=np.float32),
        "cds_2": np.ones((2, 3), dtype=np.float32),
    }
    out = tmp_path / "emb.h5"
    write_embeddings(emb, out)

    assert h5py.is_hdf5(str(out))
    with h5py.File(str(out), "r") as hf:
        assert set(hf.keys()) == {"cds_1", "cds_2"}
        np.testing.assert_array_equal(hf["cds_1"][:], emb["cds_1"])
        np.testing.assert_array_equal(hf["cds_2"][:], emb["cds_2"])
