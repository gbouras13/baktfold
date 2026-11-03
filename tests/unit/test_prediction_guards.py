"""Tests for the guards that keep one bad CDS from killing a whole genome.

``baktfold predict`` used to filter unusable ProstT5 predictions with a bare
``len(v[0]) > 0``. That handles a zero-length array (issue #47) but **raises**
on a zero-dimensional one::

    TypeError: len() of unsized object

numpy produces a 0-d array whenever ``.squeeze()`` in pholdlib's inference loop
collapses a single-residue prediction to a scalar — ``prediction[i, :, 0:1]``
has shape ``(1, 1)``, and squeezing that gives shape ``()``. pholdlib noticed
(its own ``len()`` is wrapped in try/except) and recorded the id for fails.tsv
— but left the entry in the predictions dict, so it reached baktfold, where
the guard itself threw and took the entire genome's annotation down with it.

pholdlib now drops those entries at source, but baktfold accepts
``pholdlib >=0.1.2`` and the released 0.1.2 still returns them, so this guard
is load-bearing rather than belt-and-braces. A partial failure has to degrade
gracefully: drop the bad CDS, annotate the rest, and say what was dropped.

Run::

    pytest tests/unit/test_prediction_guards.py -v
"""
from __future__ import annotations

from unittest.mock import patch

import numpy as np
import pytest
from loguru import logger

from baktfold.features.predict_3Di import drop_unusable_predictions, write_predictions
from baktfold.subcommands import predict as predict_mod


# ─── prediction shapes the inference engine really emits ────────────────────


def _good(indices=(0, 1, 2)):
    """A normal multi-residue prediction."""
    probs = [0.9] * len(indices)
    return (np.asarray(indices, dtype=np.byte), 90.0, np.asarray([probs], dtype=np.float32))


def _unsized():
    """A single-residue prediction after .squeeze() — shape (), len() raises."""
    pred = np.zeros((1, 1), dtype=np.byte)[:, 0:1].squeeze()
    assert pred.ndim == 0, "fixture must reproduce the 0-d shape"
    with pytest.raises(TypeError):
        len(pred)
    return (pred, 90.0, np.asarray([[0.9]], dtype=np.float32))


def _zero_length():
    """The issue #47 shape — a real array with no residues."""
    return (np.zeros((0,), dtype=np.byte), 0.0, np.zeros((1, 0), dtype=np.float32))


# ─── drop_unusable_predictions ──────────────────────────────────────────────


def test_drops_unsized_prediction():
    kept = drop_unusable_predictions({"cds_1": _unsized()})
    assert kept == {}


def test_drops_zero_length_prediction():
    """Issue #47's case must keep behaving exactly as it did."""
    kept = drop_unusable_predictions({"cds_1": _zero_length()})
    assert kept == {}


def test_drops_none_prediction():
    kept = drop_unusable_predictions({"cds_1": (None, 0.0, None)})
    assert kept == {}


def test_keeps_usable_predictions_and_drops_only_the_bad_ones():
    """The whole point: one bad CDS must not cost the genome its annotation."""
    preds = {
        "cds_good_1": _good(),
        "cds_unsized": _unsized(),
        "cds_empty": _zero_length(),
        "cds_good_2": _good([4, 5]),
    }
    kept = drop_unusable_predictions(preds)
    assert sorted(kept) == ["cds_good_1", "cds_good_2"]


def test_no_warning_when_everything_is_usable():
    messages: list[str] = []
    sink = logger.add(lambda m: messages.append(str(m)), level="WARNING")
    try:
        drop_unusable_predictions({"cds_1": _good()})
    finally:
        logger.remove(sink)
    assert messages == []


def test_reports_what_was_dropped():
    """A silent drop plus a hard crash gave users nothing to act on."""
    messages: list[str] = []
    sink = logger.add(lambda m: messages.append(str(m)), level="WARNING")
    try:
        drop_unusable_predictions({"cds_bad": _unsized()}, "writing 3Di FASTA")
    finally:
        logger.remove(sink)

    assert len(messages) == 1
    assert "Dropping 1 CDS" in messages[0]
    assert "cds_bad" in messages[0]
    assert "writing 3Di FASTA" in messages[0]


def test_dropped_id_list_is_truncated():
    """Do not dump thousands of ids into the log on a catastrophic run."""
    messages: list[str] = []
    sink = logger.add(lambda m: messages.append(str(m)), level="WARNING")
    try:
        drop_unusable_predictions({f"cds_{i}": _zero_length() for i in range(25)})
    finally:
        logger.remove(sink)

    assert "Dropping 25 CDS" in messages[0]
    assert "and 15 more" in messages[0]


# ─── write_predictions: the crash that killed whole genomes ─────────────────


def test_write_predictions_survives_an_unsized_entry(tmp_path):
    """Before the fix this raised TypeError and aborted baktfold predict."""
    feats = [
        {"locus": "cds_good", "id": "cds_good"},
        {"locus": "cds_bad", "id": "cds_bad"},
    ]
    preds = {"cds_good": _good(), "cds_bad": _unsized()}
    out = tmp_path / "3di.fasta"

    write_predictions(feats, preds, out, mask_threshold=0)

    assert out.read_text() == ">cds_good\nACD\n"
    assert feats[0]["3di"] == "ACD"
    assert feats[1]["3di"] is None  # dropped, but the genome still annotates


def test_write_predictions_still_masks_normally(tmp_path):
    """The guard must not disturb ordinary confidence masking."""
    feats = [{"locus": "cds_1", "id": "cds_1"}]
    preds = {"cds_1": (
        np.asarray([0, 1, 2], dtype=np.byte), 60.0,
        np.asarray([[0.9, 0.1, 0.9]], dtype=np.float32),
    )}

    write_predictions(feats, preds, tmp_path / "3di.fasta", mask_threshold=50)

    assert feats[0]["3di"] == "AXD"  # middle residue masked


# ─── the production case ────────────────────────────────────────────────────


def test_one_single_residue_cds_does_not_cost_the_genome_its_annotation():
    """Eimeria tenella (ETH001): exactly ONE 1-residue CDS, locus
    TG0AAQ_3971_3973, aborted annotation for all 4,157 of the genome's
    hypothetical proteins. Other real triggers seen in the same rerun:
    TRDBJZ_150_155, JWHKG7_62_67 — every predict-failed protist genome held at
    least one 1-residue CDS, and the run logged 29 "prediction has length 0"
    warnings in total.

    Scaled down here, but the shape of the failure is the point: a single bad
    CDS must cost exactly itself.
    """
    preds = {f"ETH_{i:04d}": _good() for i in range(50)}
    preds["TG0AAQ_3971_3973"] = _unsized()

    kept = drop_unusable_predictions(preds, "writing 3Di FASTA")

    assert "TG0AAQ_3971_3973" not in kept
    assert len(kept) == 50  # every other CDS still annotates


# ─── subcommand_predict uses the same guard ─────────────────────────────────


def test_subcommand_predict_survives_an_unsized_entry(tmp_path):
    """The second `len(v[0]) > 0` call site — it crashed the same way."""
    cds_dict = {"cds_good": "MKV", "cds_bad": "M"}
    feats = [
        {"locus": "cds_good", "id": "cds_good"},
        {"locus": "cds_bad", "id": "cds_bad"},
    ]

    def fake_get_embeddings(hypotheticals, cds, *args, **kwargs):
        return {"cds_good": _good(), "cds_bad": _unsized()}

    with patch.object(predict_mod, "get_embeddings", fake_get_embeddings):
        predict_mod.subcommand_predict(
            feats, cds_dict, tmp_path, "genome",
            cpu=True, omit_probs=False,
            model_dir=tmp_path / "db", model_name="Rostlab/ProstT5_fp16",
            checkpoint_path=tmp_path / "cnn.pt", batch_size=1,
            save_per_residue_embeddings=False, save_per_protein_embeddings=False,
            threads=1, mask_threshold=0, has_duplicate_locus=False,
        )

    aa = (tmp_path / "genome_aa.fasta").read_text()
    # the good CDS keeps its sequence; the dropped one is fully masked
    assert ">cds_good\nMKV\n" in aa
    assert ">cds_bad\nX\n" in aa
