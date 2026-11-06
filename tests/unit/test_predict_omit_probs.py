"""
Regression test for ``--omit-probs`` on the 3Di predict path.

``--omit-probs`` used to be forwarded to pholdlib as ``output_probs=False``.
pholdlib then returns ``all_prob=None``, and 3Di masking crashed with
``TypeError: 'NoneType' object is not subscriptable`` in ``predict``,
``proteins-predict``, ``run`` and ``proteins``. The flag must only drop the
per-residue probability JSON. Masking and the mean probability CSV are
unaffected.

This drives the real path (``subcommand_predict`` -> ``get_embeddings`` ->
pholdlib's dispatcher and inference loop -> writers) on CPU, with ProstT5 and
the CNN head swapped for tiny stand-ins, so it needs no model download,
database or GPU.

Run::

    pytest tests/unit/test_predict_omit_probs.py -v
"""
from __future__ import annotations

from types import SimpleNamespace

import pytest
import torch
from pholdlib.prostt5 import inference

from baktfold.subcommands.predict import subcommand_predict

PREFIX = "baktfold"
PROTEINS = {"cds_1": "MKVLA", "cds_2": "MSTNPKQR", "cds_3": "MAG"}


class _Encoding:
    """Stand-in for the tokenizer's BatchEncoding."""

    def __init__(self, attention_mask):
        self.input_ids = attention_mask.clone()
        self.attention_mask = attention_mask

    def to(self, device):
        return self


def _fake_tokenizer(seqs, **kwargs):
    """'<AA2fold> M K V' -> prefix + one token per residue + </s>, padded."""
    n_res = [len(s.split()) - 1 for s in seqs]
    attention_mask = torch.zeros(len(seqs), max(n_res) + 2, dtype=torch.long)
    for i, n in enumerate(n_res):
        attention_mask[i, : n + 2] = 1
    return _Encoding(attention_mask)


def _fake_encoder(input_ids, attention_mask=None):
    return SimpleNamespace(last_hidden_state=torch.zeros(*input_ids.shape, 8))


def _fake_cnn(residue_embedding):
    """Logits (B x 20 x L): residue j predicts 3Di class j % 20, confidently
    (p ~ 1.0) at even j and weakly (p ~ 0.08, masked at 25) at odd j."""
    n_batch, length = residue_embedding.shape[:2]
    logits = torch.zeros(n_batch, 20, length)
    for j in range(length):
        logits[:, j % 20, j] = 10.0 if j % 2 == 0 else 0.5
    return logits


@pytest.fixture
def fake_prostt5(monkeypatch):
    """Swap the model loaders inside pholdlib for the stand-ins above."""
    monkeypatch.setattr(
        inference,
        "get_T5_model",
        lambda *args, **kwargs: (_fake_encoder, _fake_tokenizer, torch.device("cpu")),
    )
    monkeypatch.setattr(inference, "load_predictor", lambda *args, **kwargs: _fake_cnn)


@pytest.mark.parametrize("batch_size", [1, 3])
@pytest.mark.parametrize("omit_probs", [False, True])
def test_omit_probs_only_drops_the_per_residue_json(
    tmp_path, fake_prostt5, omit_probs, batch_size
):
    feats = [{"locus": k, "id": k} for k in PROTEINS]
    subcommand_predict(
        feats,
        dict(PROTEINS),
        tmp_path,
        PREFIX,
        cpu=True,
        omit_probs=omit_probs,
        model_dir=tmp_path,
        model_name="unused",
        checkpoint_path=tmp_path / "unused.pt",
        batch_size=batch_size,
        save_per_residue_embeddings=False,
        save_per_protein_embeddings=False,
        threads=1,
        mask_threshold=25,
        has_duplicate_locus=False,
    )

    # low-confidence (odd) residues are masked in both 3Di and AA outputs
    assert [f["3di"] for f in feats] == ["AXDXF", "AXDXFXHX", "AXD"]
    assert (tmp_path / f"{PREFIX}_3di.fasta").read_text() == (
        ">cds_1\nAXDXF\n>cds_2\nAXDXFXHX\n>cds_3\nAXD\n"
    )
    assert (tmp_path / f"{PREFIX}_aa.fasta").read_text() == (
        ">cds_1\nMXVXA\n>cds_2\nMXTXPXQX\n>cds_3\nMXG\n"
    )

    mean_rows = [
        line.split(",")
        for line in (tmp_path / f"{PREFIX}_prostT5_3di_mean_probabilities.csv")
        .read_text()
        .splitlines()
    ]
    assert [row[0] for row in mean_rows] == list(PROTEINS)
    assert [float(row[1]) for row in mean_rows] == pytest.approx(
        [63.14, 53.95, 69.27], abs=0.01
    )

    all_probs_json = tmp_path / f"{PREFIX}_prostT5_3di_all_probabilities.json"
    assert all_probs_json.exists() is not omit_probs
