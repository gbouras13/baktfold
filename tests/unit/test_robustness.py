"""Robustness fixes from the 2026-11-03 code audit (C1-C9).

Each test pins a failure mode that produced a crash, a silent wrong result,
or a misleading message on realistic input.

Run::

    pytest tests/unit/test_robustness.py -v
"""
from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path
from unittest.mock import patch

import pytest

from baktfold.features import run_foldseek
from baktfold.io import fasta
from baktfold.results.tophit import get_tophit
from baktfold.utils.validation import instantiate_dirs


def _baktfold():
    exe = Path(sys.executable).with_name("baktfold")
    exe = str(exe) if exe.exists() else shutil.which("baktfold")
    if exe is None:
        pytest.skip("baktfold console script not installed")
    return exe


# ── C1: duplicate FASTA IDs ──────────────────────────────────────────────────

def test_duplicate_fasta_ids_rejected(tmp_path):
    f = tmp_path / "p.faa"
    f.write_text(">p1\nMKV\n>p2\nMKL\n>p1\nMAA\n")
    with pytest.raises(ValueError, match="Duplicate sequence IDs.*p1"):
        fasta.import_sequences(f, False, False)


def test_unique_fasta_ids_accepted(tmp_path):
    f = tmp_path / "p.faa"
    f.write_text(">p1\nMKV\n>p2\nMKL\n")
    assert [s["id"] for s in fasta.import_sequences(f, False, False)] == ["p1", "p2"]


# ── C2: every ProstT5 prediction failed ─────────────────────────────────────

def _predict(tmp_path, predictions):
    from baktfold.subcommands import predict as predict_mod
    cds = {"c1": "MKVL", "c2": "MKLA"}
    feats = [{"locus": k, "id": k, "aa": v} for k, v in cds.items()]
    with patch.object(predict_mod, "get_embeddings", lambda *a, **k: predictions):
        return predict_mod.subcommand_predict(
            feats, cds, tmp_path, "p", cpu=True, omit_probs=False,
            model_dir=tmp_path, model_name="m", checkpoint_path=tmp_path / "c",
            batch_size=1, mask_threshold=25, save_per_residue_embeddings=False,
            save_per_protein_embeddings=False, threads=1, has_duplicate_locus=False,
        )


def test_all_predictions_failed_exits(tmp_path):
    with pytest.raises(SystemExit) as exc:
        _predict(tmp_path, {})
    assert exc.value.code == 1


def test_only_unusable_prediction_does_not_abort(tmp_path):
    """Audit N5: inference ran, but its only output is a single-residue
    (0-d) prediction. That protein is skipped as before (#47); the run must
    not be treated as a total failure."""
    import numpy as np
    pred = {"c1": (np.asarray(3, dtype=np.byte), 90.0, np.ones((1, 1)))}
    _predict(tmp_path, pred)
    assert (tmp_path / "p_aa.fasta").read_text() == ">c1\nXXXX\n>c2\nXXXX\n"


def test_partial_failure_continues(tmp_path):
    import numpy as np
    pred = {"c1": (np.asarray([0, 1, 2, 3], dtype=np.byte), 90.0, np.ones((1, 4)))}
    _predict(tmp_path, pred)
    assert (tmp_path / "p_aa.fasta").read_text() == ">c1\nMKVL\n>c2\nXXXX\n"


# ── C3: custom database with --foldseek-gpu ─────────────────────────────────

def test_gpu_target_uses_gpu_copy_when_present(tmp_path):
    (tmp_path / "swissprot_gpu.dbtype").write_text("")
    assert run_foldseek.gpu_target_db(tmp_path / "swissprot") == f"{tmp_path / 'swissprot'}_gpu"


def test_gpu_target_uses_given_padded_db(tmp_path):
    """docs: --custom-db custom_db/db_pad (already padded, no '_gpu' copy)."""
    (tmp_path / "db_pad.dbtype").write_text("")
    assert run_foldseek.gpu_target_db(tmp_path / "db_pad") == str(tmp_path / "db_pad")


# ── C4: --gpus under an inherited CUDA_VISIBLE_DEVICES (SLURM) ──────────────

@pytest.mark.parametrize(
    "inherited, torch_idx, expected",
    [
        ("2,3", "1", "3"),
        ("2,3", "0,1", "2,3"),
        ("", "1", "1"),          # nothing inherited: indices are physical
        ("GPU-abc,GPU-def", "1", "GPU-def"),  # UUID form
    ],
)
def test_physical_cuda_visible_devices(monkeypatch, inherited, torch_idx, expected):
    monkeypatch.setenv("CUDA_VISIBLE_DEVICES", inherited)
    assert run_foldseek.physical_cuda_visible_devices(torch_idx) == expected


def test_physical_cuda_index_out_of_range(monkeypatch):
    monkeypatch.setenv("CUDA_VISIBLE_DEVICES", "2")
    with pytest.raises(ValueError, match="outside CUDA_VISIBLE_DEVICES"):
        run_foldseek.physical_cuda_visible_devices("1")


# ── C6: -f must not delete the inputs ───────────────────────────────────────

def test_force_refuses_when_input_inside_output(tmp_path):
    out = tmp_path / "out"
    out.mkdir()
    inp = out / "baktfold.json"
    inp.write_text("{}")
    with pytest.raises(SystemExit):
        instantiate_dirs(out, force=True, protected=[inp])
    assert inp.exists()


def test_force_refuses_when_output_is_predictions_dir(tmp_path):
    pred = tmp_path / "pred"
    pred.mkdir()
    (pred / "x_3di.fasta").write_text(">a\nA\n")
    with pytest.raises(SystemExit):
        instantiate_dirs(pred, force=True, protected=[tmp_path / "in.json", pred])
    assert (pred / "x_3di.fasta").exists()


def test_force_still_works_for_unrelated_inputs(tmp_path):
    out = tmp_path / "out"
    out.mkdir()
    (out / "stale").write_text("x")
    instantiate_dirs(out, force=True, protected=[tmp_path / "in.json", None])
    assert out.exists() and not (out / "stale").exists()


def test_cli_json_force_does_not_delete_its_input(tmp_path):
    out = tmp_path / "out"
    out.mkdir()
    src = Path(__file__).resolve().parents[1] / "test_data" / "golden" / "run_all" / "baktfold.json"
    inp = out / "baktfold.json"
    shutil.copy(src, inp)
    proc = subprocess.run([_baktfold(), "json", "-i", str(inp), "-o", str(out), "-f"],
                          capture_output=True, text=True)
    assert proc.returncode != 0
    assert inp.exists()


# ── C7: numeric IDs ─────────────────────────────────────────────────────────

def test_numeric_query_and_target_ids(tmp_path):
    t = tmp_path / "hits.tsv"
    t.write_text("1\t2\t100\t0.5\t1e-10\t1\t50\t50\t1\t50\t50\n")
    df = get_tophit(t, structures=False, cath=False)
    assert df.select("query", "target").rows() == [("1", "2")]


# ── C9: missing / empty input message ───────────────────────────────────────

def test_missing_input_reports_the_real_error(tmp_path):
    proc = subprocess.run([_baktfold(), "json", "-i", str(tmp_path / "nope.json"), "-o", str(tmp_path / "o")],
                          capture_output=True, text=True)
    text = proc.stdout + proc.stderr
    assert proc.returncode != 0
    assert "not readable" in text
    assert "is not defined" not in text
