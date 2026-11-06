"""Tests for the accelerator memory release that unblocks ``foldseek --gpu 1``.

``baktfold run`` predicts 3Di with ProstT5 and then shells out to Foldseek's
GPU mode **from the same process**. PyTorch's caching allocator keeps the
model's blocks reserved against the device after the model itself has gone,
so without an explicit release Foldseek starts on a card that still looks full
and dies with::

    CUDA error: out of memory : .../libmarv/src/cudasw4.cuh, line 1276
    Error: Ungapped prefilter matching step died

These tests pin two things:
  1. ``release_gpu_memory`` empties the right allocator for the accelerator in
     use, and is a harmless no-op on a CPU-only machine,
  2. ``get_embeddings`` actually calls it, *after* inference returns (i.e.
     once the model references are unreachable), and skips both the model load
     and the release when there is nothing to predict.

No GPU, no model and no network needed — torch is stubbed out.

Run::

    pytest tests/unit/test_gpu_release.py -v
"""
from __future__ import annotations

from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
import pytest

from baktfold.features import predict_3Di
from baktfold.utils import gpu


# ─── fake torch ─────────────────────────────────────────────────────────────


def _fake_torch(*, cuda=False, mps=False, xpu=None, reserved=(8.0, 0.5)):
    """Stand-in for the torch module, recording which empty_cache() was hit.

    ``xpu=None`` models a torch build compiled without XPU support, where
    ``torch.xpu`` does not exist at all.
    """
    calls: list[str] = []
    reserved_seq = iter([int(r * 1024 ** 3) for r in reserved])

    fake = SimpleNamespace()
    fake.cuda = SimpleNamespace(
        is_available=lambda: cuda,
        empty_cache=lambda: calls.append("cuda"),
        memory_reserved=lambda: next(reserved_seq),
    )
    fake.backends = SimpleNamespace(mps=SimpleNamespace(is_available=lambda: mps))
    fake.mps = SimpleNamespace(empty_cache=lambda: calls.append("mps"))
    if xpu is not None:
        fake.xpu = SimpleNamespace(
            is_available=lambda: xpu,
            empty_cache=lambda: calls.append("xpu"),
        )
    return fake, calls


# ─── release_gpu_memory: per-accelerator dispatch ───────────────────────────


def test_release_empties_cuda_cache():
    """The bug-1 path: CUDA present -> hand the reserved blocks back."""
    fake, calls = _fake_torch(cuda=True)
    with patch.object(gpu, "torch", fake):
        gpu.release_gpu_memory("after ProstT5")
    assert calls == ["cuda"]


def test_release_prefers_cuda_over_other_backends():
    """A CUDA box should never fall through to the MPS/XPU branches."""
    fake, calls = _fake_torch(cuda=True, mps=True, xpu=True)
    with patch.object(gpu, "torch", fake):
        gpu.release_gpu_memory()
    assert calls == ["cuda"]


def test_release_empties_mps_cache():
    fake, calls = _fake_torch(cuda=False, mps=True)
    with patch.object(gpu, "torch", fake):
        gpu.release_gpu_memory()
    assert calls == ["mps"]


def test_release_empties_xpu_cache():
    fake, calls = _fake_torch(cuda=False, mps=False, xpu=True)
    with patch.object(gpu, "torch", fake):
        gpu.release_gpu_memory()
    assert calls == ["xpu"]


# ─── release_gpu_memory: CPU-only safety ────────────────────────────────────


def test_release_is_a_noop_on_cpu_only_builds():
    """No accelerator and no torch.xpu attribute — must not raise or free."""
    fake, calls = _fake_torch(cuda=False, mps=False, xpu=None)
    with patch.object(gpu, "torch", fake):
        gpu.release_gpu_memory("after ProstT5")
    assert calls == []


def test_release_is_a_noop_when_every_accelerator_reports_unavailable():
    fake, calls = _fake_torch(cuda=False, mps=False, xpu=False)
    with patch.object(gpu, "torch", fake):
        gpu.release_gpu_memory()
    assert calls == []


def test_release_collects_garbage_first():
    """nn.Module reference cycles survive refcounting, so gc runs before the
    allocator is emptied — otherwise there is nothing for it to release."""
    order: list[str] = []
    fake, _ = _fake_torch(cuda=True)
    fake.cuda.empty_cache = lambda: order.append("empty_cache")
    with patch.object(gpu, "torch", fake), \
         patch.object(gpu.gc, "collect", lambda: order.append("gc")):
        gpu.release_gpu_memory()
    assert order == ["gc", "empty_cache"]


def test_release_survives_memory_reserved_failure():
    """memory_reserved() is only for the log line — a failure there must not
    stop the actual release."""
    fake, calls = _fake_torch(cuda=True)

    def _boom():
        raise RuntimeError("no CUDA context")

    fake.cuda.memory_reserved = _boom
    with patch.object(gpu, "torch", fake):
        gpu.release_gpu_memory()
    assert calls == ["cuda"]


# ─── get_embeddings wiring ──────────────────────────────────────────────────


def _pred(indices, probs):
    """A (pred, mean_prob, all_prob) tuple shaped like the inference engine's."""
    return (
        np.asarray(indices, dtype=np.byte),
        float(np.mean(probs)) * 100,
        np.asarray([probs], dtype=np.float32),
    )


def _call_get_embeddings(tmp_path, cds_dict, hypotheticals, inference, order):
    """Drive get_embeddings with inference + release stubbed out."""
    with patch.object(predict_3Di, "run_prostt5_inference_multi_gpu", inference), \
         patch.object(predict_3Di, "release_gpu_memory",
                      lambda context="": order.append("release")), \
         patch.object(predict_3Di, "parse_gpus", lambda cpu, gpus: ["cuda:0"]):
        return predict_3Di.get_embeddings(
            hypotheticals,
            cds_dict,
            tmp_path,
            "genome",
            tmp_path / "model_dir",
            "Rostlab/ProstT5_fp16",
            tmp_path / "cnn.pt",
            tmp_path / "genome_3di.fasta",
            tmp_path / "genome_per_residue.h5",
            tmp_path / "genome_per_protein.h5",
            half_precision=True,
        )


def test_get_embeddings_releases_gpu_memory_after_inference(tmp_path):
    """Bug 1: the release must happen once inference has returned — that is
    the point at which the model, tokenizer and CNN head are unreachable."""
    order: list[str] = []

    def fake_inference(seq_dict, **kwargs):
        order.append("inference")
        return ({"cds_1": _pred([0, 1, 2], [0.9, 0.9, 0.9])}, {}, {}, [])

    feats = [{"locus": "cds_1", "id": "cds_1"}]
    preds = _call_get_embeddings(
        tmp_path, {"cds_1": "MKV"}, feats, fake_inference, order
    )

    assert order == ["inference", "release"]
    assert preds["cds_1"][0].tolist() == [0, 1, 2]
    assert (tmp_path / "genome_3di.fasta").read_text() == ">cds_1\nACD\n"


def test_get_embeddings_skips_model_load_when_nothing_to_predict(tmp_path):
    """Bug 2's genomes reach predict with zero hypotheticals. Loading ProstT5
    to predict nothing wastes a model load and several GB of VRAM — but every
    downstream file must still be written (empty)."""
    order: list[str] = []

    def fake_inference(seq_dict, **kwargs):  # pragma: no cover - must not run
        raise AssertionError("ProstT5 was loaded for an empty input")

    preds = _call_get_embeddings(tmp_path, {}, [], fake_inference, order)

    assert preds == {}
    assert order == []  # nothing allocated, so nothing to release
    assert (tmp_path / "genome_3di.fasta").read_text() == ""
    assert (tmp_path / "genome_prostT5_3di_mean_probabilities.csv").read_text() == ""
    assert (tmp_path / "genome_prostT5_3di_all_probabilities.json").read_text() == ""
    assert not (tmp_path / "fails.tsv").exists()
