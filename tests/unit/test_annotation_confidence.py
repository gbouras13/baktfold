"""
Unit tests for baktfold.bakta.annotation.assign_annotation_confidence.

Mirrors phold's confidence heuristic (high/medium/low) on the
hierarchy-selected Foldseek hit. DB/GPU-free.
"""
from __future__ import annotations

import pytest

from baktfold.bakta.annotation import assign_annotation_confidence


def _entry(qcov, tcov, fident, evalue):
    return {"query_cov": qcov, "subject_cov": tcov, "identity": fident, "evalue": evalue}


# ── ProstT5 path (structures=False) ──────────────────────────────────────────


@pytest.mark.parametrize(
    "qcov,tcov,fident,evalue,prostt5,expected",
    [
        # high: reciprocal >0.8 coverage + one of (fident>0.3 / prostt5>60 / evalue<1e-10)
        (0.9, 0.9, 0.40, 1e-3, 10, "high"),    # via fident
        (0.9, 0.9, 0.10, 1e-3, 70, "high"),    # via prostt5 > 60
        (0.9, 0.9, 0.10, 1e-12, 10, "high"),   # via evalue < 1e-10
        # medium: one of qcov/tcov >0.8 + (fident>0.3 OR 45<=prostt5<=60) + evalue<1e-5
        (0.9, 0.5, 0.10, 1e-6, 50, "medium"),  # via prostt5 in [45,60]
        (0.9, 0.5, 0.40, 1e-6, 10, "medium"),  # via fident
        # low: fails the above
        (0.9, 0.5, 0.10, 1e-3, 50, "low"),     # evalue too weak for medium
        (0.5, 0.5, 0.40, 1e-12, 90, "low"),    # neither coverage > 0.8
        (0.9, 0.9, 0.10, 1e-3, 30, "low"),     # recip cov but nothing else for high; evalue weak for medium
    ],
)
def test_prostt5_path(qcov, tcov, fident, evalue, prostt5, expected):
    assert assign_annotation_confidence(_entry(qcov, tcov, fident, evalue), prostt5, structures=False) == expected


# ── structures path (no ProstT5 confidence) ─────────────────────────────────


@pytest.mark.parametrize(
    "qcov,tcov,fident,evalue,expected",
    [
        (0.9, 0.9, 0.40, 1e-3, "high"),     # recip cov + fident
        (0.9, 0.9, 0.10, 1e-12, "high"),    # recip cov + evalue
        (0.9, 0.5, 0.40, 1e-3, "medium"),   # one cov + fident
        (0.5, 0.5, 0.40, 1e-12, "low"),     # no cov
        (0.9, 0.9, 0.10, 1e-3, "low"),      # recip cov but no fident/evalue for high; medium needs fident or evalue<1e-5
    ],
)
def test_structures_path(qcov, tcov, fident, evalue, expected):
    assert assign_annotation_confidence(_entry(qcov, tcov, fident, evalue), None, structures=True) == expected


def test_prostt5_confidence_none_falls_back_to_structures_logic():
    # a ProstT5-mode protein whose prediction failed (no confidence) must not crash
    assert assign_annotation_confidence(_entry(0.9, 0.9, 0.4, 1e-3), None, structures=False) == "high"


def test_evalue_string_is_coerced():
    assert assign_annotation_confidence(
        {"query_cov": 0.9, "subject_cov": 0.9, "identity": 0.1, "evalue": "1e-12"}, 10, structures=False
    ) == "high"
