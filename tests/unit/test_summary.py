"""
Function-level tests for baktfold.io.io.write_summary_txt_file.

Regression for the proteins-compare summary bug: the counts were derived from
``str(feat).lower()`` substring matches, so any feature whose retained ``pstc``
list held a secondary hit described "hypothetical protein" (e.g. an AFDB entry)
was miscounted as still-hypothetical even after a real product was assigned —
reporting "CDS annotated with Baktfold function: 0" and
"CDS remaining hypotheticals: <all>".

The counts now read the authoritative feature keys set by
bakta.annotation: ``feature['baktfold']`` (any PSTC hit) and
``feature['hypothetical']`` (present iff still hypothetical).

Run::

    pytest tests/unit/test_summary.py -v
"""
from __future__ import annotations

from pathlib import Path

import baktfold.bakta.constants as bc
from baktfold.io.io import write_summary_txt_file


def _summary_counts(text: str) -> dict:
    counts = {}
    for line in text.splitlines():
        if line.startswith("CDS "):
            label, _, val = line.rpartition(":")
            counts[label.strip()] = int(val.strip())
    return counts


def _annotated(i: int) -> dict:
    """A CDS that got a real Baktfold function (hypothetical key removed),
    but whose retained pstc list still mentions 'hypothetical protein'."""
    return {
        "type": bc.FEATURE_CDS,
        "id": f"cds_{i}",
        "baktfold": True,
        "product": "DNA polymerase",
        # secondary structural hit described as hypothetical protein — this is
        # what poisoned the old str(feat) substring count
        "pstc": [
            {"source": "swissprot", "id": f"P{i:05d}", "description": "DNA polymerase"},
            {"source": "afdb", "id": f"A{i:05d}", "description": "hypothetical protein"},
        ],
    }


def _remaining_hypothetical(i: int) -> dict:
    """A CDS with no usable hit — still hypothetical."""
    return {
        "type": bc.FEATURE_CDS,
        "id": f"cds_{i}",
        "hypothetical": True,
        "product": bc.HYPOTHETICAL_PROTEIN,
    }


def test_summary_counts_real_function_despite_secondary_hypothetical(tmp_path):
    feats = [_annotated(i) for i in range(3)] + [_remaining_hypothetical(i) for i in range(3, 5)]

    write_summary_txt_file(str(tmp_path), "baktfold", feats)
    counts = _summary_counts((tmp_path / "baktfold.summary.txt").read_text())

    assert counts["CDS count"] == 5
    assert counts["CDS annotated with Baktfold database hit"] == 3
    assert counts["CDS annotated with Baktfold function"] == 3  # was 0 pre-fix
    assert counts["CDS remaining hypotheticals"] == 2           # was 5 pre-fix
    assert counts["CDS beginning hypotheticals"] == 5
    # internal consistency: function + remaining == beginning
    assert (counts["CDS annotated with Baktfold function"]
            + counts["CDS remaining hypotheticals"]
            == counts["CDS beginning hypotheticals"])


def test_summary_hit_to_hypothetical_protein_is_not_a_function(tmp_path):
    """A PSTC hit whose only description is 'hypothetical protein' counts as a
    database hit but NOT a function, and stays a remaining hypothetical."""
    feat = {
        "type": bc.FEATURE_CDS,
        "id": "cds_0",
        "baktfold": True,
        "hypothetical": True,  # re-marked because product resolved to hypothetical protein
        "product": bc.HYPOTHETICAL_PROTEIN,
        "pstc": [{"source": "afdb", "id": "A1", "description": "hypothetical protein"}],
    }

    write_summary_txt_file(str(tmp_path), "baktfold", [feat])
    counts = _summary_counts((tmp_path / "baktfold.summary.txt").read_text())

    assert counts["CDS annotated with Baktfold database hit"] == 1
    assert counts["CDS annotated with Baktfold function"] == 0
    assert counts["CDS remaining hypotheticals"] == 1


def test_summary_non_cds_features_ignored(tmp_path):
    feats = [
        _annotated(0),
        {"type": "tRNA", "id": "t1", "hypothetical": True},  # not a CDS — must not count
    ]

    write_summary_txt_file(str(tmp_path), "baktfold", feats)
    counts = _summary_counts((tmp_path / "baktfold.summary.txt").read_text())

    assert counts["CDS count"] == 1
    assert counts["CDS remaining hypotheticals"] == 0
