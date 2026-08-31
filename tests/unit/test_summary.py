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

import re
from pathlib import Path

import baktfold.bakta.constants as bc
from baktfold.io.io import write_summary_txt_file

# "CDS annotated with Baktfold function: 7 (0.3% of CDS; 12.7% of beginning hypotheticals)"
#  \_______________ label _____________/  \_/  \_____________ percentages ______________/
_CDS_LINE = re.compile(r"^(CDS [^:]*):\s*(\d+)(?:\s*\((.*)\))?$")


def _summary_lines(text: str) -> dict:
    """Map each ``CDS ...`` label to its (count, percentage-text) pair.

    ``percentage-text`` is the content of the trailing parentheses, or ``""``
    for a line that carries no percentages (the ``CDS count`` denominator).
    """
    lines = {}
    for line in text.splitlines():
        m = _CDS_LINE.match(line)
        if m:
            lines[m.group(1).strip()] = (int(m.group(2)), m.group(3) or "")
    return lines


def _summary_counts(text: str) -> dict:
    return {label: count for label, (count, _) in _summary_lines(text).items()}


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


def _never_hypothetical(i: int) -> dict:
    """A CDS bakta already annotated — baktfold is never given it, so it counts
    towards the CDS total but towards no baktfold denominator."""
    return {
        "type": bc.FEATURE_CDS,
        "id": f"cds_{i}",
        "product": "recombinase RecA",
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


def test_summary_percentages_use_both_denominators(tmp_path):
    """Each baktfold count is reported as a share of all CDS *and* of the
    beginning hypotheticals — the set baktfold is actually handed."""
    feats = (
        [_annotated(i) for i in range(3)]
        + [_remaining_hypothetical(i) for i in range(3, 5)]
        + [_never_hypothetical(i) for i in range(5, 10)]
    )

    write_summary_txt_file(str(tmp_path), "baktfold", feats)
    lines = _summary_lines((tmp_path / "baktfold.summary.txt").read_text())

    # 10 CDS, 5 of them beginning hypotheticals (3 annotated + 2 remaining)
    assert lines["CDS count"] == (10, "")  # the denominator itself carries no %
    assert lines["CDS beginning hypotheticals"] == (5, "50.0% of CDS")
    assert lines["CDS annotated with Baktfold database hit"] == (
        3, "30.0% of CDS; 60.0% of beginning hypotheticals")
    assert lines["CDS annotated with Baktfold function"] == (
        3, "30.0% of CDS; 60.0% of beginning hypotheticals")
    assert lines["CDS remaining hypotheticals"] == (
        2, "20.0% of CDS; 40.0% of beginning hypotheticals")


def test_summary_percentages_with_no_cds(tmp_path):
    """No CDS (or no hypotheticals) must not divide by zero, and must still
    format as '0.0' so the golden/snapshot outputs stay stable."""
    write_summary_txt_file(str(tmp_path), "baktfold", [{"type": "tRNA", "id": "t1"}])
    lines = _summary_lines((tmp_path / "baktfold.summary.txt").read_text())

    assert lines["CDS count"] == (0, "")
    assert lines["CDS beginning hypotheticals"] == (0, "0.0% of CDS")
    assert lines["CDS annotated with Baktfold database hit"] == (
        0, "0.0% of CDS; 0.0% of beginning hypotheticals")
