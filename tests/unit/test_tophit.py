"""
Function-level tests for baktfold.results.tophit.get_tophit.

``get_tophit`` is the only Foldseek-DataFrame producer and the heart of the
pandas->polars migration (REVIEW_FINDINGS B4 / POLAR_MIGRATION.md). Each test
pins one branch; the snapshots are the byte-identity gate the polars rewrite
must reproduce.

Run::

    pytest tests/unit/test_tophit.py -v
    pytest tests/unit/test_tophit.py --snapshot-update
"""
from __future__ import annotations

import polars as pl
import pytest

from baktfold.results.tophit import get_tophit

SUBDIR = "test_tophit"


# ─── basic non-structures path ──────────────────────────────────────────────


def test_get_tophit_basic_snapshot(fixtures_dir, df_snapshot):
    """Full output snapshot — the byte-identity gate (column order, qCov/tCov
    formatting, evalue scientific notation, dedup all in one)."""
    df = get_tophit(fixtures_dir / "foldseek_basic.tsv", structures=False, cath=False)
    df_snapshot(df, "basic.tsv", subdir=SUBDIR)


def test_get_tophit_basic_dedup_keeps_top_hit(fixtures_dir):
    """query_0001 has two hits; drop_duplicates(keep='first') keeps the top
    (highest-bitscore, first-listed) row → 3 rows for 3 queries."""
    df = get_tophit(fixtures_dir / "foldseek_basic.tsv", structures=False, cath=False)
    assert len(df) == 3
    assert df["query"].to_list() == ["query_0001", "query_0002", "query_0003"]
    # the survivor for query_0001 is the 520.0-bitscore target, not the 310.0 one
    q1 = df.filter(pl.col("query") == "query_0001").row(0, named=True)
    assert q1["target"] == "AF-P12345-F1-model_v4"
    assert float(q1["bitscore"]) == 520.0


def test_get_tophit_qcov_tcov_math(fixtures_dir):
    """qCov = round((qEnd - qStart) / qLen, 2); tCov analogous. No +1."""
    df = get_tophit(fixtures_dir / "foldseek_basic.tsv", structures=False, cath=False)
    q3 = df.filter(pl.col("query") == "query_0003").row(0, named=True)
    # (48 - 1) / 60 = 0.7833.. -> 0.78 ; (50 - 3) / 55 = 0.8545.. -> 0.85
    assert q3["qCov"] == pytest.approx(0.78)
    assert q3["tCov"] == pytest.approx(0.85)


def test_get_tophit_column_order(fixtures_dir):
    """The reorder must place qCov directly after qLen and the
    tStart/tEnd/tLen/tCov block together."""
    df = get_tophit(fixtures_dir / "foldseek_basic.tsv", structures=False, cath=False)
    cols = list(df.columns)
    assert cols.index("qCov") == cols.index("qLen") + 1
    assert cols.index("tCov") == cols.index("tLen") + 1
    assert cols == [
        "query", "target", "bitscore", "fident", "evalue",
        "qStart", "qEnd", "qLen", "qCov",
        "tStart", "tEnd", "tLen", "tCov",
    ]


# ─── structures path ────────────────────────────────────────────────────────


def test_get_tophit_structures_columns(fixtures_dir):
    """structures=True adds alntmscore + lddt (kept at the end after reorder)."""
    df = get_tophit(fixtures_dir / "foldseek_structures.tsv", structures=True, cath=False)
    assert "alntmscore" in df.columns
    assert "lddt" in df.columns
    cols = list(df.columns)
    assert cols[-2:] == ["alntmscore", "lddt"]


def test_get_tophit_structures_snapshot(fixtures_dir, df_snapshot):
    df = get_tophit(fixtures_dir / "foldseek_structures.tsv", structures=True, cath=False)
    df_snapshot(df, "structures.tsv", subdir=SUBDIR)


# ─── CATH greedy (no dedup) ─────────────────────────────────────────────────


def test_get_tophit_cath_keeps_all_hits(fixtures_dir):
    """cath=True skips drop_duplicates so multi-domain queries keep every
    greedy hit (query_0001 keeps both of its domain rows)."""
    df = get_tophit(fixtures_dir / "foldseek_cath.tsv", structures=False, cath=True)
    assert len(df) == 3
    assert (df["query"] == "query_0001").sum() == 2


def test_get_tophit_cath_snapshot(fixtures_dir, df_snapshot):
    df = get_tophit(fixtures_dir / "foldseek_cath.tsv", structures=False, cath=True)
    df_snapshot(df, "cath.tsv", subdir=SUBDIR)


# ─── ~PIPE~ round-trip ──────────────────────────────────────────────────────


def test_get_tophit_pipe_replaced(fixtures_dir):
    """`~PIPE~` markers in query must be restored to `|` early."""
    df = get_tophit(fixtures_dir / "foldseek_pipe.tsv", structures=False, cath=False)
    queries = df["query"].to_list()
    assert queries == ["contig|1_0001"]
    assert not any("~PIPE~" in q for q in queries)


# ─── empty foldseek result (no hits) ────────────────────────────────────────


def test_get_tophit_empty_returns_empty_no_crash(fixtures_dir):
    """A 0-byte Foldseek result (genome with no hits) must not raise — the
    `.empty` branch logs a warning and returns the (empty) frame."""
    df = get_tophit(fixtures_dir / "foldseek_empty.tsv", structures=False, cath=False)
    assert len(df) == 0
    # empty path skips qCov/tCov derivation — only the raw 11 columns exist
    assert "qCov" not in df.columns
