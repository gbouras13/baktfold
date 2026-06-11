"""
Function-level tests for baktfold.bakta.pstc.

``parse`` is the Foldseek-DataFrame consumer (folds hits back into feature
dicts); ``fetch_sql_description`` / ``lookup_sql`` resolve PSTC accessions to
product descriptions from the SQLite DB.

These pin behaviour ahead of:
  - REVIEW_FINDINGS B1 (delete the dead double-build in ``parse``)
  - REVIEW_FINDINGS B2 (swap ``iterrows`` -> ``to_dicts``)
  - REVIEW_FINDINGS B3 (stop opening a new SQLite connection per accession)
  - the pandas->polars migration (``parse`` is the consumer side)

NOTE: MIN_PSTC_IDENTITY / QCOVERAGE / TCOVERAGE are all 0.0 today, so every
matched hit passes the filter. These tests pin that current behaviour.

Run::

    pytest tests/unit/test_pstc.py -v
"""
from __future__ import annotations

import sqlite3

import polars as pl
import pytest

from baktfold.bakta import pstc


# ─── helpers ────────────────────────────────────────────────────────────────

_HITS_COLS = ["query", "target", "bitscore", "fident", "evalue", "qCov", "tCov"]


def _hits_df(rows: list[dict]) -> pl.DataFrame:
    """Build a post-get_tophit-shaped polars DataFrame (the columns ``parse`` reads)."""
    if not rows:
        return pl.DataFrame(schema={c: pl.Utf8 for c in _HITS_COLS})
    return pl.DataFrame(rows).select(_HITS_COLS)


def _feat(locus: str, **extra) -> dict:
    f = {"locus": locus, "id": locus}
    f.update(extra)
    return f


# ─── parse: empty input ─────────────────────────────────────────────────────


def test_parse_empty_df_returns_features_untouched():
    feats = [_feat("cds_1")]
    out = pstc.parse(feats, _hits_df([]), db_name="swissprot")
    assert out is feats
    assert "pstc" not in feats[0]


# ─── parse: accession extraction per database ───────────────────────────────


def test_parse_swissprot_accession_split():
    """swissprot target 'AF-P12345-F1-model_v4' -> accession 'P12345' (split('-')[1])."""
    df = _hits_df([dict(query="cds_1", target="AF-P12345-F1-model_v4",
                        bitscore=500.0, fident=0.9, evalue=1e-40, qCov=0.95, tCov=0.95)])
    feats = [_feat("cds_1")]
    pstc.parse(feats, df, db_name="swissprot")
    assert feats[0]["pstc"] == [dict(
        source="swissprot", id="P12345",
        query_cov=0.95, subject_cov=0.95, identity=0.9, score=500.0, evalue=1e-40,
    )]


def test_parse_afdb_accession_split():
    df = _hits_df([dict(query="cds_1", target="AF-Q99999-F1-model_v4",
                        bitscore=400.0, fident=0.8, evalue=1e-30, qCov=0.9, tCov=0.9)])
    feats = [_feat("cds_1")]
    pstc.parse(feats, df, db_name="afdb")
    assert feats[0]["pstc"][0]["id"] == "Q99999"
    assert feats[0]["pstc"][0]["source"] == "afdb"


def test_parse_pdb_accession_split_first():
    """pdb uses split('-')[0]."""
    df = _hits_df([dict(query="cds_1", target="8ABC-A",
                        bitscore=300.0, fident=0.7, evalue=1e-20, qCov=0.8, tCov=0.8)])
    feats = [_feat("cds_1")]
    pstc.parse(feats, df, db_name="pdb")
    assert feats[0]["pstc"][0]["id"] == "8ABC"


def test_parse_cath_accession_is_whole_target():
    df = _hits_df([dict(query="cds_1", target="1abcA00",
                        bitscore=250.0, fident=0.6, evalue=1e-18, qCov=0.7, tCov=0.7)])
    feats = [_feat("cds_1")]
    pstc.parse(feats, df, db_name="cath")
    assert feats[0]["pstc"][0]["id"] == "1abcA00"


# ─── parse: multiple CATH hits append into a list ───────────────────────────


def test_parse_cath_multiple_hits_append():
    """A multi-domain query with two greedy hits -> pstc is a list of two."""
    df = _hits_df([
        dict(query="cds_1", target="1abcA00", bitscore=300.0, fident=0.6, evalue=1e-20, qCov=0.3, tCov=0.98),
        dict(query="cds_1", target="2defB01", bitscore=250.0, fident=0.55, evalue=1e-18, qCov=0.35, tCov=0.96),
    ])
    feats = [_feat("cds_1")]
    pstc.parse(feats, df, db_name="cath")
    assert isinstance(feats[0]["pstc"], list)
    assert [p["id"] for p in feats[0]["pstc"]] == ["1abcA00", "2defB01"]


# ─── parse: locus vs id keying ──────────────────────────────────────────────


def test_parse_has_duplicate_locus_keys_on_id():
    """has_duplicate_locus=True -> match on feat['id'], not feat['locus']."""
    df = _hits_df([dict(query="unique_id_1", target="AF-P11111-F1-model_v4",
                        bitscore=500.0, fident=0.9, evalue=1e-40, qCov=0.95, tCov=0.95)])
    feats = [{"locus": "shared_locus", "id": "unique_id_1"}]
    pstc.parse(feats, df, db_name="swissprot", has_duplicate_locus=True)
    assert feats[0]["pstc"][0]["id"] == "P11111"


def test_parse_non_matching_feature_untouched():
    df = _hits_df([dict(query="cds_1", target="AF-P12345-F1-model_v4",
                        bitscore=500.0, fident=0.9, evalue=1e-40, qCov=0.95, tCov=0.95)])
    feats = [_feat("cds_1"), _feat("cds_2")]
    pstc.parse(feats, df, db_name="swissprot")
    assert "pstc" in feats[0]
    assert "pstc" not in feats[1]


def test_parse_zero_coverage_still_passes_today():
    """Thresholds are all 0.0, so even a 0.0-coverage hit is recorded. Pins
    current behaviour — change this test deliberately if thresholds change."""
    df = _hits_df([dict(query="cds_1", target="AF-P12345-F1-model_v4",
                        bitscore=10.0, fident=0.0, evalue=1.0, qCov=0.0, tCov=0.0)])
    feats = [_feat("cds_1")]
    pstc.parse(feats, df, db_name="swissprot")
    assert feats[0]["pstc"][0]["id"] == "P12345"


# ─── fetch_sql_description: table routing ───────────────────────────────────


@pytest.fixture
def sqlite_conn():
    """In-memory SQLite with one product row per source table."""
    conn = sqlite3.connect(":memory:")
    for table in ("swissprot", "afdbclusters", "pdb", "cath"):
        conn.execute(f"CREATE TABLE {table} (id TEXT, product TEXT)")
    conn.execute("INSERT INTO swissprot VALUES ('P12345', 'DNA polymerase')")
    conn.execute("INSERT INTO afdbclusters VALUES ('Q99999', 'tail fiber protein')")
    conn.execute("INSERT INTO pdb VALUES ('8ABC', 'capsid protein')")
    conn.execute("INSERT INTO cath VALUES ('1abcA00', 'lysozyme domain')")
    conn.commit()
    yield conn
    conn.close()


@pytest.mark.parametrize("source,accession,expected", [
    ("swissprot", "P12345", "DNA polymerase"),
    ("afdb", "Q99999", "tail fiber protein"),   # afdb -> afdbclusters table
    ("pdb", "8ABC", "capsid protein"),
    ("cath", "1abcA00", "lysozyme domain"),
])
def test_fetch_sql_description_routing(sqlite_conn, source, accession, expected):
    assert pstc.fetch_sql_description(sqlite_conn, source, accession) == expected


def test_fetch_sql_description_miss_returns_none(sqlite_conn):
    assert pstc.fetch_sql_description(sqlite_conn, "swissprot", "NOPE") is None


def test_fetch_sql_description_unknown_source_returns_none(sqlite_conn):
    assert pstc.fetch_sql_description(sqlite_conn, "custom_db", "anything") is None


# ─── lookup_sql: end-to-end against a real baktfold.db file ─────────────────


@pytest.fixture
def baktfold_db_dir(tmp_path):
    """A directory containing a minimal baktfold.db (what lookup_sql opens)."""
    db_path = tmp_path / "baktfold.db"
    conn = sqlite3.connect(str(db_path))
    for table in ("swissprot", "afdbclusters", "pdb", "cath"):
        conn.execute(f"CREATE TABLE {table} (id TEXT, product TEXT)")
    conn.execute("INSERT INTO swissprot VALUES ('P12345', 'DNA polymerase')")
    conn.commit()
    conn.close()
    return tmp_path


def test_lookup_sql_fills_description_on_hit(baktfold_db_dir):
    feats = [{"locus": "cds_1", "id": "cds_1",
              "pstc": {"source": "swissprot", "id": "P12345"}}]
    pstc.lookup_sql(feats, baktfold_db_dir, threads=1)
    entry = feats[0]["pstc"]
    entry = entry[0] if isinstance(entry, list) else entry
    assert entry["description"] == "DNA polymerase"


def test_lookup_sql_miss_non_custom_is_hypothetical(baktfold_db_dir):
    feats = [{"locus": "cds_1", "id": "cds_1",
              "pstc": {"source": "swissprot", "id": "MISSING"}}]
    pstc.lookup_sql(feats, baktfold_db_dir, threads=1)
    entry = feats[0]["pstc"]
    entry = entry[0] if isinstance(entry, list) else entry
    assert entry["description"] == "hypothetical protein"


def test_lookup_sql_feature_without_pstc_skipped(baktfold_db_dir):
    feats = [{"locus": "cds_1", "id": "cds_1"}]  # no 'pstc'
    pstc.lookup_sql(feats, baktfold_db_dir, threads=1)
    assert "pstc" not in feats[0]
