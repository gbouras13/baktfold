"""Chain-level PDB labels (audit A5, option 2b with 2a fallback).

A Foldseek PDB hit is a chain ('6xgw-assembly1.cif.gz_A'), but Baktfold used
to label it with the title of the whole PDB *entry*, e.g. "Crystal structure
of LeuBAT (delta5 mutant) in complex with mazindol" - not a protein name, and
for complexes possibly not even the chain that was hit. The database now
carries a ``pdb_chain`` table (UniProt name of the chain via SIFTS, else the
PDB entity name); databases without it keep the old title behaviour.

Run::

    pytest tests/unit/test_pdb_chain_names.py -v
"""
from __future__ import annotations

import sqlite3

import pytest

from baktfold.bakta import pstc


@pytest.mark.parametrize(
    "target, chain",
    [
        ("6xgw-assembly1.cif.gz_A", "A"),
        ("201l-assembly2.cif.gz_B", "B"),
        ("1abc-assembly1.cif.gz_A-2", "A"),     # symmetry copy -> base chain
        ("7xyz-assembly1.cif.gz_AAA", "AAA"),   # multi-character chain IDs
        ("6xgw", None),
    ],
)
def test_pdb_chain_from_target(target, chain):
    assert pstc.pdb_chain_from_target(target) == chain


def _db(tmp_path, with_chain_table: bool):
    d = tmp_path / "db"
    d.mkdir()
    conn = sqlite3.connect(d / "baktfold.db")
    for t in ("swissprot", "afdbclusters", "pdb"):
        conn.execute(f"CREATE TABLE {t} (id TEXT PRIMARY KEY, product TEXT)")
    conn.execute("CREATE TABLE cath (id TEXT PRIMARY KEY, cath_class TEXT, product TEXT)")
    conn.execute("INSERT INTO pdb VALUES ('101m', 'SPERM WHALE MYOGLOBIN F46V N-BUTYL ISOCYANIDE AT PH 9.0')")
    if with_chain_table:
        conn.execute("CREATE TABLE pdb_chain (id TEXT PRIMARY KEY, product TEXT, source TEXT)")
        conn.execute("INSERT INTO pdb_chain VALUES ('101m_A', 'Myoglobin', 'uniprot')")
    conn.commit()
    conn.close()
    return d


def _feat(chain):
    hit = {"source": "pdb", "id": "101m"}
    if chain:
        hit["chain"] = chain
    return [{"type": "cds", "locus": "c1", "pstc": [hit]}]


def test_chain_name_used_when_table_present(tmp_path):
    feats = _feat("A")
    pstc.lookup_sql(feats, _db(tmp_path, True), threads=1)
    assert feats[0]["pstc"][0]["description"] == "Myoglobin"


def test_unlisted_chain_is_not_given_the_entry_title(tmp_path):
    """Conservative: an unmapped chain falls through, it does not get the title."""
    feats = _feat("B")
    pstc.lookup_sql(feats, _db(tmp_path, True), threads=1)
    assert feats[0]["pstc"][0]["description"] == "hypothetical protein"


def test_old_database_keeps_entry_title(tmp_path):
    feats = _feat("A")
    pstc.lookup_sql(feats, _db(tmp_path, False), threads=1)
    assert feats[0]["pstc"][0]["description"].startswith("SPERM WHALE MYOGLOBIN")


def test_parse_records_chain():
    import polars as pl
    df = pl.DataFrame({
        "query": ["q1"], "target": ["101m-assembly1.cif.gz_A-2"],
        "bitscore": [100.0], "fident": [0.5], "evalue": ["1e-20"],
        "qStart": [1], "qEnd": [100], "qLen": [100], "qCov": [1.0],
        "tStart": [1], "tEnd": [100], "tLen": [100], "tCov": [1.0],
    })
    feats = [{"type": "cds", "locus": "q1", "id": "q1", "aa_hexdigest": None,
              "start": 1, "stop": 300, "strand": "+"}]
    pstc.parse(feats, df, "pdb")
    hit = feats[0]["pstc"][0]
    assert (hit["id"], hit["chain"]) == ("101m", "A")


def test_uniprot_accession_recorded_when_column_present(tmp_path):
    d = _db(tmp_path, False)
    conn = sqlite3.connect(d / "baktfold.db")
    conn.execute("CREATE TABLE pdb_chain (id TEXT PRIMARY KEY, product TEXT, source TEXT, uniprot TEXT)")
    conn.execute("INSERT INTO pdb_chain VALUES ('101m_A', 'Myoglobin', 'uniprot', 'P02185')")
    conn.commit(); conn.close()
    feats = _feat("A")
    pstc.lookup_sql(feats, d, threads=1)
    hit = feats[0]["pstc"][0]
    assert (hit["description"], hit["uniprot"]) == ("Myoglobin", "P02185")


def test_labelling_pdb_hit_gets_chain_and_uniprot_xrefs():
    import baktfold.bakta.annotation as anno
    feat = {"type": "cds", "locus": "c1", "product": "hypothetical protein", "hypothetical": True,
            "db_xrefs": [],
            "pstc": [{"source": "pdb", "id": "101m", "chain": "A", "uniprot": "P02185",
                      "description": "Myoglobin", "query_cov": 0.9, "subject_cov": 0.9,
                      "identity": 0.4, "evalue": 1e-20}]}
    anno.combine_annotation(feat, fast=False)
    assert feat["product"] == "Myoglobin"
    assert "pdb:pdb_101m_A" in feat["db_xrefs"]
    assert "UniProtKB:P02185" in feat["db_xrefs"]


def test_non_labelling_pdb_hit_gets_no_uniprot_xref():
    import baktfold.bakta.annotation as anno
    feat = {"type": "cds", "locus": "c1", "product": "hypothetical protein", "hypothetical": True,
            "db_xrefs": [],
            "pstc": [{"source": "swissprot", "id": "P1", "description": "DNA ligase",
                      "query_cov": 0.9, "subject_cov": 0.9, "identity": 0.4, "evalue": 1e-20},
                     {"source": "pdb", "id": "101m", "chain": "A", "uniprot": "P02185",
                      "description": "Myoglobin", "query_cov": 0.9, "subject_cov": 0.9,
                      "identity": 0.4, "evalue": 1e-20}]}
    anno.combine_annotation(feat, fast=False)
    assert feat["product"] == "DNA ligase"
    assert "pdb:pdb_101m_A" in feat["db_xrefs"]
    assert not any(x.startswith("UniProtKB:") for x in feat["db_xrefs"])
