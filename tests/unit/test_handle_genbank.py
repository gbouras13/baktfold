"""
Function-level tests for baktfold.io.handle_genbank.

Covers:
  - get_genbank input-style detection (Pharokka / Bakta / NCBI) and the
    fall-through paths (REVIEW_FINDINGS #9 — these tests pin the CURRENT
    behaviour, which is buggy for the unrecognised-style case; update them
    when #9 is fixed).
  - get_proteins FASTA->dict, plain + gzipped (REVIEW_FINDINGS #36).
  - identify_long_ids space-stripping for >54-char IDs.

GenBank fixtures are built in-test with BioPython so no binary files are
checked in.

Run::

    pytest tests/unit/test_handle_genbank.py -v
"""
from __future__ import annotations

import gzip

import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

from baktfold.io.handle_genbank import get_genbank, get_proteins, identify_long_ids


# ─── genbank fixture builder ────────────────────────────────────────────────


def _write_genbank(path, cds_qualifiers, comment=None, record_id="contig_1", with_cds=True):
    """Write a one-record GenBank file with a single CDS carrying the given
    qualifiers (and optional COMMENT). Returns the string path."""
    feats = []
    if with_cds:
        cds = SeqFeature(FeatureLocation(0, 30, strand=1), type="CDS")
        cds.qualifiers = dict(cds_qualifiers)
        feats.append(cds)
    rec = SeqRecord(Seq("ATG" * 10), id=record_id, name=record_id, description="test", features=feats)
    rec.annotations["molecule_type"] = "DNA"
    if comment is not None:
        rec.annotations["comment"] = comment
    SeqIO.write(rec, str(path), "genbank")
    return str(path)


# ─── get_genbank: style detection ───────────────────────────────────────────


def test_get_genbank_detects_pharokka(tmp_path):
    """phrog + ID qualifiers -> Pharokka."""
    p = _write_genbank(tmp_path / "phar.gbk",
                       {"ID": ["contig_1_CDS_0001"], "phrog": ["1234"],
                        "function": ["unknown function"], "product": ["hypothetical protein"]})
    gb_dict, method = get_genbank(p)
    assert method == "Pharokka"
    assert "contig_1" in gb_dict


def test_get_genbank_detects_bakta(tmp_path):
    """'Bakta' in COMMENT + locus_tag qualifier -> Bakta."""
    p = _write_genbank(tmp_path / "bakta.gbk",
                       {"locus_tag": ["LOCUS_0001"]},
                       comment="Annotated using Bakta")
    gb_dict, method = get_genbank(p)
    assert method == "Bakta"


def test_get_genbank_detects_ncbi(tmp_path):
    """protein_id and no phrog -> NCBI."""
    p = _write_genbank(tmp_path / "ncbi.gbk",
                       {"protein_id": ["ABC12345.1"]})
    gb_dict, method = get_genbank(p)
    assert method == "NCBI"


def test_get_genbank_no_cds_returns_dict_and_none(tmp_path):
    """Valid GenBank with no CDS -> (gb_dict, None)."""
    p = _write_genbank(tmp_path / "nocds.gbk", {}, with_cds=False)
    gb_dict, method = get_genbank(p)
    assert method is None
    assert "contig_1" in gb_dict


def test_get_genbank_unrecognised_style_current_behaviour(tmp_path):
    """A valid GenBank whose CDS matches no known style. CURRENT behaviour
    (REVIEW_FINDINGS #9): the unset `method` triggers UnboundLocalError which
    the broad except swallows -> ({}, None). Update this test when #9 is fixed
    (it should then return (gb_dict, None) with the real records)."""
    p = _write_genbank(tmp_path / "weird.gbk", {"gene": ["xyz"]})
    gb_dict, method = get_genbank(p)
    assert method is None
    assert gb_dict == {}  # ← buggy: records are lost. Flip to non-empty after #9.


def test_get_genbank_non_genbank_returns_empty(tmp_path):
    """A plain non-GenBank text file -> ({}, None)."""
    p = tmp_path / "plain.txt"
    p.write_text("this is not a genbank file\n")
    gb_dict, method = get_genbank(str(p))
    assert gb_dict == {}
    assert method is None


# ─── get_proteins ───────────────────────────────────────────────────────────


def test_get_proteins_plain(fixtures_dir):
    d = get_proteins(str(fixtures_dir / "db_aa.fasta"))
    assert d == {"cds_1": "MKLV", "cds_2": "MKWY", "cds_3": "MKGG"}


def test_get_proteins_gzipped(fixtures_dir, tmp_path):
    gz = tmp_path / "aa.fasta.gz"
    with gzip.open(gz, "wt") as f:
        f.write((fixtures_dir / "db_aa.fasta").read_text())
    d = get_proteins(str(gz))
    assert d == {"cds_1": "MKLV", "cds_2": "MKWY", "cds_3": "MKGG"}


# ─── identify_long_ids ──────────────────────────────────────────────────────


def test_identify_long_ids_strips_space_in_long_id():
    """A >=54-char ID that BioPython line-wrapped (introducing a space) has the
    space removed."""
    long_id_with_space = "contig_1_CDS_" + "A" * 40 + " " + "BBBB"  # >54 chars, has a space
    assert len(long_id_with_space) >= 54
    cds = SeqFeature(FeatureLocation(0, 30, strand=1), type="CDS")
    cds.qualifiers = {"ID": [long_id_with_space]}
    rec = SeqRecord(Seq("ATG" * 10), id="contig_1", features=[cds])
    gb_dict = {"contig_1": rec}

    out = identify_long_ids(gb_dict)
    fixed = out["contig_1"].features[0].qualifiers["ID"][0]
    assert " " not in fixed
    assert fixed == long_id_with_space.replace(" ", "")
