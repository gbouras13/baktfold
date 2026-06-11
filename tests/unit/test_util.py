"""
Function-level tests for baktfold.utils.util.

Covers:
  - replace_pipe_in_fasta: ~PIPE~ -> | restore in headers (REVIEW_FINDINGS
    #4/#25 — gates the atomic-write + streaming rewrite).
  - remove_directory / remove_file (REVIEW_FINDINGS #35).
  - get_type_rank / sort_euk_feature_key: pure eukaryotic feature ordering.

Run::

    pytest tests/unit/test_util.py -v
"""
from __future__ import annotations

from Bio import SeqIO

import baktfold.bakta.constants as bc
from baktfold.utils.util import (
    get_type_rank,
    remove_directory,
    remove_file,
    replace_pipe_in_fasta,
    sort_euk_feature_key,
)


# ─── replace_pipe_in_fasta ──────────────────────────────────────────────────


def test_replace_pipe_in_fasta_restores_pipe(tmp_path):
    """`~PIPE~` in a header becomes `|`; the sequence is preserved."""
    fasta = tmp_path / "in.fasta"
    fasta.write_text(">contig~PIPE~1_CDS_0001 some desc\nMKLVMKLV\n")

    replace_pipe_in_fasta(str(fasta))

    records = list(SeqIO.parse(str(fasta), "fasta"))
    assert len(records) == 1
    assert records[0].id == "contig|1_CDS_0001"
    assert "~PIPE~" not in records[0].id
    assert str(records[0].seq) == "MKLVMKLV"


def test_replace_pipe_in_fasta_no_pipe_unchanged(tmp_path):
    fasta = tmp_path / "in.fasta"
    fasta.write_text(">clean_header\nMKWY\n")
    replace_pipe_in_fasta(str(fasta))
    records = list(SeqIO.parse(str(fasta), "fasta"))
    assert records[0].id == "clean_header"
    assert str(records[0].seq) == "MKWY"


# ─── remove_directory / remove_file ─────────────────────────────────────────


def test_remove_directory_removes_tree(tmp_path):
    d = tmp_path / "to_remove"
    d.mkdir()
    (d / "f.txt").write_text("x")
    remove_directory(d)
    assert not d.exists()


def test_remove_directory_missing_is_noop(tmp_path):
    remove_directory(tmp_path / "does_not_exist")  # must not raise


def test_remove_file_removes_and_tolerates_missing(tmp_path):
    f = tmp_path / "f.txt"
    f.write_text("x")
    remove_file(f)
    assert not f.exists()
    remove_file(f)  # second call: missing, must not raise


# ─── get_type_rank ──────────────────────────────────────────────────────────


def test_get_type_rank_base_order():
    assert get_type_rank({"type": bc.FEATURE_GENE}) == 0
    assert get_type_rank({"type": "mRNA"}) == 1
    assert get_type_rank({"type": bc.FEATURE_CDS}) == 3
    assert get_type_rank({"type": bc.FEATURE_T_RNA}) == 6


def test_get_type_rank_unknown_is_99():
    assert get_type_rank({"type": "source"}) == 99


def test_get_type_rank_utr_strand_dependent():
    # 5'UTR: 2 on +, 4 on -
    assert get_type_rank({"type": bc.FEATURE_5UTR, "strand": "+"}) == 2
    assert get_type_rank({"type": bc.FEATURE_5UTR, "strand": "-"}) == 4
    # 3'UTR: 4 on +, 2 on -
    assert get_type_rank({"type": bc.FEATURE_3UTR, "strand": "+"}) == 4
    assert get_type_rank({"type": bc.FEATURE_3UTR, "strand": "-"}) == 2


# ─── sort_euk_feature_key ───────────────────────────────────────────────────


def test_sort_euk_feature_key_orders_within_locus():
    """Within one locus, features sort gene -> mRNA -> CDS by type rank."""
    gene = {"type": bc.FEATURE_GENE, "start": 100, "stop": 400, "locus": "L1"}
    mrna = {"type": "mRNA", "start": 100, "stop": 400, "locus": "L1"}
    cds = {"type": bc.FEATURE_CDS, "start": 100, "stop": 400, "locus": "L1"}
    keys = sorted([cds, mrna, gene], key=sort_euk_feature_key)
    assert [f["type"] for f in keys] == [bc.FEATURE_GENE, "mRNA", bc.FEATURE_CDS]


def test_sort_euk_feature_key_non_locus_sorts_by_start():
    a = {"type": "source", "start": 500}
    b = {"type": "source", "start": 100}
    assert sorted([a, b], key=sort_euk_feature_key) == [b, a]
