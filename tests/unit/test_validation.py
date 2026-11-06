"""
Function-level tests for baktfold.utils.validation.

Covers the pure ``is_fasta`` detector (valid / non-fasta / gzipped / empty /
bad-character).

NOTE: ``check_dependencies`` (REVIEW_FINDINGS #7/#8 — the foldseek
UnboundLocalError and the unconditional "version is ok" log) is NOT unit-tested
here because the version logic is currently inlined and shells out to foldseek.
Once #8 extracts a pure version-parsing helper, add a test driving it with
string inputs (v10 -> ok, other major -> warn, garbage -> warn).

Run::

    pytest tests/unit/test_validation.py -v
"""
from __future__ import annotations

import gzip

from baktfold.utils.validation import is_fasta


def test_is_fasta_valid(tmp_path):
    p = tmp_path / "ok.fasta"
    p.write_text(">cds_1\nMKLVMKLV\n>cds_2\nMKWY\n")
    assert is_fasta(p) is True


def test_is_fasta_non_fasta(tmp_path):
    p = tmp_path / "no.txt"
    p.write_text("just some text\nnot fasta\n")
    assert is_fasta(p) is False


def test_is_fasta_empty(tmp_path):
    p = tmp_path / "empty.fasta"
    p.write_text("")
    assert is_fasta(p) is False


def test_is_fasta_gzipped(tmp_path):
    p = tmp_path / "ok.fasta.gz"
    with gzip.open(p, "wt") as f:
        f.write(">cds_1\nMKLV\n")
    assert is_fasta(p) is True


def test_is_fasta_rejects_bad_sequence_char(tmp_path):
    """A digit in the sequence body -> not a valid amino-acid FASTA."""
    p = tmp_path / "bad.fasta"
    p.write_text(">cds_1\nMK1V\n")  # '1' is not alpha / -*.
    assert is_fasta(p) is False


def test_is_fasta_allows_gap_and_stop_chars(tmp_path):
    """'-', '*', '.' are explicitly allowed in sequence lines."""
    p = tmp_path / "gaps.fasta"
    p.write_text(">cds_1\nMK-L*V.\n")
    assert is_fasta(p) is True
