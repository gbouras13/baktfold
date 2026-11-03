"""Tests for the empty-query short circuit in ``subcommand_compare``.

A genome whose CDS all already carry a product (e.g. an Ensembl protist
assembly annotated upstream) has **zero** hypothetical proteins, so
``baktfold predict`` writes a 0-byte ``<prefix>_3di.fasta`` and
``<prefix>_aa.fasta``. Before this guard, compare only checked that those
files *existed*, built an empty Foldseek query DB from them, and the search
died on a zero-dimension kernel launch::

    CUDA error: invalid configuration argument : .../libmarv/src/pssm.cuh, line 346
    Error: Ungapped prefilter matching step died

which is a thoroughly misleading way to say "there was nothing to search".
The CPU path did not crash but still spent the whole search on an empty DB.

These tests pin that compare now returns cleanly (no Foldseek invocation, no
exception, features handed back untouched so the caller still writes the
normal annotation outputs) — and, importantly, that the guard does *not* fire
for a genome that does have hypothetical proteins.

No database, no foldseek, no GPU needed.

Run::

    pytest tests/unit/test_empty_query.py -v
"""
from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

import pytest

from baktfold.subcommands import compare as compare_mod
from baktfold.subcommands.compare import fasta_is_empty, subcommand_compare


class _FoldseekReached(Exception):
    """Raised by the stubbed DB builder to prove we got past the guard."""


def _explode(name):
    def _stub(*args, **kwargs):  # pragma: no cover - must not run
        raise AssertionError(f"{name} was called for an empty query set")
    return _stub


def _call_compare(tmp_path, *, hypotheticals, structures=False, prefix="genome"):
    """Run subcommand_compare with every Foldseek entry point stubbed out.

    ``generate_foldseek_db_*`` raises _FoldseekReached so a test can assert
    execution *did* reach the Foldseek stage; everything downstream asserts
    loudly if it is reached at all.
    """
    output = tmp_path / "output"
    output.mkdir(exist_ok=True)
    logdir = output / "logs"
    logdir.mkdir(exist_ok=True)

    def _reached(*args, **kwargs):
        raise _FoldseekReached

    with patch.object(compare_mod, "generate_foldseek_db_from_aa_3di", _reached), \
         patch.object(compare_mod, "generate_foldseek_db_from_structures", _reached), \
         patch.object(compare_mod, "run_foldseek_search", _explode("run_foldseek_search")), \
         patch.object(compare_mod, "create_result_tsv", _explode("create_result_tsv")), \
         patch.object(compare_mod, "summarise_hits", _explode("summarise_hits")), \
         patch.object(compare_mod, "get_tophit", _explode("get_tophit")):
        return subcommand_compare(
            hypotheticals,
            output,
            threads=1,
            evalue=1e-3,
            sensitivity=9.5,
            database=tmp_path / "db",
            prefix=prefix,
            predictions_dir=tmp_path / "predictions",
            structures=structures,
            structure_dir=None,
            logdir=logdir,
            proteins_flag=False,
            max_seqs=1000,
            ultra_sensitive=False,
            extra_foldseek_params="",
            custom_db=None,
            foldseek_gpu=True,
            custom_annotations=None,
            has_duplicate_locus=False,
            fast=True,
        )


def _write_predictions(tmp_path, *, aa: str, threedi: str, prefix="genome"):
    """Lay down the AA + 3Di FASTAs baktfold predict would have written."""
    pred_dir = tmp_path / "predictions"
    pred_dir.mkdir(exist_ok=True)
    (pred_dir / f"{prefix}_aa.fasta").write_text(aa)
    (pred_dir / f"{prefix}_3di.fasta").write_text(threedi)
    return pred_dir


# ─── fasta_is_empty ─────────────────────────────────────────────────────────


def test_fasta_is_empty_on_zero_byte_file(tmp_path):
    p = tmp_path / "genome_3di.fasta"
    p.write_text("")
    assert fasta_is_empty(p) is True


def test_fasta_is_empty_on_missing_file(tmp_path):
    assert fasta_is_empty(tmp_path / "nope.fasta") is True


def test_fasta_is_not_empty_with_one_record(tmp_path):
    p = tmp_path / "genome_3di.fasta"
    p.write_text(">cds_1\nACDE\n")
    assert fasta_is_empty(p) is False


# ─── the short circuit ──────────────────────────────────────────────────────


def test_zero_hypotheticals_skips_foldseek_entirely(tmp_path):
    """The Giardia case: 6 CDS, all with real products, 0 hypothetical."""
    _write_predictions(tmp_path, aa="", threedi="")

    result = _call_compare(tmp_path, hypotheticals=[])

    assert result == []
    # short circuit happens before the query DB directory is even created
    assert not (tmp_path / "output" / "foldseek_db").exists()


def test_empty_3di_with_non_empty_aa_skips_foldseek(tmp_path):
    """ProstT5 produced nothing (every prediction OOM'd) — a Foldseek DB built
    from AA with no matching 3Di is equally useless, so skip it too."""
    _write_predictions(tmp_path, aa=">cds_1\nMKV\n", threedi="")

    feats = [{"locus": "cds_1", "id": "cds_1", "aa": "MKV"}]
    result = _call_compare(tmp_path, hypotheticals=feats)

    assert result == feats
    assert not (tmp_path / "output" / "foldseek_db").exists()


def test_structures_input_with_no_hypotheticals_skips_foldseek(tmp_path):
    """The structures path writes its own AA FASTA from the feature list, so
    zero features means a 0-byte query there too."""
    (tmp_path / "predictions").mkdir(exist_ok=True)

    result = _call_compare(tmp_path, hypotheticals=[], structures=True)

    assert result == []
    assert not (tmp_path / "output" / "foldseek_db").exists()


def test_short_circuit_leaves_features_untouched(tmp_path):
    """Nothing is annotated, so the caller must get its features back exactly
    as passed — `baktfold run` goes on to write the normal outputs from them."""
    _write_predictions(tmp_path, aa="", threedi="")
    feats = [{"locus": "cds_1", "id": "cds_1", "product": "DNA gyrase subunit A"}]

    result = _call_compare(tmp_path, hypotheticals=feats)

    assert result is feats
    assert result[0] == {
        "locus": "cds_1", "id": "cds_1", "product": "DNA gyrase subunit A",
    }


# ─── the guard must not be over-eager ───────────────────────────────────────


def test_non_empty_query_still_reaches_foldseek(tmp_path):
    """A normal genome must sail past the guard into the Foldseek stage."""
    _write_predictions(tmp_path, aa=">cds_1\nMKV\n", threedi=">cds_1\nACD\n")

    feats = [{"locus": "cds_1", "id": "cds_1", "aa": "MKV"}]
    with pytest.raises(_FoldseekReached):
        _call_compare(tmp_path, hypotheticals=feats)


def test_non_empty_structures_query_still_reaches_foldseek(tmp_path):
    (tmp_path / "predictions").mkdir(exist_ok=True)

    feats = [{"locus": "cds_1", "id": "cds_1", "aa": "MKV"}]
    with pytest.raises(_FoldseekReached):
        _call_compare(tmp_path, hypotheticals=feats, structures=True)
