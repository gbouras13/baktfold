"""``compare`` must search the same (masked) amino acids as ``run``.

``baktfold predict`` writes ``{prefix}_aa.fasta`` with low-confidence residues
masked to ``X``, matching the masked 3Di. The input parsers (``parse_json_input``
/ ``parse_protein_input``) write an *unmasked* ``{prefix}_aa.fasta`` into the
output directory before ``subcommand_compare`` runs, and compare only copied the
predictions-dir file when the output file did not already exist. So
``compare``/``proteins-compare`` built the Foldseek query DB from unmasked AA +
masked 3Di, while ``run``/``proteins`` used masked AA + masked 3Di — different
hits and E-values from the same predictions (visible in the run vs compare
goldens: consistently higher bitscores via compare).

Run::

    pytest tests/unit/test_compare_masked_aa.py -v
"""
from __future__ import annotations

from unittest.mock import patch

import pytest

from baktfold.subcommands import compare as compare_mod
from baktfold.subcommands.compare import _copy_prediction, subcommand_compare


class _Captured(Exception):
    pass


def _run_compare(tmp_path, output, predictions_dir):
    logdir = output / "logs"
    logdir.mkdir(parents=True, exist_ok=True)
    seen = {}

    def _capture(fasta_aa, fasta_3di, *args, **kwargs):
        seen["aa"] = open(fasta_aa).read()
        seen["3di"] = open(fasta_3di).read()
        raise _Captured

    with patch.object(compare_mod, "generate_foldseek_db_from_aa_3di", _capture):
        with pytest.raises(_Captured):
            subcommand_compare(
                [{"locus": "cds_1", "id": "cds_1", "aa": "MKVLA"}],
                output,
                threads=1, evalue=1e-3, sensitivity=9.5,
                database=tmp_path / "db", prefix="genome",
                predictions_dir=predictions_dir, structures=False,
                structure_dir=None, logdir=logdir, proteins_flag=False,
                max_seqs=1000, ultra_sensitive=False, extra_foldseek_params="",
                custom_db=None, foldseek_gpu=False, custom_annotations=None,
                has_duplicate_locus=False, fast=True,
            )
    return seen


def test_compare_overwrites_unmasked_aa_with_masked_prediction(tmp_path):
    pred = tmp_path / "predictions"
    pred.mkdir()
    (pred / "genome_aa.fasta").write_text(">cds_1\nMKXXA\n")      # masked by predict
    (pred / "genome_3di.fasta").write_text(">cds_1\nDVXXL\n")

    output = tmp_path / "output"
    output.mkdir()
    # what the input parser leaves behind before compare runs
    (output / "genome_aa.fasta").write_text(">cds_1\nMKVLA\n")

    seen = _run_compare(tmp_path, output, pred)

    assert seen["aa"] == ">cds_1\nMKXXA\n"
    assert seen["3di"] == ">cds_1\nDVXXL\n"


def test_run_layout_predictions_in_output_dir_untouched(tmp_path):
    """`baktfold run`: predictions_dir == output, so the copy is a no-op."""
    output = tmp_path / "output"
    output.mkdir()
    (output / "genome_aa.fasta").write_text(">cds_1\nMKXXA\n")
    (output / "genome_3di.fasta").write_text(">cds_1\nDVXXL\n")

    seen = _run_compare(tmp_path, output, output)

    assert seen["aa"] == ">cds_1\nMKXXA\n"


def test_copy_prediction_same_file_is_noop(tmp_path):
    f = tmp_path / "x.fasta"
    f.write_text(">a\nM\n")
    _copy_prediction(f, f)
    assert f.read_text() == ">a\nM\n"
