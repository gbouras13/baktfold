"""`baktfold proteins` on a Bakta proteins JSON must only annotate hypotheticals.

It used to build the ProstT5 input from *every* protein in the JSON, so a
93-protein Bakta output with 17 hypotheticals sent all 93 through ProstT5 and
Foldseek, and combine_annotation then overwrote Bakta's curated products
wherever a structural hit was found -- contradicting the log message and
docs/index.md ("Only hypothetical proteins ... will be annotated").

Run::

    pytest tests/unit/test_proteins_hypotheticals_only.py -v
"""
from __future__ import annotations

from pathlib import Path

import baktfold.bakta.annotation as anno
from baktfold.io.fasta_in import parse_protein_input, proteins_to_predict

DATA = Path(__file__).resolve().parents[1] / "test_data"
BAKTA_JSON = DATA / "assembly_bakta_proteins_output_all" / "assembly.json"
BAKTA_FAA = DATA / "assembly_bakta_proteins_output_all" / "assembly.faa"


def test_bakta_json_predicts_only_hypotheticals(tmp_path):
    aas, _ = parse_protein_input(BAKTA_JSON, tmp_path / "aa.fasta")
    to_predict = proteins_to_predict(aas)

    n_hyp = sum(1 for aa in aas if aa.get("hypothetical"))
    assert len(aas) > n_hyp > 0
    assert set(to_predict) == {aa["locus"] for aa in aas if aa.get("hypothetical")}


def test_fasta_predicts_everything(tmp_path):
    aas, _ = parse_protein_input(BAKTA_FAA, tmp_path / "aa.fasta")
    assert len(proteins_to_predict(aas)) == len(aas) > 0


def test_curated_product_survives_annotation(tmp_path):
    """A non-hypothetical protein gets no pstc hits, so its product is kept."""
    aas, _ = parse_protein_input(BAKTA_JSON, tmp_path / "aa.fasta")
    curated = next(aa for aa in aas if not aa.get("hypothetical"))
    before = curated["product"]
    anno.combine_annotation(curated, fast=False)
    assert curated["product"] == before
    assert "hypothetical" not in curated
