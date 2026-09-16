"""
Regression tests for issue #40: baktfold crashed writing the feature TSV with

    TypeError: '<' not supported between instances of 'NoneType' and 'str'

convert_rrna_feature() built ``db_xrefs`` as ``[SO_RRNA.id, specific_so]``,
where ``specific_so`` is a dict lookup on the exact product string. Anything
outside "16S/23S/5S ribosomal RNA" gave None, so the list held a None next to
a str and every writer that sorts db_xrefs blew up. The reporter's archaeal
MAG contained one ``5.8S ribosomal RNA`` (barrnap emits 5.8S in eukaryote
mode; it turns up in MAGs), out of 1119 features.

Run::

    pytest tests/unit/test_rrna_db_xrefs.py -v
"""
from __future__ import annotations

from pathlib import Path

import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord

import baktfold.bakta.constants as bc
from baktfold.io.prokka_gbk_to_json import convert_rrna_feature
from baktfold.io.tsv import write_features


def _rrna(product: str) -> dict:
    feature = SeqFeature(
        FeatureLocation(0, 120, strand=1),
        type="rRNA",
        qualifiers={"product": [product], "locus_tag": ["TEST_00001"]},
    )
    rec = SeqRecord(Seq("ACGT" * 40), id="contig_1")
    return convert_rrna_feature(feature, rec, "rrna_1")


@pytest.mark.parametrize("product", [
    "5.8S ribosomal RNA",     # the product from issue #40
    "18S ribosomal RNA",
    "28S ribosomal RNA",
    "some unclassified rRNA",
])
def test_unmapped_rrna_product_yields_no_none_xref(product):
    """An unmapped product gets the generic rRNA term only — never a None."""
    xrefs = _rrna(product)["db_xrefs"]
    assert None not in xrefs
    assert xrefs == [bc.SO_RRNA] if hasattr(bc, "SO_RRNA") else xrefs == ["SO:0000252"]
    sorted(xrefs)   # must not raise — this is what crashed


@pytest.mark.parametrize("product,expected", [
    ("16S ribosomal RNA", "SO:0001000"),
    ("23S ribosomal RNA", "SO:0001001"),
    ("5S ribosomal RNA", "SO:0000652"),
])
def test_mapped_rrna_products_keep_their_specific_term(product, expected):
    xrefs = _rrna(product)["db_xrefs"]
    assert xrefs == ["SO:0000252", expected]


@pytest.mark.parametrize("product,expected", [
    ("16S ribosomal RNA (partial)", "SO:0001000"),
    ("23S ribosomal RNA (partial)", "SO:0001001"),
])
def test_partial_predictions_still_map(product, expected):
    """barrnap marks truncated genes '(partial)' — common in MAGs, and it would
    otherwise have hit the same None path."""
    assert _rrna(product)["db_xrefs"] == ["SO:0000252", expected]


def test_write_features_survives_a_none_db_xref(tmp_path):
    """Defence in depth: a malformed db_xref from any input converter must not
    take down output writing."""
    sequences = [{"id": "contig_1"}]
    features_by_sequence = {"contig_1": [{
        "sequence": "contig_1", "type": bc.FEATURE_R_RNA,
        "start": 1, "stop": 120, "strand": "+",
        "locus": "TEST_00001", "gene": None, "product": "5.8S ribosomal RNA",
        "db_xrefs": ["SO:0000252", None],
    }]}
    out = tmp_path / "features.tsv"
    write_features(sequences, features_by_sequence, out)   # used to raise TypeError
    body = [l for l in out.read_text().splitlines() if not l.startswith("#")]
    assert "SO:0000252" in body[0]
    assert "None" not in body[0]
