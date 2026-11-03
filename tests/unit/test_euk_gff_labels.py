"""Eukaryotic outputs must carry Baktfold's CDS label.

Baktfold annotates CDS features. The eukaryotic GFF3 writer emitted CDS lines
with only ID/Parent, and the mRNA line kept the input product (usually
"hypothetical protein"), so every Baktfold label was lost from eukaryotic GFF3
(and the mRNA product stayed hypothetical in GBFF/TSV too). The mRNA ``Note``
was also filled from db_xrefs, and mRNA attributes were not GFF3-escaped.

Run::

    pytest tests/unit/test_euk_gff_labels.py -v
"""
from __future__ import annotations

import io as _io

from baktfold.io import gff
from baktfold.io.io import propagate_cds_products_to_mrna


def _mrna(locus, start, stop, product="hypothetical protein", **kw):
    f = {"type": "mRNA", "sequence": "chr1", "start": start, "stop": stop,
         "strand": "+", "locus": locus, "product": product,
         "db_xrefs": ["SO:0000234"], "starts": [start], "stops": [stop]}
    f.update(kw)
    return f


def _cds(locus, start, stop, product):
    return {"type": "cds", "sequence": "chr1", "start": start, "stop": stop,
            "strand": "+", "locus": locus, "product": product,
            "db_xrefs": ["afdb_v6:swissprot_P12345", "SO:0001217"],
            "starts": [start, start + 30], "stops": [start + 20, stop]}


def test_cds_product_propagates_to_mrna():
    m, c = _mrna("g1", 100, 400), _cds("g1", 110, 390, "DNA ligase")
    assert propagate_cds_products_to_mrna([m, c]) == 1
    assert m["product"] == "DNA ligase"


def test_duplicate_locus_transcripts_pair_in_order():
    """Two transcripts under one locus tag each get their own CDS's label."""
    m1, m2 = _mrna("g1", 100, 400), _mrna("g1", 100, 900)
    c1, c2 = _cds("g1", 110, 390, "Kinase A"), _cds("g1", 500, 880, "Kinase B")
    propagate_cds_products_to_mrna([m1, c1, m2, c2])
    assert (m1["product"], m2["product"]) == ("Kinase A", "Kinase B")


def test_unmatched_mrna_untouched():
    m = _mrna("g2", 100, 400)
    propagate_cds_products_to_mrna([m, _cds("g1", 110, 390, "DNA ligase")])
    assert m["product"] == "hypothetical protein"


def test_euk_cds_lines_carry_label_and_dbxref():
    fh = _io.StringIO()
    gff.write_euk_cds_feature(fh, "chr1", _cds("g1", 110, 390, "DNA ligase; ATP-dependent"))
    lines = fh.getvalue().strip().split("\n")
    assert len(lines) == 2  # one line per CDS part
    for line in lines:
        attrs = line.split("\t")[8]
        assert "product=DNA ligase%3B ATP-dependent" in attrs  # escaped
        assert "Dbxref=afdb_v6:swissprot_P12345,SO:0001217" in attrs
        assert "locus_tag=g1" in attrs


def test_mrna_note_and_escaping():
    fh = _io.StringIO()
    m = _mrna("g1", 100, 400, product="Ser/Thr kinase, putative", note="from EggNOG")
    gff.write_mrna_feature(fh, "chr1", m)
    mrna_line = fh.getvalue().split("\n")[0]
    attrs = mrna_line.split("\t")[8]
    assert "product=Ser/Thr kinase%2C putative" in attrs
    assert "Note=from EggNOG" in attrs
    assert "Dbxref=SO:0000234" in attrs
