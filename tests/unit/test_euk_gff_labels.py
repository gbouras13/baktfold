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


def _mrna(locus, start, stop, product="hypothetical protein", parts=None, **kw):
    parts = parts or [(start, stop)]
    f = {"type": "mRNA", "sequence": "chr1", "start": start, "stop": stop,
         "strand": "+", "locus": locus, "product": product,
         "db_xrefs": ["SO:0000234"], "starts": [p[0] for p in parts], "stops": [p[1] for p in parts]}
    f.update(kw)
    return f


def _cds(locus, start, stop, product, parts=None):
    parts = parts or [(start, start + 20), (start + 30, stop)]
    return {"type": "cds", "sequence": "chr1", "start": start, "stop": stop,
            "strand": "+", "locus": locus, "product": product,
            "db_xrefs": ["afdb_v6:swissprot_P12345", "SO:0001217"],
            "starts": [p[0] for p in parts], "stops": [p[1] for p in parts]}


def test_cds_product_propagates_to_mrna():
    # mRNA exons 100-130 / 140-400 share the CDS's intron 130-140
    m = _mrna("g1", 100, 400, parts=[(100, 130), (140, 400)])
    c = _cds("g1", 110, 390, "DNA ligase")
    assert propagate_cds_products_to_mrna([m, c]) == 1
    assert m["product"] == "DNA ligase"


def test_alternative_transcripts_paired_by_splice_structure():
    """Shape of a real Aspergillus gene (I1I15_012813): two transcripts with
    the same span, three vs two exons; each CDS must label its own mRNA,
    whatever the feature order."""
    m3 = _mrna("g", 100, 1000, parts=[(100, 300), (400, 600), (700, 1000)])
    m2 = _mrna("g", 100, 1000, parts=[(100, 600), (700, 1000)])
    c3 = _cds("g", 150, 900, "Kinase three-exon", parts=[(150, 300), (400, 600), (700, 900)])
    c2 = _cds("g", 150, 900, "Kinase two-exon", parts=[(150, 600), (700, 900)])
    for order in ([m2, m3, c3, c2], [m3, c2, m2, c3]):
        for f in (m2, m3):
            f["product"] = "hypothetical protein"
        propagate_cds_products_to_mrna(order)
        assert (m3["product"], m2["product"]) == ("Kinase three-exon", "Kinase two-exon")


def test_nested_single_exon_transcripts_not_swapped():
    """Audit N2: the old 'first CDS inside' rule gave mRNA 100-2000 the other
    transcript's label. Both CDS fit both mRNAs, with different products:
    ambiguous, so neither mRNA may be given a possibly wrong label."""
    a = _mrna("g", 100, 2000)
    b = _mrna("g", 200, 2000)
    ca = _cds("g", 400, 1900, "Kinase A", parts=[(400, 1900)])
    cb = _cds("g", 250, 1900, "hypothetical protein", parts=[(250, 1900)])
    propagate_cds_products_to_mrna([a, b, cb, ca])
    # b (200-2000) can hold either CDS; a (100-2000) can hold either too:
    # no unique pairing exists, so both keep their own product
    assert a["product"] == "hypothetical protein"
    assert b["product"] == "hypothetical protein"


def test_ambiguous_but_same_product_is_assigned():
    a = _mrna("g", 100, 2000)
    b = _mrna("g", 200, 2000)
    propagate_cds_products_to_mrna([a, b,
        _cds("g", 400, 1900, "Kinase A", parts=[(400, 1900)]),
        _cds("g", 250, 1900, "Kinase A", parts=[(250, 1900)])])
    assert a["product"] == b["product"] == "Kinase A"


def test_cds_with_intron_not_in_mrna_is_not_paired():
    m = _mrna("g", 100, 1000)  # unspliced
    c = _cds("g", 150, 900, "Kinase", parts=[(150, 400), (500, 900)])
    propagate_cds_products_to_mrna([m, c])
    assert m["product"] == "hypothetical protein"


def test_unmatched_mrna_untouched():
    m = _mrna("g2", 100, 400, parts=[(100, 130), (140, 400)])
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
