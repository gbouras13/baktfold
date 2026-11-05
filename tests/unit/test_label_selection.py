"""Tests for hierarchy label selection in ``combine_annotation``.

Baktfold picks one product per CDS from its per-database top hits, in the order
custom > Swiss-Prot > AFDB > PDB > CATH. It used to take the top-ranked hit
unconditionally and only afterwards rename "Uncharacterized protein ..." to
"hypothetical protein". That:

  1. discarded informative hits from lower-precedence databases,
  2. left the CDS *unmarked* as hypothetical, so ``.summary.txt`` counted it
     as "annotated with Baktfold function",
  3. still gave it an ``annotation_confidence`` tier, and
  4. was case-sensitive, so "Putative uncharacterized protein" passed as a
     function.

Seen in the ``run_all`` golden (MEGJMN_036/057/058/059/070).

Run::

    pytest tests/unit/test_label_selection.py -v
"""
from __future__ import annotations

import pytest

import baktfold.bakta.annotation as anno


def _hit(source, description, **extra):
    hit = {
        "source": source, "id": f"{source}_X", "description": description,
        "query_cov": 0.95, "subject_cov": 0.95, "identity": 0.5, "evalue": 1e-30,
    }
    hit.update(extra)
    return hit


def _cds(pstc, product="hypothetical protein", hypothetical=True):
    feat = {"type": "cds", "locus": "c1", "product": product, "db_xrefs": []}
    if hypothetical:
        feat["hypothetical"] = True
    feat["pstc"] = pstc
    return feat


@pytest.mark.parametrize(
    "uninformative",
    [
        "Uncharacterized protein YqeH",
        "uncharacterized protein",
        "Putative uncharacterized protein",
        "Uncharacterised protein",
        "hypothetical protein",
        "Hypothetical protein",
        "Conserved hypothetical protein",
        "Uncharacterized conserved protein",
        "Uncharacterized oxidoreductase YdgJ",
        "Insertion element IS900 uncharacterized 42 kDa protein",
    ],
)
def test_uninformative_top_hit_falls_through_to_next_database(uninformative):
    feat = _cds([
        _hit("swissprot", uninformative),
        _hit("afdb", "IS66 family transposase"),
    ])
    anno.combine_annotation(feat, fast=False)

    assert feat["product"] == "IS66 family transposase"
    assert "hypothetical" not in feat
    assert feat["baktfold"] is True
    assert feat["annotation_confidence"] in {"high", "medium", "low"}


def test_custom_db_uninformative_does_not_mask_swissprot():
    feat = _cds([
        _hit("custom_db", "Uncharacterized protein"),
        _hit("swissprot", "DNA gyrase subunit A"),
    ])
    anno.combine_annotation(feat, fast=False)
    assert feat["product"] == "DNA gyrase subunit A"


def test_only_uninformative_hits_stay_hypothetical_without_confidence():
    feat = _cds([
        _hit("swissprot", "Uncharacterized protein YqeH"),
        _hit("afdb", "hypothetical protein"),
    ])
    anno.combine_annotation(feat, fast=False)

    assert feat["product"] == "hypothetical protein"
    assert feat["hypothetical"] is True
    # a hit was found, but it did not supply a label
    assert feat["baktfold"] is True
    assert "annotation_confidence" not in feat
    # cross-references to the hits are still kept
    assert len(feat["db_xrefs"]) == 2


def test_confidence_comes_from_the_hit_that_supplied_the_label():
    """The masked Swiss-Prot hit is a perfect match; the AFDB hit that actually
    supplies the label is weak. The tier must describe the AFDB hit."""
    feat = _cds([
        _hit("swissprot", "Uncharacterized protein"),
        _hit("afdb", "AAA+ ATPase domain-containing protein",
             query_cov=0.3, subject_cov=0.3, identity=0.1, evalue=1e-4),
    ])
    feat["prostt5_confidence"] = 30.0
    anno.combine_annotation(feat, fast=False)

    assert feat["product"] == "AAA+ ATPase domain-containing protein"
    assert feat["annotation_confidence"] == "low"


def test_fast_mode_still_skips_afdb():
    feat = _cds([
        _hit("swissprot", "Uncharacterized protein"),
        _hit("afdb", "IS66 family transposase"),
        _hit("pdb", "Transposase"),
    ])
    anno.combine_annotation(feat, fast=True)
    assert feat["product"] == "Transposase"


def test_informative_top_hit_unchanged():
    """Ordinary case: precedence is untouched when the top hit is informative."""
    feat = _cds([
        _hit("swissprot", "Replication protein RepA"),
        _hit("afdb", "Something else"),
    ])
    anno.combine_annotation(feat, fast=False)
    assert feat["product"] == "Replication protein RepA"
    assert "hypothetical" not in feat


def test_existing_product_not_downgraded_by_uninformative_hit():
    """--all-proteins: a curated input product must not be replaced by a hit
    that only says 'hypothetical/uncharacterized protein'."""
    feat = _cds(
        [_hit("swissprot", "Uncharacterized protein")],
        product="DNA gyrase subunit A",
        hypothetical=False,
    )
    anno.combine_annotation(feat, fast=False)
    assert feat["product"] == "DNA gyrase subunit A"
    assert "hypothetical" not in feat


def test_no_hits_stays_hypothetical():
    feat = _cds(None)
    anno.combine_annotation(feat, fast=False)
    assert feat["product"] == "hypothetical protein"
    assert feat["hypothetical"] is True
    assert "baktfold" not in feat


@pytest.mark.parametrize(
    "product, expected",
    [
        ("DUF1234 domain-containing protein", True),
        ("Uncharacterized protein YqeH", False),
        ("Putative uncharacterized protein", False),
        ("hypothetical protein", False),
        ("", False),
        (None, False),
        ("IS66 family transposase", True),
        ("UPF0234 protein YajQ", True),
        ("Uncharacterized oxidoreductase YdgJ", False),
        ("Uncharacterized ABC transporter ATP-binding protein YbbA", False),
        ("Conserved hypothetical protein", False),
        ("Hypothetically-named protein", True),  # word boundary: not the word itself
    ],
)
def test_is_informative_product(product, expected):
    assert anno.is_informative_product(product) is expected


def test_suspect_product_falls_through():
    """revise_cds_product blanks contig-name-like products; those must not win."""
    feat = _cds([
        _hit("swissprot", "NODE_12_length_5000 protein"),
        _hit("afdb", "Transposase"),
    ])
    anno.combine_annotation(feat, fast=False)
    assert feat["product"] == "Transposase"


def test_input_product_kept_even_if_it_mentions_uncharacterized():
    """--all-proteins: the strict rule applies to Foldseek hits, not to the
    product the input annotation tool assigned."""
    feat = _cds(
        [_hit("swissprot", "hypothetical protein")],
        product="Uncharacterized oxidoreductase YdgJ",
        hypothetical=False,
    )
    anno.combine_annotation(feat, fast=False)
    assert feat["product"] == "Uncharacterized oxidoreductase YdgJ"
    assert "hypothetical" not in feat


def test_input_uncharacterized_protein_renamed_is_flagged():
    """Audit N7: with --all-proteins an input "Uncharacterized protein X" with
    no usable hit is renamed by Bakta's cleanup; it must then be flagged."""
    feat = _cds([], product="Uncharacterized protein YdgJ", hypothetical=False)
    anno.combine_annotation(feat, fast=False)
    assert feat["product"] == "hypothetical protein"
    assert feat["hypothetical"] is True
