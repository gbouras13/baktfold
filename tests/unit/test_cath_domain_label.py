"""CATH hits are labelled at domain level (audit A5, option 3).

A CATH superfamily match is evidence that a protein contains a domain of that
superfamily, not that the whole protein has the named function. Previously the
CATH name was copied into /product verbatim (e.g. "Vaccinia Virus protein VP39"
on a methyltransferase) - reviewer concern R2 #6, fold-level similarity
presented as protein function.

Run::

    pytest tests/unit/test_cath_domain_label.py -v
"""
from __future__ import annotations

import pytest

import baktfold.bakta.annotation as anno


@pytest.mark.parametrize(
    "name, expected",
    [
        ("SH2 domain", "SH2 domain-containing protein"),
        ("Tetratricopeptide repeat domain", "Tetratricopeptide repeat domain-containing protein"),
        ("Enoyl acyl carrier protein reductase", "Enoyl acyl carrier protein reductase domain-containing protein"),
        ("Osmoprotection protein (prox); domain 2", "Osmoprotection protein (prox) domain-containing protein"),
        ("Potassium Channel Kv1.1; Chain A", "Potassium Channel Kv1.1 domain-containing protein"),
        ("outer membrane protein from Thermus thermophilus HB27.", "outer membrane protein domain-containing protein"),
        ("AF1437-like", "AF1437-like domain-containing protein"),
        ("Death domains", "Death domain-containing protein"),
    ],
)
def test_cath_domain_label(name, expected):
    assert anno.cath_domain_label(name) == expected


def _hit(source, description):
    return {"source": source, "id": f"{source}_X", "description": description,
            "query_cov": 0.9, "subject_cov": 0.9, "identity": 0.2, "evalue": 1e-12}


def test_cath_only_hit_gets_domain_level_product():
    feat = {"type": "cds", "locus": "c1", "product": "hypothetical protein",
            "hypothetical": True, "db_xrefs": [],
            "pstc": [_hit("cath", "Vaccinia Virus protein VP39")]}
    anno.combine_annotation(feat, fast=False)
    assert feat["product"] == "Vaccinia Virus protein VP39 domain-containing protein"
    assert "hypothetical" not in feat


def test_higher_precedence_hit_unaffected():
    feat = {"type": "cds", "locus": "c1", "product": "hypothetical protein",
            "hypothetical": True, "db_xrefs": [],
            "pstc": [_hit("swissprot", "mRNA cap 0 methyltransferase"),
                     _hit("cath", "Vaccinia Virus protein VP39")]}
    anno.combine_annotation(feat, fast=False)
    assert feat["product"] == "mRNA cap 0 methyltransferase"


def test_uninformative_cath_name_still_skipped():
    feat = {"type": "cds", "locus": "c1", "product": "hypothetical protein",
            "hypothetical": True, "db_xrefs": [],
            "pstc": [_hit("cath", "Hypothetical PUA domain-like; domain 1")]}
    anno.combine_annotation(feat, fast=False)
    assert feat["product"] == "hypothetical protein"
    assert feat["hypothetical"] is True
