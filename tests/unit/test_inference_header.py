"""Tests that the inference TSV header names the values actually written.

``build_inference_header`` and ``map_aa_columns`` have to stay in lockstep.
They drifted: the duplicate-locus header *substituted* 'ID' for 'Length' while
the row *appended* 'ID' and went on writing the length. The row therefore
carried one more field than the header named, and every column from Product
rightwards was labelled with the value to its left::

    Locus            ID              Product  Annotation_Confidence  Swissprot ...
    OYNLMN_1920_2375 KFOAEIDSCJ_41   152      Progestin and adipo... swissprot_Q6TCH7

A ``csv.DictReader`` on that returns "152" for Product — non-empty and not
"hypothetical protein", so every CDS scored as functionally annotated. Across
802 protist genomes this read as 91.2% functionally annotated against a true
~70%, with 121 genomes apparently at 100%. The annotations themselves were
fine; only the header was wrong, so the numbers derived from it were not.

Eukaryotic inputs are the affected ones — duplicate locus tags are common in
Ensembl eukaryotic GenBank, and ``has_duplicate_locus`` is what turns this on.

The invariant pinned here is cheap and would have caught it: for every
combination of flags, the header must have exactly as many names as the row
has values, and each name must sit over its own value.

Run::

    pytest tests/unit/test_inference_header.py -v
"""
from __future__ import annotations

import itertools

import pytest

from baktfold.io.tsv import build_inference_header, map_aa_columns


FLAG_MATRIX = list(itertools.product([False, True], repeat=4))
FLAG_IDS = [
    f"dup={d}-custom={c}-fast={f}-struct={s}" for d, c, f, s in FLAG_MATRIX
]


def _feature():
    """A CDS feature carrying every optional field the row can emit."""
    return {
        'locus': 'OYNLMN_1920_2375',
        'id': 'KFOAEIDSCJ_41',
        'length': 152,
        'product': 'Progestin and adipoQ receptor family member 3',
        'annotation_confidence': 'high',
        'db_xrefs': [
            'swissprot_Q6TCH7',
            'afdbclusters_A0AA40IA44',
            'pdb:1abc',
            'cath:af_Q9VFG7_111_330_1.20',
            'custom:custom_thing',
        ],
        'tmscore': 0.81,
        'lddt': 0.77,
    }


def _header_and_row(dup, custom, fast, struct, first_column='Locus'):
    header = build_inference_header(first_column, custom, dup, fast, struct)
    row = map_aa_columns(_feature(), custom, dup, fast, struct)
    return header, row


# ─── the invariant ──────────────────────────────────────────────────────────


@pytest.mark.parametrize("dup,custom,fast,struct", FLAG_MATRIX, ids=FLAG_IDS)
def test_header_has_one_name_per_value(dup, custom, fast, struct):
    header, row = _header_and_row(dup, custom, fast, struct)
    assert len(header) == len(row), (
        f"header names {len(header)} columns but the row writes {len(row)} "
        f"values\n  header: {header}\n  row:    {row}"
    )


@pytest.mark.parametrize("dup,custom,fast,struct", FLAG_MATRIX, ids=FLAG_IDS)
def test_each_name_sits_over_its_own_value(dup, custom, fast, struct):
    """The count matching is necessary but not sufficient — check the mapping."""
    header, row = _header_and_row(dup, custom, fast, struct)
    parsed = dict(zip(header, row))
    feat = _feature()

    assert parsed['Locus'] == feat['locus']
    assert parsed['Length'] == str(feat['length'])
    assert parsed['Product'] == feat['product']
    assert parsed['Annotation_Confidence'] == feat['annotation_confidence']
    if dup:
        assert parsed['ID'] == feat['id']
    else:
        assert 'ID' not in header


# ─── the specific corruption ────────────────────────────────────────────────


def test_duplicate_locus_product_is_not_the_length():
    """The exact misread: Product returning "152" made every CDS in an affected
    eukaryotic genome score as functionally annotated."""
    header, row = _header_and_row(dup=True, custom=False, fast=False, struct=False)
    parsed = dict(zip(header, row))

    assert parsed['Product'] != '152'
    assert parsed['Product'] == 'Progestin and adipoQ receptor family member 3'


def test_duplicate_locus_adds_a_column_rather_than_replacing_one():
    """'ID' must be inserted, not substituted for 'Length'."""
    plain = build_inference_header('Locus', False, False, False, False)
    dup = build_inference_header('Locus', False, True, False, False)

    assert len(dup) == len(plain) + 1
    assert set(plain) < set(dup)          # nothing was displaced
    assert dup.index('ID') == 1           # straight after Locus, as the row emits it
    assert 'Length' in dup


def test_no_value_falls_off_the_end():
    """CATH used to be written with no header name at all."""
    header, row = _header_and_row(dup=True, custom=True, fast=False, struct=True)
    assert header[-1] == 'LDDT'
    assert len(row) == len(header)
    assert dict(zip(header, row))['CATH'] == 'af_Q9VFG7_111_330_1.20'


# ─── column order is part of the contract ───────────────────────────────────


def test_default_header_order_is_unchanged():
    """The ordinary (non-eukaryotic) header is what every existing parser and
    the golden snapshots expect — the fix must not move it."""
    assert build_inference_header('Locus', False, False, False, False) == [
        'Locus', 'Length', 'Product', 'Annotation_Confidence',
        'Swissprot', 'AFDBClusters', 'PDB', 'CATH',
    ]


@pytest.mark.parametrize("dup,custom,fast,struct", FLAG_MATRIX, ids=FLAG_IDS)
def test_annotation_confidence_follows_product(dup, custom, fast, struct):
    """It is positioned relative to Product; a wrong header put it elsewhere."""
    header = build_inference_header('Locus', custom, dup, fast, struct)
    assert header[header.index('Product') + 1] == 'Annotation_Confidence'


def test_fast_drops_only_afdbclusters():
    full = build_inference_header('Locus', False, False, False, False)
    fast = build_inference_header('Locus', False, False, True, False)
    assert set(full) - set(fast) == {'AFDBClusters'}


# ─── proteins mode ──────────────────────────────────────────────────────────


@pytest.mark.parametrize("custom,fast,struct", list(itertools.product([False, True], repeat=3)))
def test_proteins_mode_header_matches_its_row(custom, fast, struct):
    """proteins mode always writes has_duplicate_locus=False, so it never had
    the substitution bug — pin that it stays that way."""
    header, row = _header_and_row(False, custom, fast, struct, first_column='ID')
    assert header[0] == 'ID'
    assert len(header) == len(row)
