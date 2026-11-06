"""Repeated sequence IDs (RERUN_INVENTORY.md section 5).

Ensembl's Emiliania huxleyi flatfile truncates ACCESSION, so 820 of its 7,794
records share a Biopython ``record.id`` with another (EMIHUscaffold_10000,
_10004, ... all read '...:EMIHUscaffold_1000'). convert-euk keyed sequences and
features on that ID, so:

  * features were attached to the wrong sequence, and
  * every writer that emits each sequence's features once per sequence entry
    (run/compare/json flattening, GFF3, GBFF/EMBL, TSV) repeated them once
    per entry: 38,552 CDS became 42,179 after one run, and 66,763 (k^2)
    after a second run seeded with that output.

convert-euk now makes record IDs unique; the writers emit each feature once
even for an older JSON that still repeats IDs.

Run::

    pytest tests/unit/test_duplicate_sequence_ids.py -v
"""
from __future__ import annotations

import json as _json
from pathlib import Path

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from click.testing import CliRunner

from baktfold import main_cli
from baktfold.io.eukaryotic_to_json import eukaryotic_gbk_to_json, make_record_ids_unique
from baktfold.utils.util import iter_sequences_once

FIXTURE = Path(__file__).parent / "fixtures" / "euk_partial_cds.gbk"
GENOME_JSON = Path(__file__).resolve().parents[2] / "tests" / "test_data" / "GCF_002368115.json"


def _rec(id_, name):
    return SeqRecord(Seq("ACGT"), id=id_, name=name)


# ── convert-euk ──────────────────────────────────────────────────────────────


def test_unique_ids_are_left_alone():
    recs = [_rec("a", "A"), _rec("b", "B")]
    assert make_record_ids_unique(recs) == 0
    assert [r.id for r in recs] == ["a", "b"]


def test_every_record_sharing_an_id_takes_its_locus_name():
    recs = [_rec("acc:scaf_1000", "scaf_10000"), _rec("acc:scaf_1000", "scaf_1000"),
            _rec("other", "other_locus")]
    assert make_record_ids_unique(recs) == 2
    assert [r.id for r in recs] == ["scaf_10000", "scaf_1000", "other"]


def test_falls_back_to_a_suffix_when_names_collide_too():
    recs = [_rec("x", "dup"), _rec("x", "dup"), _rec("dup_name_taken", "y")]
    make_record_ids_unique(recs)
    ids = [r.id for r in recs]
    assert len(set(ids)) == 3 and ids[0].startswith("x_") and ids[1].startswith("x_")


def test_convert_euk_keeps_features_on_their_own_sequence(tmp_path):
    records = list(SeqIO.parse(FIXTURE, "genbank"))
    assert len(records) >= 2
    owner = {}  # CDS locus tag -> LOCUS name of the record it is on
    for rec in records:
        for f in rec.features:
            if f.type == "CDS" and "locus_tag" in f.qualifiers:
                owner[f.qualifiers["locus_tag"][0]] = rec.name
    for rec in records[1:]:
        rec.id = records[0].id  # every record now shares one ID

    out = tmp_path / "out.json"
    eukaryotic_gbk_to_json(records, out, verbose=False)
    data = _json.loads(out.read_text())

    seq_ids = [s["id"] for s in data["sequences"]]
    assert len(set(seq_ids)) == len(seq_ids)
    assert set(seq_ids) == {r.name for r in records}
    cds = [f for f in data["features"] if f["type"] == "cds"]
    assert cds and all(f["sequence"] == owner[f["locus"]] for f in cds if f["locus"] in owner)


# ── writers, for JSONs that still repeat IDs ─────────────────────────────────


def test_iter_sequences_once_flags_repeats():
    seqs = [{"id": "a"}, {"id": "b"}, {"id": "a"}]
    assert [first for _, first in iter_sequences_once(seqs)] == [True, True, False]


@pytest.fixture(scope="module")
def reconstructed(tmp_path_factory):
    """`baktfold json` on the Bakta fixture as is, and with its first sequence
    entry repeated (same ID), as in an E. huxleyi JSON written by <=0.4.0."""
    if not GENOME_JSON.exists():
        pytest.skip("genome fixture JSON not present")
    tmp = tmp_path_factory.mktemp("dup_seq_ids")
    data = _json.loads(GENOME_JSON.read_text())
    dup = dict(data)
    dup["sequences"] = data["sequences"] + [dict(data["sequences"][0])]
    (tmp / "dup.json").write_text(_json.dumps(dup))
    for name, src in (("plain", GENOME_JSON), ("dup", tmp / "dup.json")):
        r = CliRunner().invoke(main_cli, ["json", "-i", str(src), "-o", str(tmp / name), "-p", "g", "-f"],
                               catch_exceptions=False)
        assert r.exit_code == 0, r.output
    return tmp / "plain", tmp / "dup"


def _count_cds(out: Path):
    j = sum(1 for f in _json.loads((out / "g.json").read_text())["features"] if f["type"] == "cds")
    gff = sum(1 for l in (out / "g.gff3").read_text().splitlines()
              if not l.startswith("#") and l.count("\t") >= 8 and l.split("\t")[2] == "CDS")
    gbff = sum(1 for l in (out / "g.gbff").read_text().splitlines() if l.startswith("     CDS "))
    tsv = sum(1 for l in (out / "g.tsv").read_text().splitlines() if "\tcds\t" in l)
    return j, gff, gbff, tsv


def test_repeated_sequence_id_does_not_repeat_features(reconstructed):
    plain, dup = reconstructed
    assert _count_cds(dup) == _count_cds(plain)


def test_repeated_sequence_is_still_written(reconstructed):
    plain, dup = reconstructed
    n = lambda p: (p / "g.fna").read_text().count(">")
    assert n(dup) == n(plain) + 1
