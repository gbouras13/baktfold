"""Eukaryotic CDS: partial ends, /codon_start and GFF3 phase (audit E2/E3).

Fixture: three Entamoeba contigs from tests/test_data/protist.gbk.gz with
partial CDS, including EHI_065120
``complement(join(<1..47,140..329,402..>1012))`` with ``/codon_start=2``.

Before the fix:
  * GBFF/EMBL wrote ``/codon_start=1`` for every CDS and dropped the ``<``/``>``
    markers, so the location no longer matched ``/translation``;
  * GFF3 phase was ``offset % 3`` (bases consumed, not bases to skip), the
    minus strand was reversed although the parts are already stored 5'->3',
    and ``codon_start`` was ignored.

Run::

    pytest tests/unit/test_euk_partial_cds.py -v
"""
from __future__ import annotations

import io as _io
import json
import shutil
import subprocess
import sys
from pathlib import Path

import pytest
from Bio import SeqIO

from baktfold.io import gff

FIXTURE = Path(__file__).parent / "fixtures" / "euk_partial_cds.gbk"


def _baktfold():
    exe = Path(sys.executable).with_name("baktfold")
    exe = str(exe) if exe.exists() else shutil.which("baktfold")
    if exe is None:
        pytest.skip("baktfold console script not installed")
    return exe


@pytest.fixture(scope="module")
def regenerated(tmp_path_factory):
    """convert-euk the fixture, then `baktfold json` all outputs from it."""
    tmp = tmp_path_factory.mktemp("euk_partial")
    exe = _baktfold()
    subprocess.run([exe, "convert-euk", "-i", str(FIXTURE), "-o", str(tmp / "in.json")],
                   check=True, capture_output=True)
    data = json.loads((tmp / "in.json").read_text())
    data.setdefault("baktfold_run", {})["euk"] = True
    (tmp / "in_euk.json").write_text(json.dumps(data))
    subprocess.run([exe, "json", "-i", str(tmp / "in_euk.json"), "-o", str(tmp / "out")],
                   check=True, capture_output=True)
    return data, tmp / "out"


def _input_cds():
    out = {}
    for rec in SeqIO.parse(FIXTURE, "genbank"):
        for f in rec.features:
            if f.type == "CDS":
                out[f.qualifiers["locus_tag"][0]] = (str(f.location), f.qualifiers.get("codon_start", ["1"])[0])
    return out


def test_converter_keeps_codon_start_and_partial_flags(regenerated):
    data, _ = regenerated
    cds = {f["locus"]: f for f in data["features"] if f["type"] == "cds"}
    f = cds["EHI_065120"]
    assert f["codon_start"] == 2
    assert f["partial_start"] is True and f["partial_stop"] is True
    assert cds["EHI_034110"].get("partial_stop") is True
    assert "partial_start" not in cds["EHI_034110"]


def test_gbff_location_and_codon_start_round_trip(regenerated):
    _, out = regenerated
    expected = _input_cds()
    seen = 0
    for rec in SeqIO.parse(out / "baktfold.gbff", "genbank"):
        for f in rec.features:
            if f.type == "CDS":
                lt = f.qualifiers["locus_tag"][0]
                assert (str(f.location), f.qualifiers["codon_start"][0]) == expected[lt], lt
                seen += 1
    assert seen == len(expected) == 3


def test_gff3_phases(regenerated):
    _, out = regenerated
    rows = [l.split("\t") for l in (out / "baktfold.gff3").read_text().split("##FASTA")[0].splitlines()
            if not l.startswith("#") and l.count("\t") == 8 and l.split("\t")[2] == "CDS"]
    ehi = [(int(r[3]), int(r[4]), r[7]) for r in rows if "locus_tag=EHI_065120" in r[8]]
    # 5'->3' on the minus strand: 402..1012 first, then 140..329, then 1..47.
    # codon_start=2 -> phase 1 on the first part; 611 - 1 = 610 coding bases
    # -> 610 % 3 = 1 -> phase 2; 610 + 190 = 800 -> 800 % 3 = 2 -> phase 1.
    assert ehi == [(402, 1012, "1"), (140, 329, "2"), (1, 47, "1")]
    assert all("partial=true" in r[8] for r in rows if "EHI_065120" in r[8])


def test_phase_minus_strand_complete_cds():
    """Synthetic 3-exon minus-strand CDS, stored 5'->3' as the converter does."""
    feat = {"locus": "g1", "strand": "-", "starts": [700, 400, 100], "stops": [799, 509, 151],
            "start": 100, "stop": 799, "product": "x"}
    fh = _io.StringIO()
    gff.write_euk_cds_feature(fh, "chr", feat)
    phases = [(int(l.split("\t")[3]), int(l.split("\t")[7])) for l in fh.getvalue().strip().split("\n")]
    # lengths 100, 110, 52: phase 0; 100 % 3 = 1 -> 2; 210 % 3 = 0 -> 0
    assert phases == [(700, 0), (400, 2), (100, 0)]
