"""
End-to-end *golden* regression tests for baktfold's output writers, driven
through ``baktfold json``.

``baktfold json`` is a pure, deterministic, dependency-free transform of a JSON
into every standard output format, exercising the exact same writers used by
``run``/``compare`` (gff/gbff/embl/tsv/inference.tsv/faa/ffn/summary). These
tests pin the *bytes* of those outputs against committed snapshots so any drift
— e.g. the INSDC ``/transl_table`` and tRNA-inference regressions — lights up in
CI without needing a database, ProstT5, Foldseek or a GPU.

The annotation-bearing portion of each format is snapshotted; the raw nucleotide
bulk (GenBank ORIGIN, EMBL SQ, GFF ##FASTA, .fna/.ffn bodies) is stripped so the
goldens stay small and diffs stay readable. Volatile version strings in the
headers are normalised.

Snapshots live under tests/unit/snapshots/test_json_golden/. Bootstrap/refresh
them (after an *intended* output change) with::

    pytest tests/unit/test_json_golden.py --snapshot-update

Run::

    pytest tests/unit/test_json_golden.py -v
"""
from __future__ import annotations

import re
from pathlib import Path

import pytest
from click.testing import CliRunner

from baktfold import main_cli

REPO_ROOT = Path(__file__).resolve().parents[2]
TEST_DATA = REPO_ROOT / "tests" / "test_data"
SUBDIR = "test_json_golden"

# Small, real, tracked fixtures (no DB/foldseek needed to reconstruct from them).
GENOME_JSON = TEST_DATA / "assembly_bakta_output" / "assembly.json"
PROTEINS_JSON = TEST_DATA / "assembly_bakta_proteins_output" / "assembly.hypotheticals.json"

# baktfold version is the only volatile content in these formats (the GenBank
# LOCUS date is the fixed 01-JAN-1980 placeholder).
_VERSION_RE = re.compile(r"((?:Software|Database): v)\S+")


def _normalize_version(text: str) -> str:
    return _VERSION_RE.sub(r"\1<VERSION>", text)


def _strip_after(text: str, marker: str) -> str:
    """Drop everything from the first occurrence of ``marker`` onward."""
    idx = text.find(marker)
    return text if idx == -1 else text[:idx] + "\n"


def _headers_only(text: str) -> str:
    """Keep only FASTA ``>`` header lines (feature identity/order/products)."""
    return "\n".join(l for l in text.splitlines() if l.startswith(">")) + "\n"


# ext -> function turning the raw file body into the snapshotted text
_GENOME_PROCESSORS = {
    "gff3": lambda t: _normalize_version(_strip_after(t, "\n##FASTA")),
    "gbff": lambda t: _normalize_version(_strip_after(t, "\nORIGIN")),
    "embl": lambda t: _normalize_version(_strip_after(t, "\nSQ   ")),
    "tsv": _normalize_version,
    "inference.tsv": _normalize_version,
    "faa": _normalize_version,
    "ffn": _headers_only,
    "summary.txt": _normalize_version,
}
_PROTEINS_PROCESSORS = {
    "tsv": _normalize_version,
    "faa": _normalize_version,
    "summary.txt": _normalize_version,
}


def _reconstruct(out_dir: Path, json_path: Path, prefix: str) -> Path:
    result = CliRunner().invoke(
        main_cli,
        ["json", "-i", str(json_path), "-o", str(out_dir), "-p", prefix, "-f"],
        catch_exceptions=False,
    )
    assert result.exit_code == 0, f"baktfold json failed:\n{result.output}"
    return out_dir


@pytest.fixture(scope="module")
def genome_out(tmp_path_factory):
    if not GENOME_JSON.exists():
        pytest.skip("genome fixture JSON not present")
    return _reconstruct(tmp_path_factory.mktemp("golden_genome"), GENOME_JSON, "g")


@pytest.fixture(scope="module")
def proteins_out(tmp_path_factory):
    if not PROTEINS_JSON.exists():
        pytest.skip("proteins fixture JSON not present")
    return _reconstruct(tmp_path_factory.mktemp("golden_proteins"), PROTEINS_JSON, "p")


@pytest.mark.parametrize("ext", sorted(_GENOME_PROCESSORS))
def test_genome_output_golden(ext, genome_out, text_snapshot):
    raw = (genome_out / f"g.{ext}").read_text()
    text_snapshot(_GENOME_PROCESSORS[ext](raw), f"genome.{ext}", subdir=SUBDIR)


@pytest.mark.parametrize("ext", sorted(_PROTEINS_PROCESSORS))
def test_proteins_output_golden(ext, proteins_out, text_snapshot):
    raw = (proteins_out / f"p.{ext}").read_text()
    text_snapshot(_PROTEINS_PROCESSORS[ext](raw), f"proteins.{ext}", subdir=SUBDIR)
