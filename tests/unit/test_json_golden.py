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
from baktfold.utils.util import get_version

REPO_ROOT = Path(__file__).resolve().parents[2]
TEST_DATA = REPO_ROOT / "tests" / "test_data"
SUBDIR = "test_json_golden"

# Small, real, tracked fixtures (no DB/foldseek needed to reconstruct from them).
GENOME_JSON = TEST_DATA / "assembly_bakta_output" / "assembly.json"
PROTEINS_JSON = TEST_DATA / "assembly_bakta_proteins_output" / "assembly.hypotheticals.json"

# The baktfold version is the only volatile content in these formats (the
# GenBank LOCUS date is the fixed 01-JAN-1980 placeholder). Normalise every
# context it appears in so the goldens stay version-independent and do NOT need
# refreshing on a version bump:
#   - "# Software: v<x>" / "Database: v<x>" headers (gff/gbff/embl/tsv/summary)
#   - "#Annotated with Baktfold (v<x>)" (inference.tsv / proteins.tsv)
#   - "ab initio prediction:Bakta:<major.minor>" (sORF inference, gbff/embl)
_SOFTWARE_DB_RE = re.compile(r"((?:Software|Database): v)\S+")
_VERSION = get_version().strip()
_VERSION_MM = ".".join(_VERSION.split(".")[:2]) if _VERSION else ""


def _normalize_version(text: str) -> str:
    text = _SOFTWARE_DB_RE.sub(r"\1<VERSION>", text)
    if _VERSION:
        text = text.replace(f"(v{_VERSION})", "(v<VERSION>)")
        text = text.replace(f"Bakta:{_VERSION_MM}", "Bakta:<VERSION>")
    return text


def _strip_after(text: str, marker: str) -> str:
    """Drop everything from the first occurrence of ``marker`` onward.

    Used for GFF3, which emits a single trailing ``##FASTA`` block after the
    features of *all* sequences.
    """
    idx = text.find(marker)
    return text if idx == -1 else text[:idx] + "\n"


def _strip_genbank_sequence(text: str) -> str:
    """Replace each record's ``ORIGIN..//`` nucleotide block with a stub.

    Operates per record so the feature tables of *every* contig are retained
    (a multi-record GenBank file has one ``ORIGIN`` per sequence).
    """
    return re.sub(r"(?ms)^ORIGIN.*?^//\s*$", "ORIGIN\n//", text)


def _strip_embl_sequence(text: str) -> str:
    """Replace each record's ``SQ..//`` nucleotide block with a stub (per record)."""
    return re.sub(r"(?ms)^SQ   .*?^//\s*$", "//", text)


def _headers_only(text: str) -> str:
    """Keep only FASTA ``>`` header lines (feature identity/order/products)."""
    return "\n".join(l for l in text.splitlines() if l.startswith(">")) + "\n"


# ext -> function turning the raw file body into the snapshotted text
_GENOME_PROCESSORS = {
    "gff3": lambda t: _normalize_version(_strip_after(t, "\n##FASTA")),
    "gbff": lambda t: _normalize_version(_strip_genbank_sequence(t)),
    "embl": lambda t: _normalize_version(_strip_embl_sequence(t)),
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
