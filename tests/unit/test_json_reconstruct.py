"""
Function-level tests for the ``baktfold json`` reconstruction command
(baktfold.io.json_in.parse_baktfold_json_for_reconstruction + the CLI wiring
in baktfold.__init__).

``baktfold json`` takes a baktfold/bakta JSON and re-emits every non-Foldseek
output (.gff3/.gbff/.embl/.tsv/.inference.tsv/.faa/.ffn/.fna/.summary.txt/.json)
by calling the very same writers used by ``run``/``compare``. Because it is a
pure, deterministic transform of the JSON it is idempotent once the
self-describing ``baktfold_run`` provenance block is present.

These tests also pin the fix for the INSDC ``/transl_table`` bug: the call site
in io/io.py once passed ``other_genbank``/``translation_table`` in swapped order,
and the qualifier key was ``translation_table`` rather than ``transl_table`` —
together yielding ``/translation_table=False`` on every CDS and mis-routed RNA
inference strings. See io/io.py:write_bakta_outputs and io/insdc.py.

They need NO database, NO foldseek, NO GPU and NO network, and run against
committed (git-tracked) bakta JSON fixtures.

Run::

    pytest tests/unit/test_json_reconstruct.py -v
"""
from __future__ import annotations

import json as _json
from pathlib import Path

import pytest
from click.testing import CliRunner

import baktfold.bakta.constants as bc
from baktfold import main_cli
from baktfold.io.json_in import (
    _coerce_translation_table,
    _detect_duplicate_locus,
    _first_set,
    _infer_custom_db,
    _infer_euk,
)

REPO_ROOT = Path(__file__).resolve().parents[2]
TEST_DATA = REPO_ROOT / "tests" / "test_data"

# Tracked bakta genome fixture with CDS *and* tRNA/tmRNA features — exercises
# both the /transl_table qualifier and the RNA inference routing.
GENOME_JSON = TEST_DATA / "GCF_002368115.json"
# Tracked bakta-proteins fixture.
PROTEINS_JSON = TEST_DATA / "assembly_bakta_proteins_output" / "assembly.hypotheticals.json"

GENOME_OUTPUTS = [
    "gff3", "tsv", "gbff", "embl", "faa", "ffn", "fna", "inference.tsv", "summary.txt", "json",
]
PROTEINS_OUTPUTS = ["tsv", "faa", "summary.txt", "json"]


def _run_json(args):
    """Invoke ``baktfold json`` in-process; fail loudly with captured output."""
    result = CliRunner().invoke(main_cli, ["json", *args], catch_exceptions=False)
    assert result.exit_code == 0, f"baktfold json failed:\n{result.output}"
    return result


@pytest.fixture(scope="module")
def genome_recon(tmp_path_factory):
    """Reconstruct the genome fixture twice: pass 1 from the bakta JSON, pass 2
    from pass 1's regenerated JSON (which now carries the provenance block).
    Shared across the genome tests so the writers only run twice per module."""
    if not GENOME_JSON.exists():
        pytest.skip("genome fixture JSON not present")
    first = tmp_path_factory.mktemp("genome_a")
    second = tmp_path_factory.mktemp("genome_b")
    _run_json(["-i", str(GENOME_JSON), "-o", str(first), "-p", "g", "-f"])
    _run_json(["-i", str(first / "g.json"), "-o", str(second), "-p", "g", "-f"])
    return first, second


@pytest.fixture(scope="module")
def proteins_recon(tmp_path_factory):
    if not PROTEINS_JSON.exists():
        pytest.skip("proteins fixture JSON not present")
    first = tmp_path_factory.mktemp("proteins_a")
    second = tmp_path_factory.mktemp("proteins_b")
    _run_json(["-i", str(PROTEINS_JSON), "-o", str(first), "-p", "p", "-f"])
    _run_json(["-i", str(first / "p.json"), "-o", str(second), "-p", "p", "-f"])
    return first, second


# ── pure helper unit tests ──────────────────────────────────────────────────


def test_first_set_precedence():
    assert _first_set(None, None, True) is True
    assert _first_set(False, True) is False        # an explicit False wins over later values
    assert _first_set(None, None, None) is None


def test_infer_euk_from_feature_types():
    assert _infer_euk([{"type": bc.FEATURE_CDS}, {"type": bc.FEATURE_MRNA}]) is True
    assert _infer_euk([{"type": bc.FEATURE_GENE}]) is True
    assert _infer_euk([{"type": bc.FEATURE_CDS}, {"type": bc.FEATURE_T_RNA}]) is False


def test_infer_custom_db_from_dbxrefs():
    assert _infer_custom_db([{"db_xrefs": ["pdb:pdb_1abc", "custom:custom_42"]}]) is True
    assert _infer_custom_db([{"db_xrefs": ["afdb_v6:swissprot_x"]}, {}]) is False


def test_detect_duplicate_locus():
    feats = [
        {"type": bc.FEATURE_CDS, "hypothetical": True, "locus": "L1"},
        {"type": bc.FEATURE_CDS, "hypothetical": True, "locus": "L1"},
    ]
    assert _detect_duplicate_locus(feats) is True
    feats_unique = [
        {"type": bc.FEATURE_CDS, "hypothetical": True, "locus": "L1"},
        {"type": bc.FEATURE_CDS, "hypothetical": True, "locus": "L2"},
    ]
    assert _detect_duplicate_locus(feats_unique) is False


def test_coerce_translation_table():
    assert _coerce_translation_table(None) == 11
    assert _coerce_translation_table("4") == 4
    assert _coerce_translation_table(25) == 25


# ── genome-mode end-to-end reconstruction ───────────────────────────────────


def test_genome_reconstruction_produces_all_outputs(genome_recon):
    first, _ = genome_recon
    for ext in GENOME_OUTPUTS:
        assert (first / f"g.{ext}").exists(), f"missing g.{ext}"


def test_genome_reconstruction_persists_provenance(genome_recon):
    first, _ = genome_recon
    prov = _json.loads((first / "g.json").read_text()).get("baktfold_run")
    assert prov is not None and prov["mode"] == "genome"
    for key in ("euk", "custom_db", "fast", "translation_table", "prokka", "other_genbank"):
        assert key in prov


def test_genome_reconstruction_is_idempotent(genome_recon):
    """A 2nd pass over the regenerated JSON (carrying provenance, no overrides)
    reproduces all outputs byte-for-byte — including summary.txt and json."""
    first, second = genome_recon
    for ext in GENOME_OUTPUTS:
        assert (first / f"g.{ext}").read_bytes() == (second / f"g.{ext}").read_bytes(), \
            f"g.{ext} not idempotent"


def test_transl_table_qualifier_is_correct(genome_recon):
    """Regression: the INSDC genetic-code qualifier must be /transl_table=<int>
    (matching Bakta), NOT the non-standard /translation_table, and never the
    boolean that the swapped call-site args used to produce."""
    first, _ = genome_recon
    for ext in ("gbff", "embl"):
        text = (first / f"g.{ext}").read_text()
        assert "/transl_table=11" in text, f"/transl_table=11 missing from g.{ext}"
        assert "/translation_table" not in text, f"stale /translation_table qualifier in g.{ext}"
        assert "/transl_table=False" not in text and "/transl_table=True" not in text, \
            f"boolean translation table leaked into g.{ext}"


def test_bakta_rna_inference_routing(genome_recon):
    """Regression: for bakta input each tRNA /inference must be profile:tRNAscan:2.0
    (Bakta predicts tRNAs with tRNAscan-SE, matching the GFF source column) and
    each tmRNA must be profile:aragorn:1.2 (Aragorn). Guards against (a) the old
    swapped call-site forcing the other_genbank program strings, and (b) the
    tRNA/tmRNA inference values being swapped with the prokka branch."""
    fixture = _json.loads(GENOME_JSON.read_text())
    n_trna = sum(1 for f in fixture["features"] if f["type"] == bc.FEATURE_T_RNA)
    n_tmrna = sum(1 for f in fixture["features"] if f["type"] == bc.FEATURE_TM_RNA)
    assert n_trna > 0 and n_tmrna > 0, "fixture must contain tRNA and tmRNA to be meaningful"

    first, _ = genome_recon
    gbff = (first / "g.gbff").read_text()
    assert gbff.count("profile:tRNAscan:2.0") == n_trna     # every bakta tRNA
    assert gbff.count("profile:aragorn:1.2") == n_tmrna     # every bakta tmRNA
    assert "profile:tRNAscan-SE:2.0.12" not in gbff         # other_genbank program string must not leak
    assert "profile:INFERNAL:1.1.5" not in gbff


# ── proteins-mode end-to-end reconstruction ─────────────────────────────────


def test_proteins_reconstruction(proteins_recon):
    first, _ = proteins_recon
    for ext in PROTEINS_OUTPUTS:
        assert (first / f"p.{ext}").exists(), f"missing p.{ext}"
    for ext in ("gff3", "gbff", "embl", "ffn", "fna"):
        assert not (first / f"p.{ext}").exists(), f"unexpected p.{ext} in proteins mode"
    prov = _json.loads((first / "p.json").read_text()).get("baktfold_run")
    assert prov is not None and prov["mode"] == "proteins"


def test_proteins_reconstruction_is_idempotent(proteins_recon):
    first, second = proteins_recon
    for ext in PROTEINS_OUTPUTS:
        assert (first / f"p.{ext}").read_bytes() == (second / f"p.{ext}").read_bytes(), \
            f"p.{ext} not idempotent"
