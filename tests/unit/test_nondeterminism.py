"""
Regression tests for the two run-to-run instabilities that made the golden
output tests flake independently of any code change.

1. Foldseek emits hits per query in descending-bitscore order, but the order
   *within* a run of equal bitscores is not stable. ``get_tophit`` took the
   first row per query, so a tie between two equally-scoring targets resolved
   differently from run to run — the goldens flipped between the PDB entries
   "MsDpo4-DNA complex 1" and "MsDpo4-DNA complex 2".

2. ProstT5 runs in fp16 on the GPU, so per-residue probabilities move by ~0.5
   (0-100 scale) between identical runs. A residue whose confidence sits on
   --mask-threshold is masked to 'X' in one run and left as itself in the next,
   so the masked AA/3Di FASTAs differ by a single character.

Run::

    pytest tests/unit/test_nondeterminism.py -v
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

import polars as pl

from baktfold.results.tophit import get_tophit

_spec = importlib.util.spec_from_file_location(
    "_cmp", Path(__file__).parent.parent / "compare_outputs.py"
)
_cmp = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_cmp)


# ── 1. deterministic tie-breaking ───────────────────────────────────────────

def _write_tsv(tmp_path: Path, rows: list) -> Path:
    """Foldseek convertalis output: query target bits fident evalue q/t coords."""
    p = tmp_path / "foldseek_results.tsv"
    p.write_text("".join("\t".join(str(c) for c in r) + "\n" for r in rows))
    return p


def _row(query: str, target: str, bits: int):
    return [query, target, bits, 0.5, 1e-10, 1, 100, 100, 1, 100, 100]


def test_tied_bitscores_resolve_to_the_same_target_either_way(tmp_path):
    """A tie must pick the same target regardless of the order Foldseek emitted."""
    forward = get_tophit(_write_tsv(tmp_path, [
        _row("q1", "MsDpo4-DNA-complex-1", 500),
        _row("q1", "MsDpo4-DNA-complex-2", 500),
    ]), structures=False)
    reversed_ = get_tophit(_write_tsv(tmp_path, [
        _row("q1", "MsDpo4-DNA-complex-2", 500),
        _row("q1", "MsDpo4-DNA-complex-1", 500),
    ]), structures=False)

    assert forward["target"].to_list() == reversed_["target"].to_list()
    assert forward["target"][0] == "MsDpo4-DNA-complex-1"   # tie broken by target


def test_higher_bitscore_still_wins_over_the_tie_break(tmp_path):
    """The tie-break must not outrank the score itself."""
    df = get_tophit(_write_tsv(tmp_path, [
        _row("q1", "zzz_high_score", 900),
        _row("q1", "aaa_low_score", 100),
    ]), structures=False)
    assert df["target"][0] == "zzz_high_score"


def test_query_order_is_preserved(tmp_path):
    """Queries keep Foldseek's original order — only ties inside a query move."""
    df = get_tophit(_write_tsv(tmp_path, [
        _row("zq", "t1", 500),
        _row("aq", "t2", 500),
    ]), structures=False)
    assert df["query"].to_list() == ["zq", "aq"]   # not alphabetised


def test_cath_keeps_all_hits_but_orders_ties(tmp_path):
    """CATH retains every greedy hit; ties within a query are still ordered."""
    df = get_tophit(_write_tsv(tmp_path, [
        _row("q1", "b_target", 500),
        _row("q1", "a_target", 500),
    ]), structures=False, cath=True)
    assert df.height == 2
    assert df["target"].to_list() == ["a_target", "b_target"]


# ── 2. borderline mask flips ────────────────────────────────────────────────

def test_single_mask_flip_is_tolerated():
    """One residue crossing the confidence cutoff is noise, not a difference."""
    dev = [">MEGJMN_078", "MNTNLXLTADXVHISMPAGAYLXVXXRXYXHIP"]
    ref = [">MEGJMN_078", "MNTNLXLTADXVHISMPAGAYLXVXXXXYXHIP"]
    assert _cmp._fasta_differ(dev, ref) == []


def test_isolated_substitution_is_within_budget():
    """ProstT5 flips the predicted state at a handful of positions run to run
    (measured: ~0.005% of residues, and fp32 does not fix it), so one or two
    substitutions per record are budgeted rather than failed."""
    dev = [">MEGJMN_078", "MNTNLXLTADXVHISMPAGAYLXVXXRXYXHIP"]
    ref = [">MEGJMN_078", "MNTNLXLTADXVHISMPAGAYLXVXXKXYXHIP"]
    assert _cmp._fasta_differ(dev, ref) == []


def test_substitutions_beyond_the_budget_are_reported():
    """Past the per-record budget it is no longer noise."""
    dev = [">s", "ACDEFGHIKL" * 10]
    ref = [">s", "WWWWWWWWWW" + ("ACDEFGHIKL" * 10)[10:]]
    diffs = _cmp._fasta_differ(dev, ref)
    assert diffs and "substitution(s)" in diffs[0]


def test_file_level_guard_catches_systematic_drift():
    """Many records each within their own budget still fail in aggregate —
    that pattern is a model/dependency change, not run-to-run noise."""
    dev = [line for i in range(200) for line in (f">s{i}", "ACDEFGHIKL" * 5)]
    ref = [line for i in range(200) for line in (f">s{i}", "W" + ("ACDEFGHIKL" * 5)[1:])]
    diffs = _cmp._fasta_differ(dev, ref)
    assert any("file tolerance" in d for d in diffs)


# ── 3. near-tied hit descriptions in the annotation JSON ────────────────────

def _ann(description: str, product: str = "DNA polymerase IV") -> dict:
    return {"features": [{"locus": "L1", "product": product,
                          "pstc": [{"source": "pdb", "description": description}]}]}


def _write_json(tmp_path: Path, name: str, obj: dict) -> Path:
    import json as _json
    p = tmp_path / name
    p.write_text(_json.dumps(obj))
    return p


def test_near_tied_hit_description_is_budgeted(tmp_path):
    """Foldseek scores wobble with the 3Di strings, so which of two near-equal
    hits wins can flip — the accession is already normalised away, and the
    description betrays the same flip."""
    fd = _write_json(tmp_path, "dev.json", _ann("MsDpo4-DNA complex 1"))
    fr = _write_json(tmp_path, "ref.json", _ann("MsDpo4-DNA complex 2"))
    assert _cmp._compare_annotation_json(fd, fr, normalize=True) == []


def test_changed_product_is_still_reported(tmp_path):
    """The annotation itself changing is never noise."""
    fd = _write_json(tmp_path, "dev.json", _ann("MsDpo4-DNA complex 1", product="DNA polymerase IV"))
    fr = _write_json(tmp_path, "ref.json", _ann("MsDpo4-DNA complex 1", product="hypothetical protein"))
    diffs = _cmp._compare_annotation_json(fd, fr, normalize=True)
    assert diffs and "product" in diffs[0]


def test_many_description_changes_are_reported(tmp_path):
    """Past the budget, wholesale hit changes are real drift."""
    dev = {"features": [{"locus": f"L{i}", "pstc": [{"source": "pdb", "description": f"hit {i}"}]}
                        for i in range(20)]}
    ref = {"features": [{"locus": f"L{i}", "pstc": [{"source": "pdb", "description": f"other {i}"}]}
                        for i in range(20)]}
    fd, fr = _write_json(tmp_path, "d.json", dev), _write_json(tmp_path, "r.json", ref)
    diffs = _cmp._compare_annotation_json(fd, fr, normalize=True)
    assert diffs and "pstc description changes" in diffs[0]


def test_systematic_mask_change_is_reported():
    """Wholesale masking changes must still fail — that is real drift."""
    dev = [">s", "ACDEFGHIKLMNPQRSTVWY" * 5]
    ref = [">s", "XXXXXXXXXX" + ("ACDEFGHIKLMNPQRSTVWY" * 5)[10:]]
    diffs = _cmp._fasta_differ(dev, ref)
    assert diffs and "mask flips" in diffs[0]


def test_length_and_header_changes_are_reported():
    assert _cmp._fasta_differ([">a", "ACDE"], [">b", "ACDE"])      # header
    assert _cmp._fasta_differ([">a", "ACDE"], [">a", "ACD"])       # length
