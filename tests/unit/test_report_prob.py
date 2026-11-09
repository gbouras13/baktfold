"""Tests for the --report-prob plumbing (Foldseek homology probability).

Foldseek's `prob` column is the homology probability referred to in the
Foldseek paper's fold-level calibration. It is opt-in because emitting it
requires the backtrace CIGAR (`-a 1`) on the search — without that, convertalis
aborts with "Backtrace cigar is missing in the alignment result" — which costs
extra alignment time and result-DB size.

These tests pin three things:
  1. the Foldseek command gains `-a 1` when prob is requested (and only then),
  2. the format-output column order matches what get_tophit parses,
  3. the parsed value reaches the feature dict via pstc.parse.

Run::

    pytest tests/unit/test_report_prob.py -v
"""
from __future__ import annotations

import shutil
import subprocess
from pathlib import Path
from unittest.mock import patch

import polars as pl
import pytest

from baktfold.features import run_foldseek
from baktfold.results.tophit import get_tophit


# ─── Foldseek command construction ──────────────────────────────────────────


def _captured_params(**kwargs):
    """Run the given foldseek helper with ExternalTool stubbed; return params."""
    captured = {}

    class FakeTool:
        def __init__(self, **kw):
            captured.update(kw)

        @staticmethod
        def run_tool(tool):
            return None

    with patch.object(run_foldseek, "ExternalTool", FakeTool):
        kwargs["fn"](**kwargs["call"])
    return captured["params"]


def _search_params(*, structures: bool, report_prob: bool):
    return _captured_params(
        fn=run_foldseek.run_foldseek_search,
        call=dict(
            query_db="qdb", target_db="tdb", result_db="res", temp_db="tmp",
            threads=4, logdir="log", evalue=1e-3, sensitivity=9.5, max_seqs=1000,
            ultra_sensitive=False, extra_foldseek_params="", foldseek_gpu=False,
            structures=structures, report_prob=report_prob,
        ),
    )


def _convertalis_params(*, structures: bool, report_prob: bool):
    return _captured_params(
        fn=run_foldseek.create_result_tsv,
        call=dict(
            query_db="qdb", target_db="tdb", result_db="res", result_tsv="out.tsv",
            logdir="log", foldseek_gpu=False, structures=structures, threads=4,
            report_prob=report_prob,
        ),
    )


class TestSearchCommand:
    def test_default_search_has_no_alignment_flag(self):
        """Default ProstT5 path stays cheap — no -a 1."""
        assert "-a 1" not in _search_params(structures=False, report_prob=False)

    def test_report_prob_forces_alignment_flag(self):
        """prob needs the backtrace CIGAR, so -a 1 must be added."""
        assert "-a 1" in _search_params(structures=False, report_prob=True)

    def test_structures_still_get_alignment_flag(self):
        """Pre-existing behaviour: TM-score/LDDT already required -a 1."""
        assert "-a 1" in _search_params(structures=True, report_prob=False)

    def test_structures_and_prob_add_flag_once(self):
        params = _search_params(structures=True, report_prob=True)
        assert params.count("-a 1") == 1


class TestConvertalisColumns:
    def test_default_columns_unchanged(self):
        params = _convertalis_params(structures=False, report_prob=False)
        assert "query,target,bits,fident,evalue,qstart,qend,qlen,tstart,tend,tlen" in params
        assert "prob" not in params

    def test_prob_appended_last(self):
        params = _convertalis_params(structures=False, report_prob=True)
        assert params.endswith("--threads 4")
        assert "tlen,prob" in params

    def test_structures_prob_order_is_tmscore_lddt_prob(self):
        """get_tophit appends alntmscore, lddt then prob — order must match."""
        params = _convertalis_params(structures=True, report_prob=True)
        assert "tlen,alntmscore,lddt,prob" in params


# ─── parsing ────────────────────────────────────────────────────────────────


class TestGetTophitProb:
    def test_prob_column_parsed(self, fixtures_dir):
        df = get_tophit(
            fixtures_dir / "foldseek_prob.tsv", structures=False, cath=False,
            report_prob=True,
        )
        assert "prob" in df.columns
        assert len(df) == 3
        q1 = df.filter(pl.col("query") == "query_0001").row(0, named=True)
        assert float(q1["prob"]) == 1.0

    def test_prob_is_last_column(self, fixtures_dir):
        """The reorder in get_tophit must leave prob in the tail."""
        df = get_tophit(
            fixtures_dir / "foldseek_prob.tsv", structures=False, cath=False,
            report_prob=True,
        )
        assert df.columns[-1] == "prob"

    def test_structures_prob_column_order(self, fixtures_dir):
        df = get_tophit(
            fixtures_dir / "foldseek_structures_prob.tsv", structures=True,
            cath=False, report_prob=True,
        )
        assert df.columns[-3:] == ["alntmscore", "lddt", "prob"]

    def test_without_flag_no_prob_column(self, fixtures_dir):
        """Default parse of a default TSV is untouched by this feature."""
        df = get_tophit(
            fixtures_dir / "foldseek_basic.tsv", structures=False, cath=False
        )
        assert "prob" not in df.columns


# ─── pstc hand-off ──────────────────────────────────────────────────────────


class TestPstcCarriesProb:
    def _feature(self):
        return [{
            "type": "cds", "locus": "q1", "id": "q1",
            "aa_hexdigest": None, "start": 1, "stop": 300, "strand": "+",
        }]

    def test_prob_attached_to_pstc_hit(self):
        from baktfold.bakta import pstc

        df = pl.DataFrame({
            "query": ["q1"], "target": ["AF-P12345-F1-model_v4"],
            "bitscore": [520.0], "fident": [0.95], "evalue": ["1.2e-50"],
            "qStart": [1], "qEnd": [100], "qLen": [100], "qCov": [0.99],
            "tStart": [1], "tEnd": [100], "tLen": [100], "tCov": [0.99],
            "prob": [0.981],
        })
        feats = self._feature()
        pstc.parse(feats, df, "swissprot")
        hits = feats[0].get("pstc")
        assert hits, "expected a pstc hit"
        hit = hits[0] if isinstance(hits, list) else hits
        assert hit["prob"] == pytest.approx(0.981)

    def test_absent_prob_leaves_key_off(self):
        from baktfold.bakta import pstc

        df = pl.DataFrame({
            "query": ["q1"], "target": ["AF-P12345-F1-model_v4"],
            "bitscore": [520.0], "fident": [0.95], "evalue": ["1.2e-50"],
            "qStart": [1], "qEnd": [100], "qLen": [100], "qCov": [0.99],
            "tStart": [1], "tEnd": [100], "tLen": [100], "tCov": [0.99],
        })
        feats = self._feature()
        pstc.parse(feats, df, "swissprot")
        hits = feats[0].get("pstc")
        assert hits, "expected a pstc hit"
        hit = hits[0] if isinstance(hits, list) else hits
        assert "prob" not in hit


# ─── ProstT5 path rejected up front ─────────────────────────────────────────


class TestReportProbNeedsStructures:
    """prob needs query C-alpha coordinates, which ProstT5 query DBs lack (0.4.0 bug)."""

    def test_validator(self):
        import click

        from baktfold.utils.validation import check_report_prob

        check_report_prob(False, structures=False)
        check_report_prob(True, structures=True)
        with pytest.raises(click.UsageError, match="--structure-dir"):
            check_report_prob(True, structures=False)

    @pytest.mark.parametrize("argv", [
        ["run", "-i", "in.json"],
        ["proteins", "-i", "in.faa"],
        ["compare", "-i", "in.json", "--predictions-dir", "preds"],
        ["proteins-compare", "-i", "in.faa", "--predictions-dir", "preds"],
    ])
    def test_cli_rejects_before_any_work(self, tmp_path, argv):
        from click.testing import CliRunner

        from baktfold import main_cli

        out = tmp_path / "out"
        r = CliRunner().invoke(main_cli, argv + ["-o", str(out), "--report-prob"])
        assert r.exit_code == 2, r.output
        assert "--structure-dir" in r.output
        assert not out.exists()


# ─── real Foldseek ──────────────────────────────────────────────────────────

PDBS = Path(__file__).parents[1] / "test_data" / "pdbs"


@pytest.fixture
def structure_db(tmp_path):
    if shutil.which("foldseek") is None:
        pytest.skip("foldseek not on PATH")
    pdbs = sorted(PDBS.glob("*.pdb"))[:4]
    d = tmp_path / "pdbs"
    d.mkdir()
    for p in pdbs:
        shutil.copy(p, d / p.name)
    db = tmp_path / "qdb"
    subprocess.run(["foldseek", "createdb", str(d), str(db)], check=True, capture_output=True)
    return db, [p.stem for p in pdbs]


def _search_and_convert(tmp_path, db, *, structures, report_prob):
    logdir = tmp_path / "logs"
    logdir.mkdir(exist_ok=True)
    res, tsv = tmp_path / "res", tmp_path / "res.tsv"
    run_foldseek.run_foldseek_search(
        query_db=db, target_db=db, result_db=res, temp_db=tmp_path / "tmp", threads=1,
        logdir=logdir, evalue=1e-3, sensitivity=9.5, max_seqs=1000, ultra_sensitive=False,
        extra_foldseek_params="", foldseek_gpu=False, structures=structures, report_prob=report_prob)
    run_foldseek.create_result_tsv(db, db, res, tsv, logdir, False, structures, 1, report_prob=report_prob)
    return tsv


class TestRealFoldseek:
    def test_structures_report_prob(self, tmp_path, structure_db):
        """The supported path: structure queries, -a 1, prob column parsed by get_tophit."""
        db, ids = structure_db
        tsv = _search_and_convert(tmp_path, db, structures=True, report_prob=True)
        df = get_tophit(tsv, structures=True, cath=False, report_prob=True)
        assert df.columns[-3:] == ["alntmscore", "lddt", "prob"]
        assert set(df["query"]) <= set(ids) and len(df) > 0
        probs = df["prob"].cast(pl.Float64)
        assert probs.is_between(0, 1).all()
        # each query's best hit is itself
        assert (probs > 0.9).all()

    def test_query_without_ca_cannot_report_prob(self, tmp_path, structure_db):
        """Why the validator exists: a query DB without <db>_ca (as built from ProstT5
        3Di) fails in Foldseek when prob is requested."""
        db, _ = structure_db
        for f in tmp_path.glob("qdb_ca*"):
            f.unlink()
        # control: the same DB searches fine without prob, so the failure below is prob's
        assert _search_and_convert(tmp_path, db, structures=False, report_prob=False).stat().st_size > 0
        for f in tmp_path.glob("res*"):
            f.unlink() if f.is_file() else shutil.rmtree(f)
        with pytest.raises(SystemExit):
            _search_and_convert(tmp_path, db, structures=False, report_prob=True)
