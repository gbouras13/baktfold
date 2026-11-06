"""
Function-level tests for baktfold.features.create_foldseek_db.

The foldseek subprocess calls are monkeypatched out (``ExternalTool.run_tool``
-> no-op, and ``remove_file`` -> no-op so the generated TSVs survive for
inspection). What we actually test is the pure logic around them:

  - generate_foldseek_db_from_aa_3di: the aa/3di TSV generation + the
    mismatch-drop (REVIEW_FINDINGS #27) and ignore-extra-3Di behaviour.
  - generate_foldseek_db_from_structures: the per-CDS structure matching
    (REVIEW_FINDINGS #19) and the no-structure warning path (#30).

Run::

    pytest tests/unit/test_create_foldseek_db.py -v
"""
from __future__ import annotations

import pytest
from loguru import logger

import baktfold.features.create_foldseek_db as cfdb
from baktfold.utils.external_tools import ExternalTool


@pytest.fixture
def no_foldseek(monkeypatch):
    """Stub out the foldseek subprocess and the TSV cleanup, capturing the
    ExternalTool invocations so tests can assert on the commands built."""
    calls = []
    monkeypatch.setattr(ExternalTool, "run_tool", staticmethod(lambda tool, ctx=None: calls.append(tool)))
    # keep the temp TSVs on disk so we can read them back
    monkeypatch.setattr(cfdb, "remove_file", lambda p: None)
    return calls


@pytest.fixture
def capture_logs():
    """Capture loguru messages into a list for assertion."""
    msgs = []
    sink_id = logger.add(lambda m: msgs.append(m.record["message"]), level="WARNING")
    yield msgs
    logger.remove(sink_id)


# ─── generate_foldseek_db_from_aa_3di ───────────────────────────────────────


def _read_tsvs(db_dir):
    return {
        "aa": (db_dir / "aa.tsv").read_text(),
        "3di": (db_dir / "3di.tsv").read_text(),
        "header": (db_dir / "header.tsv").read_text(),
    }


def test_aa_3di_tsvs_snapshot(fixtures_dir, tmp_path, no_foldseek, text_snapshot):
    """cds_3 (no 3Di) is dropped; cds_X (3Di only) is ignored. The three TSVs
    end up with cds_1 + cds_2 only, 1-indexed and aligned."""
    cfdb.generate_foldseek_db_from_aa_3di(
        fixtures_dir / "db_aa.fasta",
        fixtures_dir / "db_3di.fasta",
        tmp_path,
        tmp_path,
        prefix="testdb",
    )
    tsvs = _read_tsvs(tmp_path)
    text_snapshot(tsvs["aa"], "aa.tsv", subdir="test_create_foldseek_db")
    text_snapshot(tsvs["3di"], "3di.tsv", subdir="test_create_foldseek_db")
    text_snapshot(tsvs["header"], "header.tsv", subdir="test_create_foldseek_db")


def test_aa_3di_mismatch_drops_unpaired_cds(fixtures_dir, tmp_path, no_foldseek):
    """cds_3 has no 3Di string -> dropped from all three TSVs (the #27 path)."""
    cfdb.generate_foldseek_db_from_aa_3di(
        fixtures_dir / "db_aa.fasta",
        fixtures_dir / "db_3di.fasta",
        tmp_path,
        tmp_path,
        prefix="testdb",
    )
    header = (tmp_path / "header.tsv").read_text()
    assert "cds_1" in header and "cds_2" in header
    assert "cds_3" not in header  # dropped (no 3Di)
    assert "cds_X" not in header  # ignored (no aa)


def test_aa_3di_calls_foldseek_three_times(fixtures_dir, tmp_path, no_foldseek):
    """One tsv2db per TSV: aa, 3di, header."""
    cfdb.generate_foldseek_db_from_aa_3di(
        fixtures_dir / "db_aa.fasta", fixtures_dir / "db_3di.fasta",
        tmp_path, tmp_path, prefix="testdb",
    )
    assert len(no_foldseek) == 3


# ─── generate_foldseek_db_from_structures ───────────────────────────────────


@pytest.fixture
def structure_dir(fixtures_dir, tmp_path):
    """A structures dir: cds_1.pdb + cds_2.cif present, cds_3 absent."""
    d = tmp_path / "structures"
    d.mkdir()
    (d / "cds_1.pdb").write_text("dummy pdb")
    (d / "cds_2.cif").write_text("dummy cif")
    return d


def test_structures_matches_pdb_and_cif(fixtures_dir, tmp_path, structure_dir, no_foldseek):
    """cds_1.pdb and cds_2.cif both match -> createdb is invoked (num_structures>0)."""
    cfdb.generate_foldseek_db_from_structures(
        fixtures_dir / "db_aa.fasta",
        tmp_path / "fsdb",
        structure_dir,
        tmp_path,
        prefix="testdb",
        proteins_flag=True,
    )
    # exactly one createdb call
    assert len(no_foldseek) == 1
    assert "createdb" in no_foldseek[0].command_as_str


def test_structures_missing_warns(fixtures_dir, tmp_path, structure_dir, no_foldseek, capture_logs):
    """cds_3 has no .pdb/.cif -> a 'No structure found' warning naming it."""
    (tmp_path / "fsdb").mkdir(exist_ok=True)
    cfdb.generate_foldseek_db_from_structures(
        fixtures_dir / "db_aa.fasta",
        tmp_path / "fsdb",
        structure_dir,
        tmp_path,
        prefix="testdb",
        proteins_flag=True,
    )
    assert any("cds_3" in m and "No structure found" in m for m in capture_logs)
    # cds_1 and cds_2 matched, so they are NOT warned about
    assert not any("cds_1" in m and "No structure found" in m for m in capture_logs)
