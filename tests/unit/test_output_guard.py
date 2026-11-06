"""An existing output must stop the run unless -f/--force is given.

instantiate_dirs / validate_outfile / check_genbank_and_prokka run *before*
begin_baktfold installs the exit-on-ERROR log sink, so after the sink-based
refactor (REVIEW_FINDINGS B7) their logger.error calls no longer stopped
anything: `baktfold json -o existing_dir` (no -f) logged an ERROR, wrote into
the old directory alongside stale files, and exited 0.

Run::

    pytest tests/unit/test_output_guard.py -v
"""
from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from baktfold.utils.validation import (
    check_genbank_and_prokka, instantiate_dirs, validate_outfile,
)

DATA = Path(__file__).resolve().parents[1] / "test_data"


def test_existing_dir_without_force_exits(tmp_path):
    out = tmp_path / "out"
    out.mkdir()
    (out / "stale.tsv").write_text("old\n")
    with pytest.raises(SystemExit) as exc:
        instantiate_dirs(out, force=False)
    assert exc.value.code == 1
    assert (out / "stale.tsv").exists()  # nothing touched


def test_existing_dir_with_force_is_recreated(tmp_path):
    out = tmp_path / "out"
    out.mkdir()
    (out / "stale.tsv").write_text("old\n")
    instantiate_dirs(out, force=True)
    assert out.exists() and not (out / "stale.tsv").exists()


def test_new_dir_is_created(tmp_path):
    instantiate_dirs(tmp_path / "new", force=False)
    assert (tmp_path / "new").is_dir()


def test_existing_outfile_without_force_exits(tmp_path):
    f = tmp_path / "out.json"
    f.write_text("{}")
    with pytest.raises(SystemExit) as exc:
        validate_outfile(f, force=False)
    assert exc.value.code == 1
    assert f.read_text() == "{}"


def test_non_genbank_input_exits(tmp_path):
    f = tmp_path / "not.gbk"
    f.write_text(">seq\nMKV\n")
    with pytest.raises(SystemExit) as exc:
        check_genbank_and_prokka(f, euk=False)
    assert exc.value.code == 1


def test_cli_json_refuses_existing_output(tmp_path):
    """End to end: the exit status is what a pipeline/Snakemake rule sees."""
    src = DATA / "golden" / "run_all" / "baktfold.json"
    out = tmp_path / "out"
    out.mkdir()
    (out / "marker").write_text("x")
    exe = Path(sys.executable).with_name("baktfold")
    exe = str(exe) if exe.exists() else shutil.which("baktfold")
    if exe is None:
        pytest.skip("baktfold console script not installed")
    proc = subprocess.run(
        [exe, "json", "-i", str(src), "-o", str(out)],
        capture_output=True, text=True,
    )
    assert proc.returncode != 0
    assert (out / "marker").exists()
    assert not (out / "baktfold.json").exists()
