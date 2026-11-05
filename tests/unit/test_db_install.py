"""Database install: versioned download names and release detection.

* Hugging Face downloads used a fixed filename, so publishing a new database
  release under it would make every older Baktfold fail its md5 check. The
  filename now comes from the release's VERSION_DICTIONARY entry.
* An install was judged complete by file names alone, so a 0.1.0 database
  (every file present, no pdb_chain table) would never be updated to 0.2.0.

Run::

    pytest tests/unit/test_db_install.py -v
"""
from __future__ import annotations

import sqlite3
from pathlib import Path
from unittest.mock import patch

from baktfold.databases import db


def _sqlite(path: Path, tables):
    conn = sqlite3.connect(path)
    for t in tables:
        conn.execute(f"CREATE TABLE {t} (id TEXT)")
    conn.commit(); conn.close()


def test_hf_download_uses_versioned_filename(tmp_path):
    blob = tmp_path / "blob"
    blob.write_bytes(b"x")
    seen = {}

    def fake_hf(**kw):
        seen.update(kw)
        return str(blob)

    with patch.object(db, "hf_hub_download", fake_hf):
        db.download(tmp_path / "out" / "baktfold_db_v0.2.0.tar.gz", tmp_path)
    assert seen["filename"] == "baktfold_db_v0.2.0.tar.gz"
    assert (tmp_path / "out" / "baktfold_db_v0.2.0.tar.gz").read_bytes() == b"x"


def test_required_tables_detect_old_release(tmp_path, monkeypatch):
    monkeypatch.setitem(db.VERSION_DICTIONARY, "9.9.9", {"required_tables": ["pdb_chain"]})
    monkeypatch.setattr(db, "CURRENT_DB_VERSION", "9.9.9")
    old = tmp_path / "old.db"; _sqlite(old, ["swissprot", "pdb"])
    new = tmp_path / "new.db"; _sqlite(new, ["swissprot", "pdb", "pdb_chain"])
    assert db.missing_required_tables(old) == ["pdb_chain"]
    assert db.missing_required_tables(new) == []
    assert db.missing_required_tables(tmp_path / "absent.db") == ["pdb_chain"]


def test_no_required_tables_for_current_release_is_noop(tmp_path):
    # the release in use must either list its tables or require none
    required = db.VERSION_DICTIONARY[db.CURRENT_DB_VERSION].get("required_tables", [])
    path = tmp_path / "baktfold.db"
    _sqlite(path, ["swissprot", "afdbclusters", "pdb", "cath"] + list(required))
    assert db.missing_required_tables(path) == []
