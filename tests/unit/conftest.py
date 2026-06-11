"""
Shared fixtures for baktfold's function-level unit tests.

These tests are the fast, dependency-free safety net for the bug-fix pass
(see ../../REVIEW_FINDINGS.md) and the planned pandas->polars migration
(see ../../POLAR_MIGRATION.md). They pin the *current* behaviour of every
analytical / IO helper so any rewrite that drifts shape / dtype / value /
output-bytes lights up immediately.

They need NO database, NO foldseek, NO GPU, NO network — unlike the
end-to-end suite in tests/test_integration.py.

NOTE: importing any ``baktfold.*`` submodule currently triggers
``baktfold/__init__.py``, which eagerly imports torch (REVIEW_FINDINGS #29),
so test *collection* pays a ~6 s one-time import. Fixing #29 (lazy imports)
drops that to ~1 s. The tests themselves run in milliseconds.

Run::

    pytest tests/unit -v
    pytest tests/unit --snapshot-update   # bootstrap / accept new snapshots
"""
from __future__ import annotations

from pathlib import Path

import pytest
from loguru import logger

UNIT_DIR = Path(__file__).parent
FIXTURES = UNIT_DIR / "fixtures"
SNAPSHOTS = UNIT_DIR / "snapshots"


# ── logger isolation ────────────────────────────────────────────────────────


@pytest.fixture(autouse=True)
def _isolate_loguru():
    """Strip inherited loguru sinks before every unit test.

    tests/test_integration.py installs a global ``logger.add(sys.exit,
    level="ERROR")`` sink at import time. When the whole suite is collected
    (``pytest tests/``) that sink is live, so a *legitimate* ``logger.error``
    inside a function under test (e.g. ``get_genbank`` on a no-CDS GenBank)
    would raise ``SystemExit`` and fail the unit test. Removing all sinks
    around each unit test makes behaviour deterministic whether you run
    ``pytest tests/unit`` or ``pytest tests/``. Tests that need to assert on
    log output add their own sink inside the test (see test_create_foldseek_db).
    """
    logger.remove()
    yield
    logger.remove()


# ── path-shaped fixtures ────────────────────────────────────────────────────


@pytest.fixture(scope="session")
def fixtures_dir() -> Path:
    """Directory holding hand-crafted fixture inputs (TSV / faa / genbank)."""
    return FIXTURES


# ── snapshot helper: canonical TSV + per-module snapshot dir + assertion ────


def _df_to_tsv(df) -> str:
    """Render a DataFrame to canonical tab-separated text, no row index.

    Works for **both** pandas (current) and polars (post-migration) frames
    by duck-typing on ``write_csv``. This is deliberate: the snapshot encodes
    the exact output bytes today (pandas ``to_csv``); the polars migration
    (POLAR_MIGRATION.md, risk #1 — float formatting) must reproduce those
    same bytes or these tests fail, which is exactly the gate we want.
    """
    if hasattr(df, "write_csv"):  # polars.DataFrame
        return df.write_csv(separator="\t")
    return df.to_csv(sep="\t", index=False)  # pandas.DataFrame


@pytest.fixture
def df_snapshot(snapshot):
    """Snapshot-compare a DataFrame as canonical TSV.

    Snapshots live under ``tests/unit/snapshots/<subdir>/<name>``.

    Usage::

        def test_thing(df_snapshot):
            out = some_function(...)
            df_snapshot(out, "expected.tsv", subdir="test_thing_module")
    """

    def _snapshot(df, name: str, subdir: str) -> None:
        snapshot.snapshot_dir = SNAPSHOTS / subdir
        snapshot.assert_match(_df_to_tsv(df), name)

    return _snapshot


@pytest.fixture
def text_snapshot(snapshot):
    """Snapshot-compare arbitrary text (e.g. a generated TSV/FASTA file body)."""

    def _snapshot(text: str, name: str, subdir: str) -> None:
        snapshot.snapshot_dir = SNAPSHOTS / subdir
        snapshot.assert_match(text, name)

    return _snapshot
