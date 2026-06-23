"""
Shared pytest configuration:
- command line options for the integration suite
- gating of the database-backed integration tests

baktfold's integration tests drive the real CLI against the full baktfold
database installed at ``tests/test_data/baktfold_db`` (~22 GB). That database
is far too large to provision in hosted CI, so the DB-backed tests are gated on
its presence: they run wherever the database exists (a dev machine or a
self-hosted runner) and skip gracefully where it does not. The DB-free
integration tests (``json`` / ``convert`` / ``download``) and the unit + golden
suites always run, so every push still gets a meaningful regression signal.

Set ``BAKTFOLD_SKIP_DB_TESTS=1`` to force the DB-backed tests to skip even when
the database is present (useful for a fast, DB-free local run).
"""
import inspect
import os
from pathlib import Path

import pytest

# Must match ``database_dir`` in tests/test_integration.py.
_DB_DIR = Path(__file__).parent / "test_data" / "baktfold_db"


def pytest_addoption(parser):
    parser.addoption("--gpu-available", action="store_true", dest="gpu_available")
    parser.addoption("--nvidia", action="store_true", dest="nvidia")
    parser.addoption("--threads", action="store", default=1, dest="threads")
    parser.addoption("--euks", action="store_true", dest="euks")


def _database_available() -> bool:
    """True when the baktfold test database is present (and not force-skipped)."""
    if os.environ.get("BAKTFOLD_SKIP_DB_TESTS"):
        return False
    return _DB_DIR.is_dir() and any(_DB_DIR.iterdir())


def pytest_collection_modifyitems(config, items):
    """Skip database-backed integration tests when the database is unavailable.

    A test is considered database-backed iff its source references
    ``database_dir`` (every ``run``/``compare``/``predict``/``proteins``/
    ``install``/``autotune`` integration test does; the ``json``/``convert``/
    ``download`` tests and the unit + golden tests do not).
    """
    if _database_available():
        return  # database present -> run everything, including DB-backed tests

    skip_db = pytest.mark.skip(
        reason=(
            f"baktfold database not found at {_DB_DIR} - skipping DB-backed "
            f"integration tests (DB-free, unit and golden tests still run). "
            f"Install it there, or unset BAKTFOLD_SKIP_DB_TESTS, to run them."
        )
    )
    for item in items:
        func = getattr(item, "function", None)
        if func is None:
            continue
        try:
            source = inspect.getsource(func)
        except (OSError, TypeError):
            continue
        if "database_dir" in source:
            item.add_marker(skip_db)
