"""db-scripts/init-pdb-chains.py builds the chain-level PDB name table.

Synthetic inputs in tests/unit/fixtures/pdb_chains cover: a SIFTS-mapped chain
named from Swiss-Prot; a symmetry copy ('_A-2') falling back to its base
chain's PDB entity name (all-caps name recased); a chain mapped to two UniProt
accessions (longest span wins, named from afdbclusters); uninformative names
at both levels (no row); and a chain in neither source (no row).

Run::

    pytest tests/unit/test_db_pdb_chains.py -v
"""
from __future__ import annotations

import shutil
import sqlite3
import subprocess
import sys
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[2] / "db-scripts" / "init-pdb-chains.py"
FIX = Path(__file__).parent / "fixtures" / "pdb_chains"


def test_build_pdb_chain_table(tmp_path):
    db = tmp_path / "baktfold.db"
    c = sqlite3.connect(db)
    for t in ("swissprot", "afdbclusters"):
        c.execute(f"CREATE TABLE {t} (id TEXT PRIMARY KEY, product TEXT)")
    c.execute("INSERT INTO swissprot VALUES ('P02185', 'Myoglobin')")
    c.execute("INSERT INTO afdbclusters VALUES ('Q00002', 'Immunoglobulin heavy constant gamma 1')")
    c.commit(); c.close()

    subprocess.run(
        [sys.executable, str(SCRIPT), "--seqres", str(FIX / "seqres.txt"),
         "--sifts", str(FIX / "sifts.tsv"), "--db", str(db),
         "--lookup", str(FIX / "pdb.lookup"), "--uniprot-names", str(FIX / "names.tsv")],
        check=True, capture_output=True,
    )
    rows = dict((r[0], tuple(r[1:])) for r in sqlite3.connect(db).execute(
        "SELECT id, product, source, uniprot FROM pdb_chain"))
    assert rows == {
        "101m_A": ("Myoglobin", "uniprot", "P02185"),
        "1abc_A": ("DNA polymerase III subunit beta", "pdb_entity", None),
        "1abc_B": ("Immunoglobulin heavy constant gamma 1", "uniprot", "Q00002"),
    }
