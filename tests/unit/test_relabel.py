"""``baktfold json --relabel``: rerun label selection on the hits stored in a JSON.

Outputs written by earlier versions carry their Foldseek hits (``pstc``) but
were labelled with older rules: "Uncharacterized protein ..." from a higher
database masked functional lower-database hits (audit A1), CATH superfamily
names were used verbatim, and PDB hits were labelled with the entry title.
--relabel undoes the earlier labelling and reruns combine_annotation on the
stored hits, so no search is needed.

Run::

    pytest tests/unit/test_relabel.py -v
"""
from __future__ import annotations

import copy
import json as _json
import sqlite3
from pathlib import Path

import pytest
from click.testing import CliRunner

import baktfold.bakta.annotation as anno
from baktfold import main_cli
from baktfold.subcommands.relabel import relabel_features

REPO_ROOT = Path(__file__).resolve().parents[2]
GENOME_JSON = REPO_ROOT / "tests" / "test_data" / "GCF_002368115.json"

BAKTA_XREFS = ["SO:0001217", "UniRef:UniRef50_A0A845YRR5"]
PDB_TITLE = "SPERM WHALE MYOGLOBIN F46V N-BUTYL ISOCYANIDE AT PH 9.0"


def _old_v010_feature(pstc, product, xrefs_from_hits, hypothetical=False, locus="L1"):
    """A CDS as v0.1.0 wrote it: Baktfold-marked, with xrefs from its hits."""
    f = {
        "type": "cds", "id": f"ID_{locus}", "locus": locus, "product": product,
        "db_xrefs": sorted(BAKTA_XREFS + xrefs_from_hits), "pstc": pstc, "baktfold": True,
    }
    if hypothetical:
        f["hypothetical"] = True
    return f


def _uncharacterised_over_afdb():
    # v0.1.0: Swiss-Prot "Uncharacterized protein" won, was renamed to
    # "hypothetical protein", and the CDS was left unmarked as hypothetical.
    return _old_v010_feature(
        [
            {"source": "swissprot", "id": "P1", "description": "Uncharacterized protein YxeA"},
            {"source": "afdb", "id": "A1", "description": "IS66 family transposase"},
        ],
        "hypothetical protein",
        ["afdb_v6:swissprot_P1", "afdb_v6:afdbclusters_A1"],
    )


def _db(tmp_path, chain_table=True):
    d = tmp_path / "db"
    d.mkdir()
    conn = sqlite3.connect(d / "baktfold.db")
    for t in ("swissprot", "afdbclusters", "pdb"):
        conn.execute(f"CREATE TABLE {t} (id TEXT PRIMARY KEY, product TEXT)")
    conn.execute("CREATE TABLE cath (id TEXT PRIMARY KEY, cath_class TEXT, product TEXT)")
    conn.execute("INSERT INTO swissprot VALUES ('P1', 'Uncharacterized protein YxeA')")
    conn.execute("INSERT INTO afdbclusters VALUES ('A1', 'IS66 family transposase')")
    conn.execute("INSERT INTO pdb VALUES ('101m', 'SPERM WHALE MYOGLOBIN F46V N-BUTYL ISOCYANIDE AT PH 9.0')")
    if chain_table:
        conn.execute("CREATE TABLE pdb_chain (id TEXT PRIMARY KEY, product TEXT, source TEXT)")
        conn.execute("INSERT INTO pdb_chain VALUES ('101m_A', 'Myoglobin', 'uniprot')")
    conn.commit()
    conn.close()
    return d


def _pdb_feature():
    return _old_v010_feature(
        [{"source": "pdb", "id": "101m", "description": PDB_TITLE}],
        anno.revise_cds_product(PDB_TITLE),
        ["pdb:pdb_101m"],
    )


# ── reset ────────────────────────────────────────────────────────────────────


def test_reset_removes_only_baktfold_xrefs():
    f = _uncharacterised_over_afdb()
    assert anno.reset_baktfold_annotation(f)
    assert f["db_xrefs"] == BAKTA_XREFS
    assert f["product"] == "hypothetical protein"
    assert f["hypothetical"] is True
    assert "baktfold" not in f and "annotation_confidence" not in f
    assert len(f["pstc"]) == 2  # hits are kept


def test_reset_keeps_input_gene_fields():
    f = _uncharacterised_over_afdb()
    f.update(gene="yxeA", genes=["yxeA"])
    anno.reset_baktfold_annotation(f)
    assert f["gene"] == "yxeA" and f["genes"] == ["yxeA"]


def test_reset_leaves_unlabelled_features_alone():
    f = {"type": "cds", "product": "DNA gyrase subunit A", "db_xrefs": list(BAKTA_XREFS)}
    before = copy.deepcopy(f)
    assert not anno.reset_baktfold_annotation(f)
    assert f == before


def test_reset_then_combine_matches_combine_on_fresh_hits():
    """Relabelling an old feature gives what run would give on the same hits."""
    old = _uncharacterised_over_afdb()
    fresh = {
        "type": "cds", "id": old["id"], "locus": old["locus"], "product": "hypothetical protein",
        "hypothetical": True, "db_xrefs": list(BAKTA_XREFS), "pstc": copy.deepcopy(old["pstc"]),
    }
    anno.reset_baktfold_annotation(old)
    for f in (old, fresh):
        anno.combine_annotation(f, fast=False)
    assert old == fresh


# ── relabel_features ─────────────────────────────────────────────────────────


def test_uninformative_hit_no_longer_masks_a_functional_one():
    f = _uncharacterised_over_afdb()
    info = relabel_features([f], fast=False, has_duplicate_locus=False)
    assert f["product"] == "IS66 family transposase"
    assert "hypothetical" not in f
    assert f["db_xrefs"] == sorted(BAKTA_XREFS + ["afdb_v6:swissprot_P1", "afdb_v6:afdbclusters_A1"])
    assert info["labels_changed"] == 1 and info["functional"] == 1


def test_only_uninformative_hits_stay_hypothetical():
    f = _old_v010_feature(
        [{"source": "afdb", "id": "A2", "description": "Uncharacterized protein"}],
        "hypothetical protein", ["afdb_v6:afdbclusters_A2"],
    )
    relabel_features([f], fast=False, has_duplicate_locus=False)
    assert f["product"] == "hypothetical protein"
    assert f["hypothetical"] is True
    assert "annotation_confidence" not in f


def test_cath_only_hit_becomes_domain_label():
    f = _old_v010_feature(
        [{"source": "cath", "id": "1abcA01", "description": "SH2 domain"}], "SH2 domain", ["cath:cath_1abcA01"],
    )
    relabel_features([f], fast=False, has_duplicate_locus=False)
    assert f["product"] == "SH2 domain-containing protein"


def test_relabel_is_idempotent():
    f = _uncharacterised_over_afdb()
    relabel_features([f], fast=False, has_duplicate_locus=False)
    once = copy.deepcopy(f)
    relabel_features([f], fast=False, has_duplicate_locus=False)
    assert f == once


def test_pdb_without_chain_keeps_entry_title(tmp_path):
    f = _pdb_feature()
    info = relabel_features([f], fast=False, has_duplicate_locus=False, database=_db(tmp_path))
    assert f["product"] == anno.revise_cds_product(PDB_TITLE)  # Bakta cleanup: "9.0" -> "9-0"
    assert info["pdb_no_chain_kept"] == 1


def test_pdb_chain_from_tophit_gives_chain_name(tmp_path):
    tophit = tmp_path / "baktfold_pdb_tophit.tsv"
    tophit.write_text("query\ttarget\tbitscore\nL1\t101m-assembly1.cif.gz_A\t100\n")
    f = _pdb_feature()
    info = relabel_features([f], fast=False, has_duplicate_locus=False,
                            database=_db(tmp_path), pdb_tophit=tophit)
    assert f["product"] == "Myoglobin"
    assert f["pstc"][0]["chain"] == "A"
    assert "pdb:pdb_101m_A" in f["db_xrefs"] and "pdb:pdb_101m" not in f["db_xrefs"]
    assert info["pdb_chains_added"] == 1


def test_unlisted_pdb_chain_falls_through(tmp_path):
    f = _pdb_feature()
    f["pstc"][0]["chain"] = "Z"  # not in pdb_chain
    relabel_features([f], fast=False, has_duplicate_locus=False, database=_db(tmp_path))
    assert f["product"] == "hypothetical protein" and f["hypothetical"] is True


def test_database_without_chain_table_keeps_titles(tmp_path):
    f = _pdb_feature()
    relabel_features([f], fast=False, has_duplicate_locus=False, database=_db(tmp_path, chain_table=False))
    assert f["product"] == anno.revise_cds_product(PDB_TITLE)  # Bakta cleanup: "9.0" -> "9-0"


def test_custom_db_descriptions_survive_database_refresh(tmp_path):
    f = _old_v010_feature(
        [{"source": "custom_db", "id": "arc_0001", "description": "S-layer protein"}],
        "S-layer protein", ["custom:custom_arc_0001"],
    )
    relabel_features([f], fast=False, has_duplicate_locus=False, database=_db(tmp_path))
    assert f["product"] == "S-layer protein"


def test_confidence_needs_prostt5_confidence(tmp_path):
    hit = {"source": "afdb", "id": "A1", "description": "IS66 family transposase",
           "query_cov": 0.95, "subject_cov": 0.95, "identity": 0.5, "evalue": 1e-30}
    f = _old_v010_feature([dict(hit)], "IS66 family transposase", ["afdb_v6:afdbclusters_A1"])
    info = relabel_features([f], fast=False, has_duplicate_locus=False)
    assert "annotation_confidence" not in f and info["confidence_dropped"] == 1

    probs = tmp_path / "baktfold_prostT5_3di_mean_probabilities.csv"
    probs.write_text("L1,72.0\n")
    f = _old_v010_feature([dict(hit)], "IS66 family transposase", ["afdb_v6:afdbclusters_A1"])
    relabel_features([f], fast=False, has_duplicate_locus=False, prostt5_probs=probs)
    assert f["annotation_confidence"] == "high"


def test_structure_hits_keep_confidence():
    hit = {"source": "afdb", "id": "A1", "description": "IS66 family transposase", "tmscore": 0.8,
           "query_cov": 0.95, "subject_cov": 0.95, "identity": 0.5, "evalue": 1e-30}
    f = _old_v010_feature([hit], "IS66 family transposase", ["afdb_v6:afdbclusters_A1"])
    relabel_features([f], fast=False, has_duplicate_locus=False)
    assert f["annotation_confidence"] == "high" and f["tmscore"] == 0.8


# ── CLI ──────────────────────────────────────────────────────────────────────


@pytest.fixture
def old_genome_json(tmp_path):
    """The tracked Bakta fixture with one hypothetical CDS given v0.1.0-style hits."""
    if not GENOME_JSON.exists():
        pytest.skip("genome fixture JSON not present")
    data = _json.loads(GENOME_JSON.read_text())
    cds = next(f for f in data["features"] if f["type"] == "cds" and f.get("hypothetical"))
    old = _uncharacterised_over_afdb()
    cds.update(pstc=old["pstc"], baktfold=True, product="hypothetical protein",
               db_xrefs=sorted(cds["db_xrefs"] + ["afdb_v6:swissprot_P1", "afdb_v6:afdbclusters_A1"]))
    cds.pop("hypothetical")
    p = tmp_path / "old.json"
    p.write_text(_json.dumps(data))
    return p, cds["locus"]


def _run(args):
    return CliRunner().invoke(main_cli, ["json", *args], catch_exceptions=False)


def test_cli_relabel_updates_outputs_and_records_provenance(tmp_path, old_genome_json):
    src, locus = old_genome_json
    out = tmp_path / "out"
    r = _run(["-i", str(src), "-o", str(out), "-p", "g", "--relabel", "-d", str(_db(tmp_path)), "-f"])
    assert r.exit_code == 0, r.output
    data = _json.loads((out / "g.json").read_text())
    cds = next(f for f in data["features"] if f.get("locus") == locus)
    assert cds["product"] == "IS66 family transposase"
    rel = data["baktfold_run"]["relabel"]
    assert rel["labels_changed"] == 1 and rel["database"]
    inference = (out / "g.inference.tsv").read_text()
    assert f"{locus}\t" in inference and "IS66 family transposase" in inference

    # a plain reconstruction of the relabelled JSON keeps the record and the labels
    out2 = tmp_path / "out2"
    assert _run(["-i", str(out / "g.json"), "-o", str(out2), "-p", "g", "-f"]).exit_code == 0
    data2 = _json.loads((out2 / "g.json").read_text())
    assert data2["baktfold_run"]["relabel"] == rel
    assert next(f for f in data2["features"] if f.get("locus") == locus)["product"] == "IS66 family transposase"


def test_cli_without_relabel_leaves_labels(tmp_path, old_genome_json):
    src, locus = old_genome_json
    out = tmp_path / "out"
    assert _run(["-i", str(src), "-o", str(out), "-p", "g", "-f"]).exit_code == 0
    cds = next(f for f in _json.loads((out / "g.json").read_text())["features"] if f.get("locus") == locus)
    assert cds["product"] == "hypothetical protein"


@pytest.mark.parametrize("extra", [["-d", "x"], ["--pdb-tophit", "x"], ["--prostt5-probs", "x"]])
def test_cli_relabel_options_need_relabel(tmp_path, old_genome_json, extra):
    src, _ = old_genome_json
    r = CliRunner().invoke(main_cli, ["json", "-i", str(src), "-o", str(tmp_path / "o"), *extra])
    assert r.exit_code == 2 and "--relabel" in r.output


def test_cli_relabel_rejects_database_without_baktfold_db(tmp_path, old_genome_json):
    src, _ = old_genome_json
    r = CliRunner().invoke(main_cli, ["json", "-i", str(src), "-o", str(tmp_path / "o"),
                                      "--relabel", "-d", str(tmp_path)])
    assert r.exit_code == 2 and "baktfold.db" in r.output


def test_relabel_keeps_curated_label_that_run_kept():
    """--all-proteins: a curated CDS whose only hit is uninformative keeps its
    input label in run (no tier). Relabel used to reset it to "hypothetical
    protein" (golden run_all MEGJMN_009, "Chromate transporter")."""
    from baktfold.subcommands.relabel import relabel_features

    curated = {
        "type": "cds", "locus": "A_1", "id": "A_1", "product": "Chromate transporter",
        "baktfold": True, "db_xrefs": [],
        "pstc": [{"source": "afdb", "id": "AF-X", "description": "hypothetical protein",
                  "evalue": 5e-43, "bitscore": 300.0, "qcov": 0.9, "tcov": 0.9, "fident": 0.4}],
    }
    tiered = {
        "type": "cds", "locus": "A_2", "id": "A_2", "product": "Transposase",
        "baktfold": True, "annotation_confidence": "high", "db_xrefs": [],
        "pstc": [{"source": "swissprot", "id": "P1", "description": "Transposase",
                  "evalue": 1e-50, "bitscore": 400.0, "qcov": 0.95, "tcov": 0.95, "fident": 0.5}],
    }
    relabel_features([curated, tiered], fast=False, has_duplicate_locus=False)
    assert curated["product"] == "Chromate transporter"
    assert "hypothetical" not in curated
