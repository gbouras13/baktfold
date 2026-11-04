"""Old vs new PDB / CATH labels, for review before a database release.

Writes two CSVs:

  pdb_labels_old_vs_new.csv   one row per chain in the Foldseek PDB database
      entry, chain, old_title (entry title, the label used until now),
      old_label (as Baktfold wrote it), new_name (pdb_chain table), new_source
      (uniprot / pdb_entity / none), new_label (as Baktfold will write it),
      uniprot_accession (SIFTS), changed
  cath_labels_old_vs_new.csv  one row per distinct CATH superfamily name
      cath_name, n_domains, old_label, new_label, skipped_as_uninformative

"as Baktfold writes it" = after the hierarchy's uninformative-hit rule and
Bakta's product cleanup (revise_cds_product); an empty label means the hit is
skipped and the next database in the hierarchy is used.

Usage:
  python compare-pdb-cath-labels.py --db baktfold_with_pdb_chain.db \
      --lookup pdb_gpu.lookup --sifts pdb_chain_uniprot.tsv.gz \
      --baktfold-src /path/to/baktfold/src --outdir .
"""
import argparse
import csv
import os
import re
import sqlite3
import sys
import types
from collections import Counter, defaultdict
from pathlib import Path

from xopen import xopen

parser = argparse.ArgumentParser()
parser.add_argument("--db", required=True, help="baktfold.db that already has the pdb_chain table")
parser.add_argument("--lookup", required=True, help="Foldseek PDB .lookup file")
parser.add_argument("--sifts", required=True, help="pdb_chain_uniprot.tsv(.gz)")
parser.add_argument("--baktfold-src", required=True, help="baktfold src/ dir (label rules are imported from it)")
parser.add_argument("--outdir", default=".")
args = parser.parse_args()

# import only baktfold's label code (the package __init__ pulls in the CLI)
pkg = types.ModuleType("baktfold"); pkg.__path__ = [os.path.join(args.baktfold_src, "baktfold")]
sys.modules["baktfold"] = pkg
import baktfold.bakta.annotation as anno  # noqa: E402
from loguru import logger  # noqa: E402
logger.remove()


def as_written(desc: str, cath: bool = False) -> str:
    """Label Baktfold writes for a hit with this description ('' = skipped)."""
    if not anno.is_informative_product(desc):
        return ""
    if cath:
        desc = anno.cath_domain_label(desc)
    return anno.revise_cds_product(desc) or ""


def old_written(desc: str) -> str:
    """Label under the previous rule: top hit taken verbatim, then cleanup."""
    if not desc or "hypothetical protein" in desc.lower():
        return ""
    r = anno.revise_cds_product(desc)
    return "" if (not r or "hypothetical protein" in r.lower()) else r


out = Path(args.outdir)
conn = sqlite3.connect(f"file:{args.db}?mode=ro", uri=True)

# ---------------------------------------------------------------- PDB
titles = dict(conn.execute("SELECT id, product FROM pdb"))
chain_names = {r[0]: (r[1], r[2]) for r in conn.execute("SELECT id, product, source FROM pdb_chain")}
span = defaultdict(Counter)
with xopen(args.sifts) as fh:
    header = None
    for line in fh:
        if line.startswith("#"):
            continue
        c = line.rstrip("\n").split("\t")
        if header is None:
            header = {k: i for i, k in enumerate(c)}; continue
        try:
            n = int(c[header["SP_END"]]) - int(c[header["SP_BEG"]]) + 1
        except ValueError:
            n = 1
        span[f"{c[header['PDB']].lower()}_{c[header['CHAIN']]}"][c[header["SP_PRIMARY"]]] += max(n, 1)

chains = set()
with open(args.lookup) as fh:
    for line in fh:
        target = line.split("\t")[1]
        entry = target.split("-")[0]
        chains.add((entry, re.sub(r"-\d+$", "", target.rsplit("_", 1)[1])))

tally = Counter()
with (out / "pdb_labels_old_vs_new.csv").open("w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["entry", "chain", "old_title", "old_label", "new_name", "new_source",
                "new_label", "uniprot_accession", "changed"])
    for entry, chain in sorted(chains):
        key = f"{entry}_{chain}"
        title = titles.get(entry, "")
        name, source = chain_names.get(key, ("", "none"))[:2]
        old_l, new_l = old_written(title), as_written(name)
        acc = span[key].most_common(1)[0][0] if key in span else ""
        w.writerow([entry, chain, title, old_l, name, source, new_l, acc, int(old_l != new_l)])
        tally[source] += 1
        tally["old_labelled"] += bool(old_l); tally["new_labelled"] += bool(new_l)
print("PDB chains:", len(chains), dict(tally))

# ---------------------------------------------------------------- CATH
n_dom = Counter(r[0] for r in conn.execute("SELECT product FROM cath") if r[0])
with (out / "cath_labels_old_vs_new.csv").open("w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["cath_name", "n_domains", "old_label", "new_label", "skipped_as_uninformative"])
    for name, n in sorted(n_dom.items(), key=lambda x: -x[1]):
        new_l = as_written(name, cath=True)
        w.writerow([name, n, old_written(name), new_l, int(not new_l)])
print("CATH names:", len(n_dom))
