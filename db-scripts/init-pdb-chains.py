"""Build the chain-level PDB name table (``pdb_chain``) in baktfold.db.

Foldseek PDB hits are chains (``6xgw-assembly1.cif.gz_A``). The original
``pdb`` table names whole entries by their title ("Crystal structure of ... in
complex with ..."), which is not a protein name. This table names each chain:

  1. UniProt name of the chain via SIFTS (pdb_chain_uniprot.tsv), taken from
     the swissprot / afdbclusters tables already in baktfold.db or, failing
     that, from --uniprot-names (see fetch-uniprot-names.py);          source=uniprot
  2. else the PDB entity (molecule) name from pdb_seqres.txt;          source=pdb_entity
  3. else no row: Baktfold then treats the hit as unlabelled and moves on
     to the next database (it never falls back to the entry title).

Names that are uninformative under Baktfold's rule (mention hypothetical /
uncharacterised) are skipped at each step, so step 2 can still apply.

Inputs (wwPDB and EBI SIFTS, both public):
  https://files.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt
  https://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/tsv/pdb_chain_uniprot.tsv.gz

Two-step usage (step 1 writes the accessions that still need names):
  python init-pdb-chains.py --seqres pdb_seqres.txt --sifts pdb_chain_uniprot.tsv.gz \
      --db baktfold.db --lookup pdb_gpu.lookup --missing-out missing_accessions.txt --dry-run
  python fetch-uniprot-names.py --accessions missing_accessions.txt --out uniprot_names.tsv
  python init-pdb-chains.py --seqres pdb_seqres.txt --sifts pdb_chain_uniprot.tsv.gz \
      --db baktfold.db --lookup pdb_gpu.lookup --uniprot-names uniprot_names.tsv
"""
import argparse
import re
import sqlite3
from collections import Counter, defaultdict
from pathlib import Path

from xopen import xopen

parser = argparse.ArgumentParser(description="Create the pdb_chain table in baktfold.db.")
parser.add_argument("--seqres", required=True, help="pdb_seqres.txt (optionally gzipped)")
parser.add_argument("--sifts", required=True, help="pdb_chain_uniprot.tsv(.gz) from SIFTS")
parser.add_argument("--db", required=True, help="baktfold.db (modified in place unless --dry-run)")
parser.add_argument("--lookup", help="Foldseek PDB .lookup file: restrict to, and report coverage of, chains actually searched")
parser.add_argument("--uniprot-names", help="accession<TAB>name TSV from fetch-uniprot-names.py")
parser.add_argument("--missing-out", help="Write SIFTS accessions with no name available yet")
parser.add_argument("--dry-run", action="store_true", help="Report only; do not write the table")
args = parser.parse_args()

# Baktfold's skip rule (bakta/annotation.py RE_UNINFORMATIVE_PRODUCT)
RE_UNINFORMATIVE = re.compile(r"\b(hypothetical|uncharacteri[sz]ed)\b", re.IGNORECASE)


def informative(name: str) -> bool:
    return bool(name) and bool(name.strip()) and not RE_UNINFORMATIVE.search(name)


# short words that are ordinary words, not acronyms, in all-caps PDB names
SHORT_WORDS = {"alpha", "beta", "gamma", "delta", "zeta", "theta", "kappa", "sigma", "omega",
               "chain", "like", "type", "unit", "form", "and", "of", "the", "for", "with",
               "from", "to", "in", "on", "by", "its", "acid", "rich", "site", "head", "tail",
               "half", "core", "heme", "zinc", "iron", "beta", "light", "heavy"}


def recase(name: str) -> str:
    """Old PDB entity names are all caps ("DNA POLYMERASE III"). Lower-case
    long alphabetic words and keep short tokens / tokens with digits, which are
    usually acronyms or gene names: -> "DNA polymerase III"."""
    if not name.isupper():
        return name
    words = [
        w.lower() if (w.isalpha() and (len(w) > 4 or w.lower() in SHORT_WORDS)) else w
        for w in name.split(" ")
    ]
    out = " ".join(words)
    return out[:1].upper() + out[1:]


def chain_key(entry: str, chain: str) -> str:
    return f"{entry.lower()}_{chain}"


# --- chains actually in the Foldseek PDB database --------------------------
searched = None
if args.lookup:
    searched = set()
    with open(args.lookup) as fh:
        for line in fh:
            target = line.split("\t")[1]
            entry = target.split("-")[0]
            chain = re.sub(r"-\d+$", "", target.rsplit("_", 1)[1])
            searched.add(chain_key(entry, chain))
    print(f"Foldseek PDB database: {len(searched):,} distinct chains")

# --- PDB entity names (pdb_seqres.txt) --------------------------------------
entity = {}
with xopen(args.seqres) as fh:
    for line in fh:
        if not line.startswith(">"):
            continue
        # >101m_A mol:protein length:154  MYOGLOBIN
        head, _, rest = line[1:].rstrip("\n").partition(" ")
        if "mol:protein" not in rest:
            continue
        entry, _, chain = head.partition("_")
        name = re.sub(r"^mol:\S+\s+length:\d+\s*", "", rest).strip()
        entity[chain_key(entry, chain)] = recase(name)
print(f"pdb_seqres: {len(entity):,} protein chains")

# --- SIFTS chain -> UniProt (largest mapped span wins) ----------------------
span = defaultdict(Counter)
with xopen(args.sifts) as fh:
    header = None
    for line in fh:
        if line.startswith("#"):
            continue
        cols = line.rstrip("\n").split("\t")
        if header is None:
            header = {c: i for i, c in enumerate(cols)}
            continue
        try:
            n = int(cols[header["SP_END"]]) - int(cols[header["SP_BEG"]]) + 1
        except ValueError:
            n = 1
        span[chain_key(cols[header["PDB"]], cols[header["CHAIN"]])][cols[header["SP_PRIMARY"]]] += max(n, 1)
sifts = {k: c.most_common(1)[0][0] for k, c in span.items()}
print(f"SIFTS: {len(sifts):,} chains mapped to UniProt")

# --- UniProt names -----------------------------------------------------------
names = {}
with sqlite3.connect(f"file:{args.db}?mode=ro", uri=True) as ro:
    for table in ("afdbclusters", "swissprot"):  # swissprot last: curated names win
        for acc, product in ro.execute(f"SELECT id, product FROM {table}"):
            if product:
                names[acc] = product
if args.uniprot_names:
    with open(args.uniprot_names) as fh:
        for line in fh:
            acc, _, name = line.rstrip("\n").partition("\t")
            if name and acc not in names:
                names[acc] = name

keys = searched if searched is not None else set(entity) | set(sifts)
missing = sorted({sifts[k] for k in keys if k in sifts and sifts[k] not in names})
if args.missing_out:
    Path(args.missing_out).write_text("".join(f"{a}\n" for a in missing))
print(f"UniProt accessions still without a name: {len(missing):,}")

# --- resolve -----------------------------------------------------------------
rows, tally = [], Counter()
for k in sorted(keys):
    acc = sifts.get(k)
    if acc and informative(names.get(acc, "")):
        rows.append((k, names[acc], "uniprot")); tally["uniprot"] += 1
    elif informative(entity.get(k, "")):
        rows.append((k, entity[k], "pdb_entity")); tally["pdb_entity"] += 1
    else:
        tally["no_name" if (k in entity or k in sifts) else "not_in_seqres_or_sifts"] += 1
total = len(keys)
for label in ("uniprot", "pdb_entity", "no_name", "not_in_seqres_or_sifts"):
    print(f"  {label:24} {tally[label]:>9,}  ({tally[label] / max(total, 1):.1%})")

if args.dry_run:
    print("dry run: table not written")
else:
    with sqlite3.connect(args.db) as conn:
        conn.execute("DROP TABLE IF EXISTS pdb_chain")
        conn.execute("CREATE TABLE pdb_chain (id TEXT PRIMARY KEY, product TEXT, source TEXT)")
        conn.executemany("INSERT INTO pdb_chain VALUES (?, ?, ?)", rows)
        conn.commit()
    print(f"wrote {len(rows):,} rows to pdb_chain in {args.db}")
