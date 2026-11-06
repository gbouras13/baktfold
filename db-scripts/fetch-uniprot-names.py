"""Fetch UniProt protein names for a list of accessions (UniProt REST API).

Used by init-pdb-chains.py for SIFTS-mapped PDB chains whose UniProt accession
is not already in the Baktfold database (swissprot / afdbclusters tables).

Writes ``accession<TAB>name`` using the recommended full name (Swiss-Prot) or
the first submission name (TrEMBL). Accessions UniProt no longer returns
(obsolete / merged) are written with an empty name, so a rerun resumes rather
than re-querying them.

Usage:
    python fetch-uniprot-names.py --accessions accs.txt --out uniprot_names.tsv
"""
import argparse
import json
import time
import urllib.error
import urllib.request
from pathlib import Path

API = "https://rest.uniprot.org/uniprotkb/accessions?accessions={}&fields=accession,protein_name&format=json"
BATCH = 500

parser = argparse.ArgumentParser(description="Fetch UniProt protein names.")
parser.add_argument("--accessions", required=True, help="One accession per line.")
parser.add_argument("--out", required=True, help="Output TSV (appended; resumable).")
args = parser.parse_args()


def protein_name(entry: dict) -> str:
    desc = entry.get("proteinDescription", {})
    rec = desc.get("recommendedName")
    if rec:
        return rec["fullName"]["value"]
    sub = desc.get("submissionNames") or []
    return sub[0]["fullName"]["value"] if sub else ""


def fetch(batch):
    url = API.format(",".join(batch))
    for attempt in range(5):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                return json.load(r).get("results", [])
        except (urllib.error.URLError, TimeoutError) as e:
            wait = 10 * (attempt + 1)
            print(f"  retry in {wait}s ({e})")
            time.sleep(wait)
    raise RuntimeError(f"UniProt query failed repeatedly for batch starting {batch[0]}")


wanted = [l.strip() for l in open(args.accessions) if l.strip()]
out = Path(args.out)
done = set()
if out.exists():
    done = {l.split("\t", 1)[0] for l in out.read_text().splitlines() if l}
todo = [a for a in dict.fromkeys(wanted) if a not in done]
print(f"{len(wanted):,} accessions, {len(done):,} already fetched, {len(todo):,} to fetch")

with out.open("a") as fh:
    for i in range(0, len(todo), BATCH):
        batch = todo[i:i + BATCH]
        got = {e["primaryAccession"]: protein_name(e) for e in fetch(batch)}
        for acc in batch:
            fh.write(f"{acc}\t{got.get(acc, '')}\n")
        fh.flush()
        print(f"  {min(i + BATCH, len(todo)):,}/{len(todo):,}")
        time.sleep(0.5)  # be polite to the public API
