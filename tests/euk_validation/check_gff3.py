"""Structural + functional check of a Baktfold eukaryotic GFF3.

Usage: check_gff3.py <baktfold.gff3> <baktfold.json>
"""
import json
import sys
from collections import Counter, defaultdict
from urllib.parse import unquote

from Bio.Seq import Seq

gff, js = sys.argv[1], sys.argv[2]
rows, seqs, cur, in_fasta = [], {}, None, False
for line in open(gff):
    line = line.rstrip("\n")
    if line.startswith("##FASTA"):
        in_fasta = True
        continue
    if in_fasta:
        if line.startswith(">"):
            cur = line[1:].split()[0]; seqs[cur] = []
        elif cur:
            seqs[cur].append(line)
        continue
    if not line or line.startswith("#"):
        continue
    c = line.split("\t")
    if len(c) != 9:
        continue
    attrs = dict(kv.split("=", 1) for kv in c[8].split(";") if "=" in kv)
    rows.append(dict(seqid=c[0], type=c[2], start=int(c[3]), end=int(c[4]), strand=c[6],
                     phase=c[7], ID=unquote(attrs.get("ID", "")),
                     Parent=unquote(attrs.get("Parent", "")), locus=attrs.get("locus_tag")))
seqs = {k: "".join(v) for k, v in seqs.items()}

problems = Counter(); examples = defaultdict(list)
def bad(kind, msg):
    problems[kind] += 1
    if len(examples[kind]) < 3:
        examples[kind].append(msg)

by_id = defaultdict(list)
for r in rows:
    if r["ID"]:
        by_id[r["ID"]].append(r)
types = {i: {r["type"] for r in rs} for i, rs in by_id.items()}
for i, rs in by_id.items():
    if len(rs) > 1:
        if {r["type"] for r in rs} != {"CDS"}:
            bad("duplicate ID across non-CDS features", i)
        elif len({r["Parent"] for r in rs}) != 1:
            bad("CDS segments with one ID but different Parents", i)
for r in rows:
    if r["Parent"]:
        for p in r["Parent"].split(","):
            if p not in by_id:
                bad(f"{r['type']} Parent missing", f"{r['ID']} -> {p}")
            elif r["type"] in ("CDS", "exon") and not (types[p] & {"mRNA", "gene", "tRNA", "ncRNA", "rRNA", "transcript"}):
                bad(f"{r['type']} Parent not mRNA/gene", f"{r['ID']} -> {p} {types[p]}")
            elif r["type"] == "mRNA" and "gene" not in types[p]:
                bad("mRNA Parent not gene", f"{r['ID']} -> {p}")

# translate each CDS from the GFF3 alone and compare with the JSON protein
data = json.load(open(js))
aa_of = defaultdict(set); table_of = {}
for f in data["features"]:
    if f.get("type") == "cds" and f.get("aa"):
        aa_of[f.get("locus")].add(f["aa"].rstrip("*"))
        if f.get("transl_table"):
            table_of[f.get("locus")] = int(f["transl_table"])
cds_groups = defaultdict(list)
for r in rows:
    if r["type"] == "CDS":
        cds_groups[r["ID"]].append(r)
tr_ok = tr_bad = tr_skip = 0
for cid, segs in cds_groups.items():
    strand = segs[0]["strand"]
    segs = sorted(segs, key=lambda r: r["start"])  # genomic order; RC below
    s = seqs.get(segs[0]["seqid"])
    locus = segs[0]["locus"]
    if not s or locus not in aa_of:
        tr_skip += 1; continue
    nt = "".join(s[r["start"] - 1:r["end"]] for r in segs)
    nt = str(Seq(nt).reverse_complement()) if strand == "-" else nt
    first = segs[-1] if strand == "-" else segs[0]   # 5'-most segment
    nt = nt[int(first["phase"]):]
    prot = str(Seq(nt[:len(nt) // 3 * 3]).translate(table=table_of.get(locus, 1))).rstrip("*")
    # partial / non-ATG starts and terminal partial codons: compare cores
    if any(prot == a or prot[1:] == a[1:] or prot[1:len(a)] == a[1:len(prot)] for a in aa_of[locus]):
        tr_ok += 1
    else:
        tr_bad += 1
        if tr_bad <= 3:
            bad("translation mismatch", f"{cid}: {prot[:25]}… vs {next(iter(aa_of[locus]))[:25]}…")

n = Counter(r["type"] for r in rows)
print(f"features: gene {n['gene']}, mRNA {n['mRNA']}, CDS lines {n['CDS']} in {len(cds_groups)} CDS, exon {n['exon']}")
print(f"CDS translation from GFF3 matches JSON protein: {tr_ok} ok, {tr_bad} mismatch, {tr_skip} not checked")
print("structural problems:", dict(problems) if problems else "none")
for k, v in examples.items():
    print(f"   {k}: {v}")
