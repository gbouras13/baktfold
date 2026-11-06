"""Check Baktfold's GenBank/EMBL output against the original GenBank input.

Usage: check_gbff.py <input.gbk[.gz]> <baktfold.gbff> <baktfold.embl>
"""
import gzip
import sys
from collections import Counter

from Bio import SeqIO

inp, gbff, embl = sys.argv[1:4]
opener = gzip.open if inp.endswith(".gz") else open


def records(path, fmt, op=open):
    with op(path, "rt") as fh:
        return list(SeqIO.parse(fh, fmt))


src = records(inp, "genbank", opener)
out = records(gbff, "genbank")
try:
    out_embl = records(embl, "embl")
    embl_note = f"EMBL parses ({len(out_embl)} records)"
except Exception as e:  # noqa: BLE001
    out_embl = None
    embl_note = f"EMBL PARSE ERROR: {type(e).__name__}: {e}"

TYPES = ("gene", "mRNA", "CDS", "tRNA", "rRNA", "ncRNA", "misc_RNA", "exon")
cin = Counter(f.type for r in src for f in r.features if f.type in TYPES)
cout = Counter(f.type for r in out for f in r.features if f.type in TYPES)
count_diff = {t: (cin[t], cout[t]) for t in TYPES if cin[t] != cout[t]}


def cds_key(rec, f):
    return (rec.id.split(".")[0], str(f.location), f.qualifiers.get("codon_start", ["1"])[0])


def source_table(rec):
    for f in rec.features:
        if f.type == "source" and "transl_table" in f.qualifiers:
            return f.qualifiers["transl_table"][0]
    return "1"


FLAGS = ("transl_except", "ribosomal_slippage", "pseudo", "exception", "trans_splicing")
def cds_sig(rec, f):
    q = f.qualifiers
    return (cds_key(rec, f), q.get("transl_table", [source_table(rec)])[0]) + tuple(k in q for k in FLAGS)


sig_in = Counter(cds_sig(r, f) for r in src for f in r.features if f.type == "CDS")
sig_out = Counter(cds_sig(r, f) for r in out for f in r.features if f.type == "CDS")
loc_in = Counter(s[0] for s in sig_in.elements()); loc_out = Counter(s[0] for s in sig_out.elements())
raw = open(gbff).read()
bad_text = {"/pseudo=""": raw.count('/pseudo=""'),
            "non-integer /transl_table": sum(1 for l in raw.splitlines() if "/transl_table=" in l and not l.split("=", 1)[1].strip().isdigit())}

tr_ok = tr_bad = tr_none = 0; tr_examples = []
no_locus = no_product = 0
for r in out:
    for f in r.features:
        if f.type != "CDS":
            continue
        q = f.qualifiers
        cs = q.get("codon_start", ["1"])[0]
        if "locus_tag" not in q:
            no_locus += 1
        if not q.get("product"):
            no_product += 1
        if "translation" not in q:
            tr_none += 1
            continue
        table = int(q.get("transl_table", ["1"])[0])
        nt = f.location.extract(r.seq)[int(cs) - 1:]
        prot = str(nt[: len(nt) // 3 * 3].translate(table=table)).rstrip("*")
        exp = q["translation"][0]
        if prot == exp or prot[1:] == exp[1:] or prot[1:len(exp)] == exp[1:len(prot)]:
            tr_ok += 1
        else:
            tr_bad += 1
            if len(tr_examples) < 3:
                tr_examples.append((q.get("locus_tag", ["?"])[0], prot[:20], exp[:20]))

print(f"GBFF parses ({len(out)} records); {embl_note}")
print(f"feature counts vs input: {'identical' if not count_diff else count_diff}")
print(f"CDS (location, codon_start) multiset == input: {loc_in == loc_out} (only in input {sum((loc_in - loc_out).values())}, only in output {sum((loc_out - loc_in).values())})")
print(f"CDS + transl_table + {'/'.join(FLAGS)} flags == input: {sig_in == sig_out}"
      + ("" if sig_in == sig_out else f" (e.g. input {list(sig_in - sig_out)[:1]} vs output {list(sig_out - sig_in)[:1]})"))
print(f"text problems: {bad_text}")
print(f"CDS translation from written location == /translation: {tr_ok} ok, {tr_bad} mismatch, {tr_none} without /translation")
for e in tr_examples:
    print("   e.g.", e)
print(f"CDS without /locus_tag: {no_locus}; without /product: {no_product}")
