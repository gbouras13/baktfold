#!/usr/bin/env bash
# convert-euk -> baktfold json -> GFF3 and GenBank/EMBL checks, for one input.
# usage: run_source.sh <name> <input.gbk[.gz]> <workdir>
# Needs baktfold and biopython on PATH; uses GenomeTools `gt` and `gffread`
# when they are on PATH (skipped otherwise).
set -uo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
name=$1; input=$(cd "$(dirname "$2")" && pwd)/$(basename "$2"); W=$3/$name
mkdir -p "$W" && cd "$W" || exit 1
echo "### $name"
if ! baktfold convert-euk -i "$input" -o in.json -f > conv.log 2>&1; then
  echo "convert-euk FAILED:"; grep -i -E "error|Traceback" conv.log | tail -3; exit 1
fi
python -I -c "
import json
d = json.load(open('in.json')); d.setdefault('baktfold_run', {})['euk'] = True
json.dump(d, open('in_euk.json', 'w'))"
if ! baktfold json -i in_euk.json -o out -f > json.log 2>&1; then
  echo "baktfold json FAILED:"; grep -i -E "error|Traceback" json.log | tail -3; exit 1
fi
python -I "$HERE/check_gff3.py" out/baktfold.gff3 out/baktfold.json
python -I "$HERE/check_gbff.py" "$input" out/baktfold.gbff out/baktfold.embl
if command -v gt >/dev/null; then
  gt gff3validator out/baktfold.gff3 > gt.log 2>&1 && echo "gt gff3validator: OK" \
    || { echo "gt gff3validator: $(grep -c -i error gt.log) error line(s)"; head -3 gt.log; }
fi
if command -v gffread >/dev/null; then
  awk '/^##FASTA/{f=1;next} f' out/baktfold.gff3 > genome.fa
  awk '/^##FASTA/{exit} {print}' out/baktfold.gff3 > nofasta.gff3
  gffread -y gffread_prot.faa -g genome.fa nofasta.gff3 > gffread.log 2>&1 \
    && echo "gffread: translated $(grep -c '>' gffread_prot.faa) CDS (compare with the JSON proteins if needed)" \
    || { echo "gffread FAILED"; head -3 gffread.log; }
fi
