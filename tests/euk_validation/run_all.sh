#!/usr/bin/env bash
# Run every source (downloaded by fetch_sources.sh, plus the repo's own test
# genomes) and print one summary.  usage: run_all.sh <sources_dir> <workdir>
set -uo pipefail
HERE=$(cd "$(dirname "$0")" && pwd); S=${1:?sources dir}; W=${2:?workdir}
T=$HERE/../test_data
for pair in \
  yeast_refseq:$S/yeast_refseq.gbff.gz dmel_chr4_refseq:$S/dmel_chr4_refseq.gbk \
  pfal_chr1_refseq:$S/pfal_chr1_refseq.gbk celegans_chrI_ensembl:$S/celegans_chrI_ensembl.dat.gz \
  ostreococcus_chr1_ensemblplants:$S/ostreococcus_chr1_ensemblplants.dat.gz \
  pneumocystis_genbank:$S/pneumocystis_genbank.gbff malassezia_genbank:$S/malassezia_genbank.gbff \
  ecuniculi_genbank:$S/ecuniculi_genbank.gbff \
  aspergillus_funannotate_ncbi:$T/euks/aspergillus_subset.gbk entamoeba_insdc:$T/protist.gbk.gz \
  mucor_funannotate:$T/funannotate.gbk.gz; do
  f=${pair#*:}
  # aspergillus_subset.gbk is not committed to the repo; skip it if absent
  [ -e "$f" ] || { echo "### ${pair%%:*}: input not found ($f), skipped"; continue; }
  bash "$HERE/run_source.sh" "${pair%%:*}" "$f" "$W"
done
