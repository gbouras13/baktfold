#!/usr/bin/env bash
# Download the public eukaryotic GenBank/EMBL flatfiles used by run_all.sh.
# usage: fetch_sources.sh <dir>
set -euo pipefail
D=${1:?usage: fetch_sources.sh <dir>}
mkdir -p "$D"
cd "$D"

# NCBI RefSeq: S. cerevisiae S288C (whole genome, incl. mitochondrion)
curl -sSfL -o yeast_refseq.gbff.gz \
  https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.gbff.gz

# NCBI RefSeq, one chromosome each: D. melanogaster chr4, P. falciparum chr1
for pair in dmel_chr4_refseq:NC_004353.4 pfal_chr1_refseq:NC_004325.2; do
  curl -sSfL -o "${pair%%:*}.gbk" \
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${pair##*:}&rettype=gbwithparts&retmode=text"
done

# Ensembl: C. elegans chrI; Ensembl Plants: Ostreococcus lucimarinus chr1
E=https://ftp.ensembl.org/pub/current/genbank/caenorhabditis_elegans
curl -sSfL -o celegans_chrI_ensembl.dat.gz "$E/$(curl -sSfL "$E/" | grep -o 'Caenorhabditis_elegans\.WBcel235\.[0-9]*\.chromosome\.I\.dat\.gz' | head -1)"
P=https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/current/genbank/ostreococcus_lucimarinus
curl -sSfL -o ostreococcus_chr1_ensemblplants.dat.gz "$P/$(curl -sSfL "$P/" | grep -o 'Ostreococcus_lucimarinus\.ASM9206v1\.[0-9]*\.chromosome\.1\.dat\.gz' | head -1)"

# INSDC (GenBank) assemblies from three submitters, via NCBI Datasets
for pair in pneumocystis_genbank:GCA_017301755.1 malassezia_genbank:GCA_029542905.1 ecuniculi_genbank:GCA_000091225.2; do
  n=${pair%%:*}; a=${pair##*:}
  curl -sSfL -o "$a.zip" "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/$a/download?include_annotation_type=GENOME_GBFF"
  unzip -p "$a.zip" "ncbi_dataset/data/$a/genomic.gbff" > "$n.gbff"
  rm "$a.zip"
done
ls -l "$D"
