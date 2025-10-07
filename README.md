# baktfold
Rapid &amp; standardized annotation of bacterial genomes, MAGs &amp; plasmids using protein structural information

## Install

```bash
conda create -n baktfold foldseek
conda activate baktfold
pip install -e .
baktfold --help
baktfold run -i tests/assembly_bakta_output/assembly.json  -o baktfold_output -f -t 8 -d ../baktfold_db/ 
baktfold run -i tests/GCA_019351845.1_ASM1935184v1_bakta_output/GCA_019351845.1_ASM1935184v1_genomic.json  -o baktfold_output_GCA_019351845 -f -t 8 -d ../baktfold_db/ 

baktfold run -i tests/ek_isolate6_bakta_output/ek_isolate6.json  -o baktfold_output_ek_isolate6 -f -t 8 -d ../baktfold_db/ 



bakta -d ../../bakta_db/db -o ek_isolate6_bakta_output -t 4 --force ek_isolate6.fasta

```

* Where the `baktfold_db` for now is the Phold DB (for ProstT5) along with

1. Swissprot foldseek db

```bash
foldseek databases Alphafold/Swiss-Prot swissprot tmp --threads 8
```

2. AFDB 2.3M non-singleton clusters from https://www.nature.com/articles/s41586-023-06510-w https://afdb-cluster.steineggerlab.workers.dev

```bash
wget https://afdb-cluster.steineggerlab.workers.dev/1-AFDBClusters-entryId_repId_taxId.tsv.gz
zcat 1-AFDBClusters-entryId_repId_taxId.tsv.gz   |cut -f2 | uniq > AFDBCluster_reps_uniq.txt

# need to take it from full AFDB as some cluster reps are not AFDB50 member reps?

foldseek databases Alphafold/UniProt  AFDB tmp --threads 4

# https://github.com/steineggerlab/foldseek/issues/97


awk 'NR==FNR { ids[$1]=1; next }
     {
       acc=$2
       sub(/^AF-/, "", acc)
       sub(/-F[0-9]+-model_v[0-9]+$/, "", acc)
       if (acc in ids) print $1
     }' AFDBCluster_reps_uniq.txt AFDB50.lookup > subset.lookup

# does seq, ss and ca
foldseek createsubdb subset.lookup AFDB AFDBClusters
foldseek createsubdb subset.lookup AFDB_h AFDBClusters_h


# gives 9 less than AFDBCluster_reps_uniq.txt - good enough for now


# dont need this I think
# foldseek createsubdb subset.lookup AFDB_ss AFDBClusters_ss
# foldseek createsubdb subset.lookup AFDB_ca AFDBClusters_ca
rm subset.lookup

```

* Only issue that remains is the large _h and _taxonomy files - ask Milot/Martin


# TSVs

* For now, I just want the function (can add more later)
* As per Fold first ask later https://www.biorxiv.org/content/10.1101/2025.07.17.665397v1.full.pdf - AFDB/Swissport foldseek DB's out of date, so download updated functions 

```bash
awk '{IGNORECASE=1; if (match($0, />?AF-([A-Za-z0-9]+)-F/, m)) print m[1];}' AFDBClusters_h > AFDBClusters_uniprot_accessions.txt
awk '{IGNORECASE=1; if (match($0, />?AF-([A-Za-z0-9]+)-F/, m)) print m[1];}' ../baktfold_db/swissprot_h > swissprot_uniprot_accessions.txt

7 October 2025

python get_uniprot_functions.py
python get_swissprot_functions.py
```

* To get the two columns of interest

```bash
cut -f1,3 swissprot_uniprot_accessions_protein_names.tsv > ../baktfold_db/swissprot.tsv
cut -f1,3 AFDBClusters_uniprot_accessions_protein_names.tsv > ../baktfold_db/AFDBClusters.tsv

```



