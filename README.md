# baktfold
Rapid &amp; standardized annotation of bacterial genomes, MAGs &amp; plasmids using protein structural information

## Install

```bash
conda create -n baktfold foldseek
conda activate baktfold
pip install -e .
baktfold --help
```

## Test

```bash
# with nvidia gpu
baktfold run -i tests/assembly_bakta_output/assembly.json  -o baktfold_output -f -t 8 -d ../baktfold_db/   --foldseek_gpu
# without
baktfold run -i tests/assembly_bakta_output/assembly.json  -o baktfold_output -f -t 8 -d ../baktfold_db/   


baktfold run -i S1C2_bakta_1_11_4/S1C2.json  -o baktfold_S1C2 -f -t 8 -d ../baktfold_db/   --foldseek_gpu
bakta -d ../bakta_db/db -o S1C2_bakta_1_11_4 -t 8 --force S1C2/S1C2.fna
```


## Notes


```bash
conda create -n baktfold foldseek
conda activate baktfold
pip install -e .
baktfold --help
baktfold run -i tests/assembly_bakta_output/assembly.json  -o baktfold_output -f -t 8 -d ../baktfold_db/  
baktfold run -i tests/assembly_bakta_output/assembly.json  -o baktfold_output -f -t 8 -d ../baktfold_db/  --foldseek_gpu


baktfold run -i tests/assembly_bakta_output/assembly.json  -o baktfold_output -f -t 8 -d ../baktfold_db/   --custom_db ../baktfold_db/pdb
baktfold run -i tests/assembly_bakta_output/assembly.json  -o baktfold_output -f -t 8 -d ../baktfold_db/   --custom_db tests/custom_db/dummy_custom_db --custom_annotations tests/dummy_custom_db_annotations.tsv

baktfold run -i tests/GCA_019351845.1_ASM1935184v1_bakta_output/GCA_019351845.1_ASM1935184v1_genomic.json  -o baktfold_output_GCA_019351845 -f -t 8 -d ../baktfold_db/ 

baktfold run -i tests/ek_isolate6_bakta_output/ek_isolate6.json  -o baktfold_output_ek_isolate6 -f -t 8 -d ../baktfold_db/ 

# need to pad the db first 
foldseek makepaddedseqdb pdb pdb_gpu
foldseek makepaddedseqdb swissprot swissprot_gpu
foldseek makepaddedseqdb AFDBClusters AFDBClusters_gpu 
foldseek makepaddedseqdb cath cath_gpu 


    foldseek_makepaddedseqdb = ExternalTool(
        tool="foldseek",
        input=f"",
        output=f"",
        params=f"makepaddedseqdb {baktfold_db_search} {baktfold_db_search_gpu}",
        logdir=logdir,
    )

baktfold run -i tests/ek_isolate6_bakta_output/ek_isolate6.json  -o baktfold_output_ek_isolate6 -f -t 8 -d ../baktfold_db/ --foldseek_gpu

baktfold run -i tests/ek_isolate6_bakta_output/ek_isolate6.json  -o baktfold_output_ek_isolate6 -f -t 8 -d ../baktfold_db/ --foldseek_gpu 


bakta -d ../../bakta_db/db -o ek_isolate6_bakta_output -t 4 --force ek_isolate6.fasta
bakta -d ../../bakta_db/db -o bas_sterne_bakta_output -t 4 --force GCA_019351845.1_ASM1935184v1_genomic.fna

baktfold run -i tests/bas_sterne_bakta_output/GCA_019351845.1_ASM1935184v1_genomic.json  -o baktfold_output_bas_sterne -f -t 8 -d ../baktfold_db/ --foldseek_gpu 



# proteins
bakta_proteins -d ../../bakta_db/db -o assembly_bakta_proteins_output -t 8 --force assembly.hypotheticals.faa

baktfold proteins -i tests/assembly.hypotheticals.faa  -o baktfold_proteins_output -f -t 8 -d ../baktfold_db/ 
baktfold proteins-predict -i tests/assembly.hypotheticals.faa  -o baktfold_proteins_predict_output -f -t 8 -d ../baktfold_db/ 
baktfold proteins-compare -i tests/assembly.hypotheticals.faa  --predictions_dir baktfold_proteins_predict_output -o baktfold_proteins_compare_output -f -t 8 -d ../baktfold_db/ 


# compare

# using existing ProstT5 predictions (will be fine one I add predict)
baktfold compare -i tests/assembly_bakta_output/assembly.json --predictions_dir baktfold_output  -o baktfold_output_compare -f -t 8 -d ../baktfold_db/ 

# using pdbs
baktfold compare -i tests/assembly_bakta_output/assembly.json --structure_dir tests/pdbs  -o baktfold_output_compare -f -t 8 -d ../baktfold_db/ 
baktfold compare -i tests/assembly_bakta_output/assembly.json --structure_dir tests/cifs  -o baktfold_output_compare -f -t 8 -d ../baktfold_db/ 

# predict

baktfold predict -i tests/assembly_bakta_output/assembly.json  -o baktfold_output_predict -f -t 8 -d ../baktfold_db/ 
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


# 9 October 

* New v6 AFDB and Swissprot DBs available synced to Swissprot v2025_03
* Therefore, use those headers and tags as they are updated - no need to query uniprot API
* Only question now is using AFDB50 clusters - redo clustering as fastest 
* with Foldseek v `8979d230fb64c7089380b652758d8705493ed4a5` on the AFDB50 minimal
* Cluster at 90% coverage E=0.01 as per nature paper https://www.nature.com/articles/s41586-023-06510-w

```bash
./foldseek/bin/foldseek databases Alphafold/UniProt50-minimal afdb50 tmp --threads 8

wc -l afdb50_h # 66725340


# foldseek cluster at the same thresholds as nature paper
./foldseek/bin/foldseek cluster afdb50  afdb50_clusterDB tmp -c 0.9 -e 0.01

# gets the tsv
./foldseek/bin/foldseek createtsv afdb50  afdb50  afdb50_cluster_2itsDB cluster_2its.tsv


cut -f1 cluster_2its.tsv | uniq  | wc -l     
32967214
# total number of singletons
cut -f1 cluster_2its.tsv | uniq -u | wc -l
28909850

# fine to just take the non singletons I guess
cut -f1 cluster_2its.tsv | uniq -u > singletons.txt

# 
awk '!a[$1]++ {print $2}' cluster_2its.tsv > representatives.txt


grep -v -F -f singletons.txt representatives.txt > cluster_reps_2plus.txt
wc -l cluster_reps_2plus.txt
# 4057364


# 

awk 'NR==FNR { ids[$1]=1; next }
     {
       acc=$2
       if (acc in ids) print $1
     }' cluster_reps_2plus.txt afdb50.lookup > subset.lookup




# does seq, ss and ca
foldseek createsubdb subset.lookup afdb50 AFDBClustersv6
foldseek createsubdb subset.lookup afdb50_h AFDBClustersv6_h

awk '{match($1, /AF-([A-Z0-9]+)-F1-model_v6/, m); if (m[1] != "") print m[1] "\t" substr($0, index($0,$2))}' AFDBClustersv6_h > AFDBClustersv6.tsv


# then to get the Fragment

grep  "Fragment" AFDBClustersv6.tsv | wc -l 
# 971006

grep -v "Fragment" AFDBClustersv6.tsv | cut -f1 > cluster_reps_2plus_no_Fragment.txt

wc -l cluster_reps_2plus_no_Fragment.txt
# 3085778

# need to account for the full AF string

awk 'NR==FNR { ids[$1]=1; next }
     {
       acc=$2
       sub(/^AF-/, "", acc)
       sub(/-F[0-9]+-model_v[0-9]+$/, "", acc)
       if (acc in ids) print $1
     }' cluster_reps_2plus_no_Fragment.txt afdb50.lookup > subset.lookup

# does seq, ss and ca
foldseek createsubdb subset.lookup afdb50 AFDBClustersv6
foldseek createsubdb subset.lookup afdb50_h AFDBClustersv6_h

awk '{match($1, /AF-([A-Z0-9]+)-F1-model_v6/, m); if (m[1] != "") print m[1] "\t" substr($0, index($0,$2))}' AFDBClustersv6_h > AFDBClustersv6.tsv

# in the database directory
for f in AFDBClustersv6*; do
    mv "$f" "${f/v6/}"
done




# 32967214 reps -> try cluster some more with 2+2 its

# for the extra clustering
awk 'NR==FNR { ids[$1]=1; next }
     {
       acc=$2
       if (acc in ids) print $1
     }' representatives.txt afdb50.lookup > subset.lookup

# does seq, ss and ca
foldseek createsubdb subset.lookup afdb50 afdb50_2its
foldseek createsubdb subset.lookup afdb50_h afdb50_2its_h

# then cluster this some more




```

* Swissprot

```bash
foldseek databases Alphafold/Swiss-Prot swissprot tmp --threads 8
awk '{match($1, /AF-([A-Z0-9]+)-F1-model_v6/, m); if (m[1] != "") print m[1] "\t" substr($0, index($0,$2))}' swissprot_h > swissprot.tsv
```

* Pdb 
* Note only using the  reps for now as it was clustered at 100% seqid 95% coverage https://github.com/steineggerlab/foldseek/issues/258

```bash
# Clustering: Create 'targetDB_clu100' by removing redundancies from 'targetDB'
MMSEQS_FORCE_MERGE=1 $foldseek cluster pdb_seq pdb_clu tmp -c 0.95 --min-seq-id 1.0 --cov-mode 0
# Create sub-database 'pdb'
$foldseek createsubdb pdb_clu pdb_seq pdb --subdb-mode 1
```

* To get the descriptions (deduplicated as well)

```bash
foldseek databases PDB pdb tmp --threads 8
awk '{split($1, a, "-"); $1=""; print a[1] "\t" substr($0, 2)}' pdb_h > pdb.tsv
tr -cd '\11\12\15\40-\176' < pdb.tsv > pdb_clean.tsv
awk '!seen[$0]++' pdb_clean.tsv > pdb_unique.tsv
mv pdb_unique.tsv pdb.tsv
rm pdb_clean.tsv
```

* Swissprot

```bash
foldseek databases Alphafold/Swiss-Prot swissprot tmp --threads 8
awk '{match($1, /AF-([A-Z0-9]+)-F1-model_v6/, m); if (m[1] != "") print m[1] "\t" substr($0, index($0,$2))}' swissprot_h > swissprot.tsv
```


## 13 October Custom DB

* Test custom DB consists of the files in `tests/pdbs/`

```bash 
# from tests/
mkdir -p custom_db
cd custom_db
foldseek createdb ../pdbs dummy_custom_db
head dummy_custom_db_h
# MEGJMN_003
# MEGJMN_005
```

* Then create a dummy custom annotations csv


## 20 Oct

* Get a list of the AFDB Cluster accessions
* Next filter the enormous Uniref100 xml to keep only the AFDB Cluster accessions for parsing later

```bash
 cut -f1 AFDBClusters.tsv > AFDBClusters_accessions.tsv
 cut -f1 swissprot.tsv > swissprot_accessions.tsv

python ../baktfold/filter_xml.py --accessions swissprot_accessions.tsv --xml uniprot_2025_03/uniprot_sprot.xml.gz --output uniprot_2025_03/uniprot_sprot_test.xml
python ../baktfold/filter_xml.py --accessions AFDBClusters_accessions.tsv --xml uniprot_2025_03/uniprot_trembl.xml.gz --output uniprot_2025_03/uniprot_trembl_AFDBClusters.xml
```

```bash
python ../baktfold/filter_xml.py --accessions swissprot_accessions.tsv --xml uniprot_2025_03/uniprot_sprot.xml.gz --output uniprot_2025_03/uniprot_sprot_test.xml --output_tsv found_sprot_accessions.tsv

python ../baktfold/filter_xml.py --accessions AFDBClusters_accessions.tsv --xml uniprot_2025_03/uniprot_trembl.xml.gz --output uniprot_2025_03/uniprot_trembl_AFDBClusters.xml --output_tsv uniprot_2025_03/found_AFDBClusters_accessions.tsv

wc -l .
3085778

wc -l found_AFDBClusters_accessions.tsv
3084296
```

* Inspect 1482 or so accessions that didn't map

```bash
grep -Fxv -f uniprot_2025_03/found_AFDBClusters_accessions.tsv AFDBClusters_accessions.tsv   > uniprot_2025_03/missing_AFDBCluster_accessions.tsv
```

* These seem to be found in Swissprot but not trembl

```bash
zgrep Q9S747 uniprot_2025_03/uniprot_sprot.xml.gz
    <accession>Q9S747</accession>
    <dbReference type="AlphaFoldDB" id="Q9S747"/>
    <dbReference type="SMR" id="Q9S747"/>
```

* So to get these

```bash
python ../baktfold/filter_xml.py --accessions uniprot_2025_03/missing_AFDBCluster_accessions.tsv --xml uniprot_2025_03/uniprot_sprot.xml.gz --output uniprot_2025_03/missing_AFDBCluster_sprot.xml --output_tsv missing_AFDBCluster_accessions_found_swissprot.tsv
✅ Done. 1,482 matching entries written to /scratchdata1/users/a1667917/baktfold_db/uniprot_2025_03/missing_AFDBCluster_sprot.xml
✅ 1,482 accessions written to /scratchdata1/users/a1667917/baktfold_db/missing_AFDBCluster_accessions_found_swissprot.tsv
```



## 22 October CATH

* Decide to add CATH50 (after nice experiments with Foldseek)

```bash
foldseek databases CATH50 CATH50 tmp --threads 4
```

* This contains approx 195k entries. Approx 135k are from AFDB entries, with the CATH classification code helpfully in the foldseek header
* e.g.

```bash
af_A0A023GPK8_21_130_2.60.40.10
cat CATH50_h | grep -a af_ | wc -l
134491
```

* The rest are the (PDB) CATH domain only
* e.g.

```bash
9antA00
```

* For the PDB ones, we can match the domain to the CATH classification code using

```
wget --no-check-certificate https://download.cathdb.info/cath/releases/all-releases/v4_4_0/cath-classification-data/cath-domain-list-v4_4_0.txt

head -n20 cath-domain-list-v4_4_0.txt

#---------------------------------------------------------------------
# FILE NAME:    CathDomainList.v4.4.0
# FILE DATE:    16.12.2024
#
# CATH VERSION: v4.4.0
# VERSION DATE: 16.12.2024
#
# FILE FORMAT:  Cath List File (CLF) Format 2.0
#
# FILE DESCRIPTION:
# Contains all classified protein domains in CATH
# for class 1 (mainly alpha), class 2 (mainly beta),
# class 3 (alpha and beta) and class 4 (few secondary structures).
#
# See 'README.file_formats' for file format information
#---------------------------------------------------------------------
1oaiA00     1    10     8    10     1     1     1     1     1    59 1.000
1go5A00     1    10     8    10     1     1     1     1     2    69 999.000
3frhA01     1    10     8    10     2     1     1     1     1    58 1.200
3friA01     1    10     8    10     2     1     1     1     2    54 1.800
```

* Note that columns 2-5 give the CATH code as per https://download.cathdb.info/cath/releases/all-releases/v4_4_0/cath-classification-data/README-cath-list-file-format.txt

* To get the AFDB and PDB CATH entries

```bash

tr -d '\000' < CATH50_h > CATH50_clean_h

cat CATH50_clean_h | grep -a '^af_' > cath_afdb.tsv
cat CATH50_clean_h | grep -av '^af_' > cath_pdb.tsv

# to get the CATH classification for AF
awk -F'_' '{
    prefix = $1"_"$2"_"$3"_"$4
    suffix = substr($0, length(prefix)+2)  # +2 to skip the underscore
    print $0 "\t" suffix
}' cath_afdb.tsv  > cath_afdb_with_class.tsv

# to get CATH classification for cath_pdb.tsv
awk '
FNR==NR && $1 !~ /^#/ {           # read cath-domain-list first
    code = $1
    domains[code] = $2 "." $3 "." $4 "." $5
    next
}
FNR!=NR {                         # read cath_pdb.tsv second
    if($1 in domains)
        print $1 "\t" domains[$1]
    else
        print $1 "\tNA"
}
' cath-domain-list-v4_4_0.txt cath_pdb.tsv > cath_pdb_with_class.tsv

```

* There are 16 entries with no info

```bash
cat cath_pdb_with_class.tsv  | grep NA
1ohfA01	NA
1ohfB01	NA
1ohfD01	NA
1xmzB00	NA
1zrtD02	NA
1zrtQ01	NA
1zrtQ02	NA
2fbwP02	NA
2h8pD00	NA
3u0kA02	NA
4ct3D00	NA
4wxfA01	NA
5einC00	NA
5excH00	NA
5excJ00	NA
5ft8B00	NA
```

* Try and get them using v4.3.0

```bash
wget --no-check-certificate https://download.cathdb.info/cath/releases/all-releases/v4_3_0/cath-classification-data/cath-domain-list-v4_3_0.txt

cat cath_pdb_with_class.tsv | grep NA | cut -f1 > cath_pdb_16_NAs.tsv

# to get CATH classification for cath_pdb.tsv
awk '
FNR==NR && $1 !~ /^#/ {           # read cath-domain-list first
    code = $1
    domains[code] = $2 "." $3 "." $4 "." $5
    next
}
FNR!=NR {                         # read cath_pdb.tsv second
    if($1 in domains)
        print $1 "\t" domains[$1]
    else
        print $1 "\tNA"
}
' cath-domain-list-v4_3_0.txt cath_pdb_16_NAs.tsv > cath_pdb_16_NAs_with_class.tsv
```

* This worked 

```bash
cat cath_pdb_with_class.tsv | grep -v NA > cath_pdb_with_class_no_NA.tsv
cat cath_pdb_with_class_no_NA.tsv cath_pdb_16_NAs_with_class.tsv > cath_pdb_with_class_final.tsv
```


* To combine

```bash
cat cath_afdb_with_class.tsv  cath_pdb_with_class_final.tsv > cath_combined_with_class.tsv
```

* To add the descriptions

```bash
wget --no-check-certificate https://download.cathdb.info/cath/releases/all-releases/v4_4_0/cath-classification-data/cath-names-v4_4_0.txt


awk '
    # Build a mapping from CATH code to name
    FNR==NR {
        if($0 ~ /^[^#]/ && NF >= 3) {
            code = $1
            name = substr($0, index($0, $3))   # take everything from 3rd column to end
            gsub(/^:/, "", name)               # remove leading colon
            cath[code] = name
        }
        next
    }
    # Process the second file
    {
        code = $2
        name = (code in cath) ? cath[code] : "UNKNOWN"
        print $0 "\t" name
    }
' cath-names-v4_4_0.txt cath_combined_with_class.tsv > cath_combined_with_names.tsv


mv cath_combined_with_names.tsv cath.tsv

# the actual download for the final db - name 'cath'
foldseek databases CATH50 CATH50 tmp --threads 4


```