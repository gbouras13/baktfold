# #!/bin/bash

# set -e

# mkdir db
# cd db

# printf "Create Baktfold database\n"

# # download NCBI Taxonomy DB
# printf "\n6/19: download NCBI Taxonomy DB ...\n"
# mkdir taxonomy
# cd taxonomy
# wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
# tar -I pigz -xf taxdump.tar.gz
# cd ..
# mv taxonomy/nodes.dmp .
# rm -rf taxonomy


############################################################################
# Setup SQLite Bakta db
############################################################################
printf "\n7/19: setup SQLite Baktfold db ...\n"
python3 ${BAKTFOLD_DB_SCRIPTS}/init-db.py --db baktfold.db

python3 ${BAKTFOLD_DB_SCRIPTS}/init-swissprot.py --taxonomy uniprot_2025_03/nodes.dmp --xml uniprot_2025_03/uniprot_sprot_filtered.xml --fasta swissprot_check.fasta --db baktfold.db 
python3 ${BAKTFOLD_DB_SCRIPTS}/init-pdb.py --tsv pdb.tsv --db baktfold.db 
# chain-level PDB names (UniProt via SIFTS, else PDB entity name); see init-pdb-chains.py
# wget https://files.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt
# wget https://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/tsv/pdb_chain_uniprot.tsv.gz
python3 ${BAKTFOLD_DB_SCRIPTS}/init-pdb-chains.py --seqres pdb_seqres.txt --sifts pdb_chain_uniprot.tsv.gz --db baktfold.db --lookup pdb_gpu.lookup --missing-out missing_accessions.txt --dry-run
python3 ${BAKTFOLD_DB_SCRIPTS}/fetch-uniprot-names.py --accessions missing_accessions.txt --out uniprot_names.tsv
python3 ${BAKTFOLD_DB_SCRIPTS}/init-pdb-chains.py --seqres pdb_seqres.txt --sifts pdb_chain_uniprot.tsv.gz --db baktfold.db --lookup pdb_gpu.lookup --uniprot-names uniprot_names.tsv
python3 ${BAKTFOLD_DB_SCRIPTS}/init-cath.py --tsv cath.tsv --db baktfold.db 
python3 ${BAKTFOLD_DB_SCRIPTS}/init-afdbclusters.py --taxonomy uniprot_2025_03/nodes.dmp --xml uniprot_2025_03/uniprot_trembl_AFDBClusters.xml --missingxml uniprot_2025_03/missing_AFDBCluster_sprot.xml  --db baktfold.db 





