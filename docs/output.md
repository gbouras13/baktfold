# `baktfold` Output Files

## Main Outputs

* `.faa` which holds all amino acid sequences of CDSs (not only hypotheticals) (like Bakta)
* `.ffn` which holds all nucletotide sequences of CDSs (like Bakta)
* `.fna` which holds the input genome FASTA (parsed from the JSON) (like Bakta)
* `_3di.fasta` which holds all Foldseek 3Di sequences of predicted CDSs as predicted by ProstT5
* `.gbff` `.gff3` and `.embl` are GenBank, GFF3 and EMBL format output files (like Bakta) with full genome annotations including the input annotations combined with Baktfold's additional annotations 
* `{prefix}.tsv` is the same as Bakta's `.tsv output`, containing all full genome annotations including the input annotations combined with Baktfold's additional annotations 
* `.json` which is a Bakta JSON format output file with all genome annotations including the input annotations combined with Baktfold's additional annotations 
* `*_tophit.tsv` which have detailed alignment statistics for the Foldseek output tophits for each of the four consituent Baktfold databases (Swiss-Prot, AFDB clusters, CATH and PDB) 
    * For example:

```bash
query	target	bitscore	fident	evalue	qStart	qEnd	qLen	qCov	tStart	tEnd	tLen	tCov
MEGJMN_070	AF-A0A1I3V7E0-F1-model_v6	292	0.41	2.619e-06	1	91	93	0.97	1	95	99	0.95
```

* There will also be a `custom_database_tophit.tsv` if you have used `--custom-db`

* `baktfold.inference.tsv` which contains the new annotated funciton (`Function` column) and all annotated tophits (if they exist) for all proteins that were attempted to be annotated by Baktfold
* Note: This will contain an extra column `CustomDB` if you have used `--custom-db` with a custom database 
    * For example:

```bash
ID	Length	Product	Swissprot	AFDBClusters	PDB	CATH
MEGJMNBEGN_27	162	HTH-type quorum-sensing regulator RhlR	swissprot_P54292	afdbclusters_A0A9E1VSB0	pdb_5l09	cath_3sztB01
MEGJMNBEGN_30	68	hypothetical protein				
MEGJMNBEGN_70	94	hypothetical protein		afdbclusters_A0A1I3V7E0		
```

* `baktfold.summary.txt` which contains basic summary statistics showing how many CDS baktfold was able to annotate
* `CDS beginning hypotheticals` is the amount of hypotheticals parsed from the input JSON
* `CDS annotated with Baktfold database hit` is the number of CDS that have at least 1 hit to a constituent Baktfold DB
* `CDS annotated with Baktfold function` is the number of CDS whose Baktfold tophit has a non-hypothetical function transferred
* `CDS remaining hypotheticals` is the number of CDS remaining hypothetical most `baktfold`
    * For example

```bash
Annotation:
CDS count: 2635
CDS beginning hypotheticals: 55
CDS annotated with Baktfold database hit: 12
CDS annotated with Baktfold function: 7
CDS remaining hypotheticals: 48

Baktfold:
Software: v0.2.0
Database: v0.1.0
DOI: https://doi.org/10.64898/2026.03.31.715528
URL: github.com/gbouras13/baktfold
``` 

## Supplementary Outputs

*  `_prostT5_3di_mean_probabilities.csv` - contains the mean ProstT5 probability score for each CDS. These are equivalent to the probability of how similar the overall ProstT5 3Di sequence is predicted to be compared to its Alphafold2 baseline
*  `_prostT5_3di_all_probabilities.json` - contains the ProstT5 probabilities for each residue for each CDS, in the json format

## Optional Outputs

*  `_embeddings_per_protein.h5` - contains the ProstT5 embeddings for each protein in the h5 format
*  `_embeddings_per_residue.h5` - contains the ProstT5 embeddings for each residue in the h5 format

