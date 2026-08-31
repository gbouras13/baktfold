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

* `baktfold.inference.tsv` which contains the new annotated funciton (`Product` column), an `Annotation_Confidence` column and all annotated tophits (if they exist) for all proteins that were attempted to be annotated by Baktfold
* Note: This will contain an extra column `Custom_DB` if you have used `--custom-db` with a custom database, and extra `TMscore` and `LDDT` columns if you ran with user-provided structures (`--structure-dir`)
    * For example:

```bash
ID	Length	Product	Annotation_Confidence	Swissprot	AFDBClusters	PDB	CATH
MEGJMNBEGN_27	162	HTH-type quorum-sensing regulator RhlR	high	swissprot_P54292	afdbclusters_A0A9E1VSB0	pdb_5l09	cath_3sztB01
MEGJMNBEGN_30	68	hypothetical protein					
MEGJMNBEGN_70	94	hypothetical protein	low		afdbclusters_A0A1I3V7E0		
```

### Annotation confidence

Each Baktfold annotation is labelled `high`, `medium` or `low` confidence (in the `Annotation_Confidence` column of the TSV and the `annotation_confidence` field of the JSON — it is intentionally *not* written to the GFF3/GenBank/EMBL outputs) to make the structural-alignment results more interpretable. The heuristics mirror [Phold](https://github.com/gbouras13/phold):

* **High**: both query and target have ≥80% reciprocal alignment coverage, *and* one of (i) >30% amino acid identity (the "light zone" of sequence homology), (ii) mean ProstT5 confidence ≥60%, or (iii) alignment E-value < 1e-10.
* **Medium**: either the query or target has ≥80% coverage, *and* one of (i) >30% amino acid identity or (ii) ProstT5 confidence between 45–60%, *and* E-value < 1e-5.
* **Low**: all other hits below the E-value threshold that do not meet the above (low coverage, low identity, low ProstT5 confidence, or an E-value near the 0.001 threshold).

If Baktfold is run with user-provided structures (`--structure-dir`) instead of ProstT5, the ProstT5-confidence criteria are dropped and the output additionally includes Foldseek `TMscore` (template modeling score) and `LDDT` (local distance difference test) values to help assess quality.

A `low` annotation is **not** necessarily wrong — only that it should be trusted less than a `high` (extremely trustworthy) or `medium` (very trustworthy) annotation, and likely warrants a bit more manual curation.

* `baktfold.summary.txt` which contains basic summary statistics showing how many CDS baktfold was able to annotate
* `CDS beginning hypotheticals` is the amount of hypotheticals parsed from the input JSON
* `CDS annotated with Baktfold database hit` is the number of CDS that have at least 1 hit to a constituent Baktfold DB
* `CDS annotated with Baktfold function` is the number of CDS whose Baktfold tophit has a non-hypothetical function transferred
* `CDS remaining hypotheticals` is the number of CDS remaining hypothetical post `baktfold`
* Every count is also given as a percentage of the total CDS count. The three counts that describe what `baktfold` did are additionally given as a percentage of the `CDS beginning hypotheticals`, as by default these are the only CDS `baktfold` is given (unless you use `--all-proteins`) — so this is usually the more informative denominator
    * For example

```bash
Annotation:
CDS count: 2635
CDS beginning hypotheticals: 55 (2.1% of CDS)
CDS annotated with Baktfold database hit: 12 (0.5% of CDS; 21.8% of beginning hypotheticals)
CDS annotated with Baktfold function: 7 (0.3% of CDS; 12.7% of beginning hypotheticals)
CDS remaining hypotheticals: 48 (1.8% of CDS; 87.3% of beginning hypotheticals)

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

## Reconstituting outputs from the JSON

The `.json` is the authoritative record of a `baktfold` run. You can regenerate every standard output above from it (without re-running ProstT5 or Foldseek) with [`baktfold json`](run.md#baktfold-json):

* **Reconstitutable** from the `.json`: `.gff3`, `.gbff`, `.embl`, `.tsv`, `.inference.tsv`, `.faa`, `.ffn`, `.fna`, `.summary.txt` (and a fresh `.json`).
* **Not** reconstitutable (only produced by an actual `compare`/`run`): the Foldseek `*_tophit.tsv`, `foldseek_results_*.tsv`, `_3di.fasta`, the ProstT5 probability files and the embedding `.h5` files.

