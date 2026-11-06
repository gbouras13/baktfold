# Eukaryotic input validation (manual, not run in CI)

Checks that `convert-euk` + `baktfold json` turn real eukaryotic GenBank/EMBL
flatfiles from many annotation sources into valid GFF3 and GenBank/EMBL that
still match the input. It needs no database, GPU or Foldseek: the checks are on
the output writers, which is where eukaryotic input has broken before. It is
not part of CI because 8 of the 11 inputs are downloaded (~25 MB compressed)
from NCBI and Ensembl.

Run it after changing `convert-euk`, the GFF3/GenBank/EMBL writers, or the
eukaryotic transcript handling (`io/euk_transcripts.py`):

```bash
bash tests/euk_validation/fetch_sources.sh /tmp/euk_sources
bash tests/euk_validation/run_all.sh /tmp/euk_sources /tmp/euk_val | tee /tmp/euk_val/summary.txt
```

`baktfold` and Biopython must be on `PATH`. If GenomeTools (`gt`) and `gffread`
are on `PATH` too, they run as independent checks (`conda install -c bioconda
genometools-genometools gffread`).

## Sources

| Name | Annotation style | Input |
|---|---|---|
| yeast_refseq | NCBI RefSeq, incl. mitochondrion (genetic code 3) | GCF_000146045.2 |
| dmel_chr4_refseq | NCBI RefSeq, many isoforms | NC_004353.4 |
| pfal_chr1_refseq | NCBI RefSeq | NC_004325.2 |
| celegans_chrI_ensembl | Ensembl (links by `/gene`), many isoforms | WBcel235 chromosome I |
| ostreococcus_chr1_ensemblplants | Ensembl Plants | ASM9206v1 chromosome 1 |
| pneumocystis_genbank | INSDC, NIH 2021 | GCA_017301755.1 |
| malassezia_genbank | INSDC, Duke 2023 | GCA_029542905.1 |
| ecuniculi_genbank | INSDC, Genoscope; no mRNA features | GCA_000091225.2 |
| aspergillus_funannotate_ncbi | funannotate, NCBI-submitted | `tests/test_data/euks/aspergillus_subset.gbk` (not committed; skipped if absent) |
| entamoeba_insdc | INSDC, EVM-style transcript notes | `tests/test_data/protist.gbk.gz` |
| mucor_funannotate | funannotate | `tests/test_data/funannotate.gbk.gz` |

## What is checked

`check_gff3.py` (GFF3 alone):
- IDs are unique, except that all segments of one spliced CDS share one ID;
- every Parent exists; a CDS/exon hangs off a transcript or gene, an mRNA off a gene;
- each CDS, translated from its GFF3 segments and phase with its own genetic
  code, equals its protein in the JSON.

`check_gbff.py` (GenBank and EMBL against the original input):
- both parse; per-type feature counts equal the input's;
- the CDS (location, `/codon_start`) multiset equals the input's;
- per CDS, `/transl_table` and the `/transl_except`, `/ribosomal_slippage`,
  `/pseudo`, `/exception` and `/trans_splicing` qualifiers equal the input's;
- each CDS translated from its written location equals its `/translation`;
- no non-integer `/transl_table` and no `/pseudo=""`; every CDS has
  `/locus_tag` and `/product`.

Also, when installed: GenomeTools `gff3validator` must report OK, and `gffread`
must translate the GFF3.

## Expected result

Every line reports no problems, with these known exceptions:
- `dmel_chr4_refseq`: one CDS (Dmel_CG11186) differs from a plain translation at
  its `/transl_except` stop read-through, which is written to the output.
- `yeast_refseq`: 6 CDS have no `/translation`. They are `/pseudo` blocked
  reading frames with no protein in the input, so none is written.
- gffread (if compared with the JSON proteins) mistranslates yeast's 19
  mitochondrial CDS, because GFF3 has no genetic-code column.

The first run of this validation (November 2026) found and fixed, among others:
crashes at output writing on 7 of the 8 downloaded sources; lost gene/mRNA/CDS
links for Ensembl input; `/transl_table=Unknown` on every eukaryotic CDS;
invalid `/pseudo=""`; dropped `/transl_except` and `/ribosomal_slippage`; and
`/inference` claiming Prodigal for annotations Baktfold did not predict.
