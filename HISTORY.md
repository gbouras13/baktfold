# History

0.4.0 (unreleased)
------------------

* Adds `baktfold json`, which reconstitutes all non-Foldseek outputs (`.gff3`/`.gbff`/`.embl`/`.tsv`/`.inference.tsv`/`.faa`/`.ffn`/`.fna`/`.summary.txt`/`.json`) from a Baktfold or Bakta JSON — no database, ProstT5 or Foldseek required. Every JSON now carries a self-describing `baktfold_run` provenance block
* Adds an `Annotation_Confidence` column (high/medium/low) to the JSON and TSV outputs, derived from reciprocal coverage, sequence identity, ProstT5 confidence and E-value
* Adds `--report-prob`, which reports Foldseek's homology probability (`prob`) in each `*_tophit.tsv`. Off by default: it requires Foldseek to compute full alignments (`-a 1`), so searches are slower and use more disk
* Adds multi-GPU support (#33)
* Adds percentages to `.summary.txt`: every count as a share of the total CDS count, and the Baktfold hit, function and remaining-hypothetical counts additionally as a share of the `CDS beginning hypotheticals` — the set `baktfold` is actually given, so usually the more informative denominator (thanks @TlaskalV)
* Fixes three INSDC (`.gbff`/`.embl`) output bugs: `translation_table`/`other_genbank` were passed in swapped order (giving `/translation_table=False` on every CDS), the CDS genetic-code qualifier now uses the standard `transl_table`, and tRNA `/inference` was swapped between the Bakta and Prokka branches. The GFF3 output — the route used for GenBank submission via `table2asn_GFF` — was not affected
* Fixes `.summary.txt` miscounting annotated CDS as still-hypothetical when a retained secondary hit was described "hypothetical protein"
* Migrates pandas to polars (output verified byte-identical); lazy imports drop `import baktfold` from ~16 s to ~1-2 s; assorted database-lookup and TSV-write speedups
* Hardens IO: atomic writes, path-traversal-safe archive extraction, download timeouts and enforced md5 checks, narrowed bare `except:` clauses
* Adds a unit test suite and golden output regression tests covering every subcommand, plus publicly visible GitHub CI across supported Python versions
* Narrows declared Python support to >=3.10,<3.14 (the versions actually tested)
* Fixes `--foldseek-gpu` running out of GPU memory on a subset of genomes. `baktfold run` predicts 3Di and then launches Foldseek from the *same* process, and PyTorch's caching allocator kept ProstT5's weights and activations reserved against the card after the model itself was gone — so the `foldseek --gpu 1` subprocess started on a GPU that still looked full and died with `CUDA error: out of memory : .../libmarv/src/cudasw4.cuh` / `Error: Ungapped prefilter matching step died`. ProstT5's cached GPU memory is now released once 3Di prediction finishes (and once `--autotune` finishes) before Foldseek is launched. Larger genomes and longer proteins failed more often; around 4% of fungal and 9% of protist genomes were affected in testing, and the failure was silent per-genome
* Fixes a crash on genomes where every CDS already carries a product, i.e. zero hypothetical proteins. These wrote a 0-byte `<prefix>_3di.fasta` and `<prefix>_aa.fasta`; Baktfold checked only that those files existed, built an empty Foldseek query database from them, and Foldseek died on a zero-dimension kernel launch reported as `CUDA error: invalid configuration argument : .../libmarv/src/pssm.cuh` / `Error: Ungapped prefilter matching step died`. Baktfold now detects the empty query set, skips both the ProstT5 model load and the Foldseek search, and writes the normal annotation outputs for the (already fully annotated) genome
* Fixes `TypeError: len() of unsized object` aborting `baktfold predict` for an entire genome. The filter that drops empty ProstT5 predictions (#47) used `len()`, which raises rather than returning 0 on a zero-dimensional numpy array — the shape produced when a single-residue prediction is squeezed to a scalar. pholdlib records such an entry as a failure but still returns it, so one unusable CDS crashed the run and left every other protein in the genome unannotated. Unusable predictions (zero-length, zero-dimensional or None) are now dropped at both filter sites, and the count and ids are logged rather than silently discarded
* Fixes the `.inference.tsv` header being one column short on eukaryotic inputs. When an input repeats locus tags across CDS (common in Ensembl eukaryotic GenBank, and what `has_duplicate_locus` detects) each row carries an extra `ID` field, but the header *substituted* `ID` for `Length` instead of adding it — so every column from `Product` rightwards was named after the value to its left, and `CATH` had no header at all. Reading the file by column name returned the CDS length for `Product`: non-empty and not "hypothetical protein", so every CDS scored as functionally annotated. Header and row are now built by a single shared function with a unit-test invariant that they agree for every combination of `--custom-db`, `--fast` and `--structure-dir`. The annotations themselves were never affected, but **`.inference.tsv` files written by earlier versions from eukaryotic input carry the wrong header** — re-parse them with `Locus, ID, Length, Product, ...`, or regenerate them with `baktfold json`. Bacterial and `proteins`-mode outputs were never affected
* Fixes `--omit-probs` crashing 3Di output with `TypeError: 'NoneType' object is not subscriptable` in `predict`, `proteins-predict`, `run` and `proteins`. The flag was forwarded to pholdlib as `output_probs=False`, which returns no per-residue probabilities, so confidence masking had nothing to index. Per-residue probabilities are now always computed; `--omit-probs` only skips writing `*_prostT5_3di_all_probabilities.json`, and masking and the mean-probability CSV are unaffected. Long-standing rather than a v0.4.0 regression (also present in v0.3.0)
* Fixes a crash writing the feature TSV for rRNA products other than 16S/23S/5S (#40): the rRNA `db_xrefs` list carried a `None` beside a string, so every writer that sorts `db_xrefs` raised `TypeError: '<' not supported between instances of 'NoneType' and 'str'`
* Golden regression tests now tolerate ProstT5 run-to-run noise (3Di differs at a handful of residues between identical runs on ROCm), so they flag real output changes rather than inference nondeterminism
* Fixes label selection masking informative hits: the database hierarchy took the top-ranked hit unconditionally and only afterwards renamed "Uncharacterized protein ..." to "hypothetical protein", so a functional hit from a lower-precedence database was discarded, the CDS was not flagged hypothetical, and it still received an `Annotation_Confidence`. Hits whose description mentions "hypothetical" or "uncharacterised/uncharacterized" anywhere (any case), or that Bakta's product cleanup discards, are now skipped and the next database is used; the confidence tier is taken from the hit that supplied the label; an existing curated product is never replaced by an uninformative hit
* Fixes `compare` and `proteins-compare` searching *unmasked* amino acids while `run` and `proteins` search masked ones (the AA FASTA from `--predictions-dir` was not copied when the input parser had already written an unmasked one), so the same predictions gave different hits depending on the entry point
* Fixes `proteins` and `proteins-predict` with a Bakta proteins JSON predicting and searching every protein and overwriting Bakta's curated products; only Bakta's hypotheticals are annotated, as documented
* Fixes an existing output directory/file without `-f` only logging an error and continuing (exit 0, new results mixed with stale ones); it now exits with status 1. Non-GenBank input to `convert-prokka`/`convert-euk` now exits cleanly instead of crashing
* Fixes eukaryotic GFF3 output carrying no Baktfold labels: CDS lines now include `product`, `Name`, `locus_tag` and `Dbxref`, and each CDS label is copied to its parent mRNA in all outputs (GFF3/GBFF/TSV/JSON). Also fixes the mRNA `Note` being filled from db_xrefs and unescaped GFF3 mRNA attributes. Labels in earlier eukaryotic outputs can be regenerated with `baktfold json`; the CDS phase and partial-end fixes below need the input re-converted with `convert-euk` and re-annotated, because older JSONs did not record `codon_start` or partial ends
* Fixes `Annotation_Confidence` disagreeing with its documentation: thresholds documented as ≥ (80% coverage, 60% ProstT5 confidence) were applied as >, the ProstT5 "medium" band (45–60%) made the tier non-monotonic (confidence 61 scored lower than 50), and for structure input "medium" did not require E-value < 1e-5. The tiers now follow `docs/output.md` exactly and are monotonic in ProstT5 confidence
* Fixes query/target coverage (`qCov`/`tCov`) being one residue short: Foldseek coordinates are 1-based inclusive, so a full-length alignment is now reported as 1.0 rather than 0.99. This affects the `*_tophit.tsv` columns, the JSON `query_cov`/`subject_cov`, and hits at the 80% confidence-tier boundary
* Fixes partial eukaryotic CDS in GenBank/EMBL output: `/codon_start` was always written as 1 and partial ends (`<`/`>`) were dropped, so the location no longer matched `/translation`. `convert-euk` now records `codon_start` and `partial_start`/`partial_stop` on CDS, mRNA and gene features, and the writers restore them (all 8,163 CDS of the Entamoeba test genome now round-trip to the input location and `/codon_start` exactly)
* Fixes eukaryotic GFF3 CDS phase: phase was computed as bases consumed rather than bases to skip, ignored `/codon_start`, and minus-strand parts were reversed although they are already stored 5'→3' (exons were also numbered 3'→5'). Phases are now correct and partial CDS carry `partial=true`
* CATH-sourced labels are now explicitly domain-level: "<superfamily> domain-containing protein" (e.g. "SH2 domain" -> "SH2 domain-containing protein"), instead of the superfamily name verbatim, which could present a fold-level match as the whole protein's function
* PDB hits are labelled per chain when the database has the new `pdb_chain` table: the chain's UniProt name (via SIFTS), else its PDB entity name, instead of the PDB entry title ("Crystal structure of ... in complex with ..."). The PDB cross-reference now names the chain (`pdb:pdb_6iv5_A`), and when a PDB hit supplies the label from a UniProt name, `UniProtKB:<accession>` is added so the label can be traced. Unmapped chains fall through to the next database. Older databases keep entry titles (with a warning). Build scripts: `db-scripts/init-pdb-chains.py`, `db-scripts/fetch-uniprot-names.py`
* Fixes the `Database:` line in the TSV/GFF3/GBFF/EMBL headers, which printed the Baktfold *software* version; it now prints the database version (0.1.0), defined once in `src/baktfold/databases/version.py` and used by the headers, `.summary.txt`, JSON provenance and downloads (previously 0.1.0, the software version and 0.0.1 in different places)
* Fixes sORF `/inference` citing Baktfold's version as Bakta's (`ab initio prediction:Bakta:0.3`); it now uses the Bakta version recorded in the input JSON
* Robustness: duplicate FASTA IDs in protein input now stop the run (they were silently merged, both records getting the second one's label); a run where ProstT5 produced no prediction for any protein now exits non-zero instead of finishing with nothing annotated, and partial failures are counted; `--custom-db` with `--foldseek-gpu` accepts the documented padded database path; `--gpus` respects an inherited `CUDA_VISIBLE_DEVICES` (e.g. under SLURM) when passing devices to Foldseek; `proteins-compare` validates the database; `-f` refuses to delete an output directory that contains the inputs; all-numeric sequence IDs no longer crash hit parsing; a missing or empty input file reports the real error. `--foldseek-gpu` now warns that `-s` and `--ultra-sensitive` do not apply in GPU mode
* Fixes `proteins-compare` with a Bakta proteins JSON searching every protein and overwriting Bakta's curated products (as `proteins` did before); only Bakta's hypotheticals are annotated, in both commands
* Fixes eukaryotic mRNA labels being swapped between alternative transcripts of one gene: a CDS is now paired with an mRNA only if it lies inside it and their splice junctions agree, and an mRNA that still matches CDS with different products keeps its own product rather than a possibly wrong one
* A run where ProstT5 returned predictions but none were usable (e.g. a single one-residue protein) no longer aborts; only a run where ProstT5 returned nothing at all exits with an error
* With `--all-proteins`, an input product renamed to "hypothetical protein" by Bakta's product cleanup is now also flagged hypothetical
* Fixes `.summary.txt` with `--all-proteins`: "beginning hypotheticals" is now the number of hypothetical CDS in the input (recorded before annotation and kept in the JSON provenance block, so `baktfold json` reproduces it), and "Baktfold function" counts only CDS whose label Baktfold supplied, not curated CDS that merely had a hit
* Fixes the eukaryotic GFF3 hierarchy: transcripts now use the input's own transcript identifier where it has one (`/transcript_id`, `/standard_name`, "transcript X" notes, or a funannotate `-T<n>` protein ID) and otherwise `<locus>-T<n>`, unique in the file (every transcript used to be `<locus>-T1`); all segments of a spliced CDS share one ID (each segment had its own, so GFF3 readers and `table2asn` saw one CDS per exon); and each CDS's `Parent` is its own transcript. CDS are paired with mRNAs by transcript identifier, then by splice structure. Also fixes eukaryotic tRNA lines, which had no ID or Parent (their exons pointed at a missing parent), and a typo in tRNA exon IDs
* `convert-euk`/`run --euk`/`json` now handle eukaryotic GenBank from more sources. Tested on NCBI RefSeq (yeast, Drosophila, Plasmodium), Ensembl (C. elegans, Ostreococcus, Ensembl Fungi/Protists), INSDC submitter annotations (Pneumocystis, Malassezia, Encephalitozoon, Entamoeba) and funannotate. Seven of the eight new sources previously crashed when writing outputs:
  * rRNA/ncRNA features were converted without their sequence (`.ffn` writer crashed)
  * features without a product or `/locus_tag` crashed the GFF3/GenBank/TSV writers
  * features without `/locus_tag` were given a random locus, so Ensembl mRNA/CDS (linked by `/gene`) lost their gene; `/gene` is now used
  * CDS without `/product` now carry "hypothetical protein"
  * rRNA/ncRNA GFF3 lines lacked the score column; eukaryotic rRNA/ncRNA no longer reuse their gene's ID
  * gene IDs are unique even when an input gives two genes the same `/gene`
  * transcript identifiers are also read from `transcript_id=` notes and from identifiers listed in an mRNA's `/note`
  All 11 test sources produce GFF3 that passes GenomeTools `gff3validator`, and every CDS translates from the GFF3 to its protein
* The version is now defined once, in `src/baktfold/utils/VERSION`, and read by packaging, `--version` and output headers alike

0.3.0 (2026-05-05)
------------------

* Adds improved [documentation](https://gbouras13.github.io/baktfold/) with tutorial
* Adds `.summary.txt` output file with summary statisics of Baktfold's performance (thanks @mbhall88 for the suggestion)
* e.g.

```
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

0.2.0 (2026-04-08)
------------------

* Add support for `bakta_proteins` JSON output format as `baktfold proteins` input - `baktfold proteins` will automatically detect the format (`.faa` or Bakta JSON) #30
* Fix tool description in JSON (inherits Bakta's and add Baktfold's) #27
* Adds comprehensive CI tests synced to Setonix via @ignatiusm [HPCI](https://github.com/Garvan-Data-Science-Platform/hpci), getting around GitHub CI runner difficulties (and adding GPU access)

0.1.1 (2026-03-30)
------------------

* Support Uniprot/Trembl etc protein headers with baktfold proteins that Foldseek trims by default
* Fix bug where embeddings were cast to half precision with `--cpu` #28


0.1.0 (2026-02-19)
------------------

* Changes to database hosting - now on HuggingFace so should be much much faster by default to download
* Introduction of `baktfold convert-prokka` for automatically converting a Prokka GenBank to the required Bakta JSON format
* Experimental support for eukaryotes
    * Via `baktfold convert-euk` for automatically converting a eukaryotic GenBank to the required Bakta JSON format
    * Then via the `--euk` flag with `baktfold run`
    * All eukaryote CDS features should be annotated - the glue-code and standards compliance support for other features is not guaranteed
    * Please try and reach out if you have any feedback
* Otherwise, we recommend [genbank_to](https://github.com/linsalrob/genbank_to) for non-Bakta GenBank conversion to Bakta JSON format


0.0.2 (2025-11-07)
------------------
* Keeps all non-overlapping top hits for CATH, not just the tophit (as multi-domain proteins can and should have multiple different hits to CATH). This is equivalent to using `--greedy-best-hits` with Foldseek

0.0.2 (2025-11-07)
------------------
* Baktfold is currently under active development. We would welcome any and all feedback (especially bugs) via Issues
* Fixes https://github.com/gbouras13/baktfold/issues/3 allowing for bioconda integration

0.0.1 (2025-11-07)
------------------

* Initial release of baktfold
* Baktfold is currently under active development. We would welcome any and all feedback (especially bugs) via Issues