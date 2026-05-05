# History

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