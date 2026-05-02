# Running `baktfold`

`baktfold` has been split into two overall modules regardless of your input - `predict` and `compare`. 

The reason `baktfold` has been split as such was to enable optimal resource usage between GPU and CPU in cluster environments (before Foldseek-GPU).

However, now if you have an NVIDIA GPU available, please use `baktfold run`, which combines both in one as both `predict` and `compare` can utilise GPU acceleration.

If you do not have an NVIDIA GPUI available, we recommend using `baktfold predict` with GPU and `baktfold compare` with CPU.

If you do not have an available GPU, please use `baktfold run` with `--cpu`.

If you are having troubing runing `baktfold` offline after installation (e.g. on a HPC), you may need to add `TRANSFORMERS_OFFLINE=True` to your environment. 

All off the above applies for `baktfold proteins` with `baktfold proteins-predict` and `baktfold proteins-compare` as the analogous command

## All commands

```bash
Usage: baktfold [OPTIONS] COMMAND [ARGS]...

  Main command line interface for baktfold.

Options:
  -h, --help     Show this message and exit.
  -V, --version  Show the version and exit.

Commands:
  autotune          Determines optimal batch size for 3Di prediction with...
  citation          Print the citation(s) for this tool
  compare           Runs Foldseek vs baktfold db
  convert-euk       (Experimental) Converts eukaryotic GenBank to Bakta...
  convert-prokka    Converts Prokka GenBank to Bakta format json
  createdb          Creates foldseek DB from AA FASTA and 3Di FASTA input...
  install           Installs ProstT5 model and baktfold database
  predict           Uses ProstT5 to predict 3Di tokens - GPU recommended
  proteins          baktfold proteins-predict then comapare all in one -...
  proteins-compare  Runs Foldseek vs baktfold db on proteins input
  proteins-predict  Runs ProstT5 on a multi FASTA input - GPU recommended
  run               baktfold predict then comapare all in one - GPU...
```

## Input

`baktfold run` takes a [Bakta](https://github.com/oschwengers/bakta) output JSON format genome annotation output. `baktfold` supports conversion from Prokka GenBank files using `baktfold convert-prokka` and from eukaryotic GenBank files using `baktfold convert-euk`. 

**For genomes annotated using other methods, we recommend [genbank_to](https://github.com/linsalrob/genbank_to) to convert you GenBank file to the Bakta JSON format.**

`baktfold proteins` takes a `.faa` FASTA file of amino acid protein sequences as input. You can also use a `bakta_proteins` output JSON file as input if you have already annotated your protein FASTA with `bakta`.

## Subcommands

### `baktfold autotune`

If you are running baktfold on a large dataset, you might want to determine the optimal batchsize for ProstT5 for your data. 

```bash
Usage: baktfold autotune [OPTIONS]

  Determines optimal batch size for 3Di prediction with your hardware

Options:
  -h, --help             Show this message and exit.
  -V, --version          Show the version and exit.
  -i, --input PATH       Optional path to input file of proteins if you do not
                         want to use the default sample of  5,000 swissprot proteins
                         proteins
  --cpu                  Use CPU only
  -t, --threads INTEGER  Number of threads  [default: 1]
  -d, --database TEXT    Path to Baktfold's database
  --min-batch INTEGER    Minimum batch size  [default: 1]
  --step INTEGER         Batch size step increment  [default: 10]
  --max-batch INTEGER    Maximum batch size  [default: 251]
  --sample-seqs INTEGER  Subsample size of input proteins  [default: 500]
```

Example usage (assuming you have run `baktfold install`)

```bash
baktfold autotune -d baktfold_db
```

### `baktfold predict`

`predict` uses the [ProstT5](https://github.com/mheinzinger/ProstT5) protein language model to translate protein amino acid sequences to the 3Di token alphabet used by [foldseek](https://github.com/steineggerlab/foldseek). This module is greatly accelerated if you have a GPU available and is recommended.

```bash
Usage: baktfold predict [OPTIONS]

  Uses ProstT5 to predict 3Di tokens - GPU recommended

Options:
  -h, --help                     Show this message and exit.
  -V, --version                  Show the version and exit.
  -i, --input PATH               Path to input file (Bakta JSON
                                 format)  [required]
  -o, --output PATH              Output directory   [default: output_baktfold]
  -t, --threads INTEGER          Number of threads  [default: 1]
  -p, --prefix TEXT              Output files' prefix  [default: baktfold]
  -d, --database TEXT            Path to Baktfold's database
  -f, --force                    Force overwrites output directory
  --autotune                     Run autotuning to detect and set best batch
                                 size for local hardware (recommended only for
                                 large dataset, e.g. thousands of proteins)
  --batch-size INTEGER           Batch size for ProstT5 (1 is usually fastest)
                                 [default: 1]
  --cpu                          Use CPU only.
  --omit-probs                   Do not output per residue 3Di probabilities
                                 from ProstT5
  --save-per-residue-embeddings  Save ProstT5 embeddings per resuide in a H5
                                 file
  --save-per-protein-embeddings  Save ProstT5 embeddings as means per protein
                                 in a H5 file
  --mask-threshold FLOAT         Mask 3Di residues below this value of ProstT5
                                 confidence for Foldseek searches  [default:
                                 25]
  -a, --all-proteins             Annotate all proteins (not only
                                 hypotheticals)
```

* Example usage (assuming you have run `baktfold install`) using a batch size of `10` using `-b`

```bash
baktfold predict -i bakta.json -o baktfold_predict_output -b 10
```

* Use `-a` to annotate all proteins (not just hypotheticals)


### `baktfold compare`

`baktfold compare` runs [foldseek](https://github.com/steineggerlab/foldseek) to compare ProstT5 predictions generated by `baktfold predict` to the baktfold databases.

Alternatively, if you have provided pre-generated `.pdb` or `,cif` format protein structures for you proteins, you can specify those by specifiying `--structures --structure-dir <directory>`. 

`baktfold compare` does not use a GPU by default. However, if you have an NVIDIA GPU available, you can utilise Foldseek-GPU acceleration using `--foldseek-gpu`. Note that you need to make sure your also run `baktfold install` with `--foldseek-gpu` prior. Regardless of whether you use `--foldseek-gpu` or not, it is recommended to use as many CPU threads with `-t` as you can (as the GPU only accelerates Foldseek's prefilter, not the alignment step).

Example usage of `baktfold compare` following `baktfold predict`

```bash
baktfold compare -i bakta.json-o baktfold_compare_output  --predictions_dir baktfold_predict_output
```

* Example usage if you have `.pdb` or `.cif` format structures available for your input proteins (note: the .pdb file names must be called cds_id.pdb, where cds_id is the CDS ids). You can see an example in the `tests/test_data/pdbs` directory [here](https://github.com/gbouras13/baktfold/tree/main/tests/test_data/pdbs).

```bash
baktfold compare -i bakta.json -o baktfold_compare_output_custom  --structure-dir directory_with_pdbs  -t 8
```

* If you have a custom database of protein structures would would additionally like to search against, Baktfold will support this using `--custom-db`, with the 
You will first need to use `foldseek createdb` first to create a Foldseek compatible database. For example, assuming you have protein structures in a directory called `custom-structures`

```bash
mkdir -p custom_db
foldseek createdb custom_structures/ custom_db/db
# if you want to use with `--foldseek-gpu`
foldseek makepaddedseqdb custom_db/db custom_db/db_pad

baktfold compare -i bakta.json -o baktfold_compare_output_custom  --custom_db custom_db/db -t 8
# if you want to use with `--foldseek-gpu`
baktfold compare -i bakta.json -o baktfold_compare_output_custom  --custom_db custom_db/db_pad -t 8 --foldseek-gpu
```

* If you have custom annotations for your custom database, these can be specified using `--custom-annotations`. This file must be a headerless 2 column TSV file where column 1 are your DB ids, and column 2 is the annotation. 

e.g. 

```bash
MEGJMN_003	dummy annotation
MEGJMN_005	dummy annotation
MEGJMN_010	dummy annotation
MEGJMN_024	dummy annotation
MEGJMN_027	dummy annotation
MEGJMN_030	dummy annotation
MEGJMN_042	dummy annotation
```

For an example, see this (link)[https://github.com/gbouras13/baktfold/tree/main/tests/test_data/custom_db] for an example custom database, with `dummy_custom_db_annotations.tsv` representing the custom annotations file. 

* If your input genome is eukaryotic, please use `--euk` to ensure all common eukaryotic genomic features are parsed correctly.

* If you want to keep all Foldseek hits (not just the tophits), use `--keep-tmp-files`

* If you want to skip searching against the AFDB (the most expensive step), use `--fast`

* If your input genome was not annotated with Bakta or Prokka and you want to specify various specific inference tools in your output, please use `--cds-program`, `--trna-program`, `--tmrna-program`, `--ncrna-program`, `--rrna-program`

* If you would like to specify custom Foldseek parameters, use `--extra-foldseek-params` along with double quotes enclosing your extra parameters e.g. `--extra-foldseek-params "-c 0.8"`

```bash
Usage: baktfold compare [OPTIONS]

  Runs Foldseek vs baktfold db

Options:
  -h, --help                    Show this message and exit.
  -V, --version                 Show the version and exit.
  -i, --input PATH              Path to input file (Bakta JSON
                                format)  [required]
  --predictions-dir PATH        Path to output directory from Baktfold predict
  --structure-dir PATH          Path to directory with .pdb or .cif file
                                structures (IDs need to be in file names, i.e
                                id.pdb or id.cif)
  -o, --output PATH             Output directory   [default: output_baktfold]
  -t, --threads INTEGER         Number of threads  [default: 1]
  -p, --prefix TEXT             Output files' prefix  [default: baktfold]
  -d, --database TEXT           Path to Baktfold's database
  -f, --force                   Force overwrites output directory
  -e, --evalue FLOAT            Evalue threshold for Foldseek  [default: 1e-3]
  -s, --sensitivity FLOAT       Sensitivity parameter for Foldseek  [default:
                                9.5]
  --keep-tmp-files              Keep temporary intermediate files (e.g. large
                                foldseek_results.tsv of all Foldseek hits)
  --max-seqs INTEGER            Maximum results per query sequence allowed to
                                pass the prefilter (saves disk space for
                                enormous datasets)  [default: 1000]
  --ultra-sensitive             Run with maximum sensitivity by skipping
                                Foldseek prefilter (not recommended for large
                                datasets)
  --extra-foldseek-params TEXT  Extra Foldseek search params
  --custom-db TEXT              Path to custom database
  --foldseek-gpu                Enable Foldseek-GPU search acceleration
  --custom-annotations PATH     Custom Foldseek DB annotations (2 column tsv:
                                Foldseek headers, description)
  --euk                         Eukaryotic input genome.
  --fast                        Skips Foldseek search against AFDB Clusters.
  --cds-program TEXT            CDS prediction tool string for compliant
                                outputs if using non-Bakta or Prokka input
  --trna-program TEXT           tRNA prediction tool string for compliant
                                outputs if using non-Bakta or Prokka input
  --tmrna-program TEXT          tmRNA profile prediction tool string for
                                compliant outputs if using non-Bakta or Prokka
                                input
  --rrna-program TEXT           rRNA profile prediction tool string for
                                compliant outputs if using non-Bakta or Prokka
                                input
  --ncrna-program TEXT          ncRNA profile prediction tool string for
                                compliant outputs if using non-Bakta or Prokka
                                input
  -a, --all-proteins            Annotate all proteins (not only hypotheticals)
```


### `baktfold run`

* `baktfold run` runs `baktfold predict` and `baktfold compare` together in one command. Recommended if you are running `baktfold` with an NVIDIA GPU, along with `--foldseek-gpu`.
* Also recommended if you are running `baktfold` in a lower-resource environment without a GPU - you will also need to specify `--cpu`.
* See `baktfold predict` and `baktfold compare` for more details on specific parameters

Example usage where NVIDIA GPU is available:

```bash
baktfold proteins -i proteins,faa -o baktfold_proteins_output -d baktfold_db -t 8 --foldseek-gpu
```

```bash
Usage: baktfold run [OPTIONS]

  baktfold predict then comapare all in one - GPU recommended

Options:
  -h, --help                     Show this message and exit.
  -V, --version                  Show the version and exit.
  -i, --input PATH               Path to input file in Bakta JSON format  [required]
  -o, --output PATH              Output directory   [default: output_baktfold]
  -t, --threads INTEGER          Number of threads  [default: 1]
  -p, --prefix TEXT              Output files' prefix  [default: baktfold]
  -d, --database TEXT            Path to Baktfold's database
  -f, --force                    Force overwrites output directory
  --autotune                     Run autotuning to detect and set best batch
                                 size for local hardware (recommended only for
                                 large dataset, e.g. thousands of proteins)
  --batch-size INTEGER           Batch size for ProstT5 (1 is usually fastest)
                                 [default: 1]
  --cpu                          Use CPU only.
  --omit-probs                   Do not output per residue 3Di probabilities
                                 from ProstT5
  --save-per-residue-embeddings  Save ProstT5 embeddings per resuide in a H5
                                 file
  --save-per-protein-embeddings  Save ProstT5 embeddings as means per protein
                                 in a H5 file
  --mask-threshold FLOAT         Mask 3Di residues below this value of ProstT5
                                 confidence for Foldseek searches  [default:
                                 25]
  -e, --evalue FLOAT             Evalue threshold for Foldseek  [default:
                                 1e-3]
  -s, --sensitivity FLOAT        Sensitivity parameter for Foldseek  [default:
                                 9.5]
  --keep-tmp-files               Keep temporary intermediate files (e.g. large
                                 foldseek_results.tsv of all Foldseek hits)
  --max-seqs INTEGER             Maximum results per query sequence allowed to
                                 pass the prefilter (saves disk space for
                                 enormous datasets)  [default: 1000]
  --ultra-sensitive              Run with maximum sensitivity by skipping
                                 Foldseek prefilter (not recommended for large
                                 datasets)
  --extra-foldseek-params TEXT   Extra Foldseek search params
  --custom-db TEXT               Path to custom database
  --foldseek-gpu                 Enable Foldseek-GPU search acceleration
  --custom-annotations PATH      Custom Foldseek DB annotations (2 column tsv:
                                 Foldseek headers, description)
  --euk                          Eukaryotic input genome.
  --fast                         Skips Foldseek search against AFDB Clusters.
  --cds-program TEXT             CDS prediction tool string for compliant
                                 outputs if using non-Bakta or Prokka input
  --trna-program TEXT            tRNA prediction tool string for compliant
                                 outputs if using non-Bakta or Prokka input
  --tmrna-program TEXT           tmRNA profile prediction tool string for
                                 compliant outputs if using non-Bakta or
                                 Prokka input
  --rrna-program TEXT            rRNA profile prediction tool string for
                                 compliant outputs if using non-Bakta or
                                 Prokka input
  --ncrna-program TEXT           ncRNA profile prediction tool string for
                                 compliant outputs if using non-Bakta or
                                 Prokka input
  -a, --all-proteins             Annotate all proteins (not only
                                 hypotheticals)
```

### `baktfold proteins-predict`

* Identical to `baktfold predict`, but instead takes a FASTA input file of amino acid protein sequences. Useful for bulk annotation of protein catalogs.

* Note: You can also use a `bakta_proteins` output JSON file as input if you have already annotated your protein FASTA with `bakta`

* Example usage 

```bash
baktfold proteins-predict -i proteins.faa -o baktfold_proteins_predict_output -d baktfold_db -b 10
```

```bash
Usage: baktfold proteins-predict [OPTIONS]

  Runs ProstT5 on a multi FASTA input - GPU recommended

Options:
  -h, --help                     Show this message and exit.
  -V, --version                  Show the version and exit.
  -i, --input PATH               Path to input file (amino acid FASTA or Bakta 
                                  proteins output JSON format) [required]
  -o, --output PATH              Output directory   [default: output_baktfold]
  -t, --threads INTEGER          Number of threads  [default: 1]
  -p, --prefix TEXT              Output files' prefix  [default: baktfold]
  -d, --database TEXT            Path to Baktfold's database
  -f, --force                    Force overwrites output directory
  --autotune                     Run autotuning to detect and set best batch
                                 size for local hardware (recommended only for
                                 large dataset, e.g. thousands of proteins)
  --batch-size INTEGER           Batch size for ProstT5 (1 is usually fastest)
                                 [default: 1]
  --cpu                          Use CPU only.
  --omit-probs                   Do not output per residue 3Di probabilities
                                 from ProstT5
  --save-per-residue-embeddings  Save ProstT5 embeddings per resuide in a H5
                                 file
  --save-per-protein-embeddings  Save ProstT5 embeddings as means per protein
                                 in a H5 file
  --mask-threshold FLOAT         Mask 3Di residues below this value of ProstT5
                                 confidence for Foldseek searches  [default:
                                 25]
```

### `baktfold proteins-compare`

* Identical to `baktfold compare`, but instead takes a FASTA input file of amino acid protein sequences. Useful for bulk annotation of protein catalogs.

* Note: You can also use a `bakta_proteins` output JSON file as input if you have already annotated your protein FASTA with `bakta`


Example usage 

```bash
baktfold proteins-compare -i proteins.faa --predictions_dir baktfold_proteins_predict_output -o baktfold_proteins_compare_output  -t 8
```

```bash
Usage: baktfold proteins-compare [OPTIONS]

  Runs Foldseek vs baktfold db on proteins input

Options:
  -h, --help                    Show this message and exit.
  -V, --version                 Show the version and exit.
  -i, --input PATH               Path to input file (amino acid FASTA or Bakta 
                                  proteins output JSON format) [required]
  --predictions-dir PATH        Path to output directory from Baktfold
                                proteins-predict
  --structure-dir PATH          Path to directory with .pdb or .cif file
                                structures. The CDS IDs need to be in the name
                                of the file
  -o, --output PATH             Output directory   [default: output_baktfold]
  -t, --threads INTEGER         Number of threads  [default: 1]
  -p, --prefix TEXT             Output files' prefix  [default: baktfold]
  -d, --database TEXT           Path to Baktfold's database
  -f, --force                   Force overwrites output directory
  -e, --evalue FLOAT            Evalue threshold for Foldseek  [default: 1e-3]
  -s, --sensitivity FLOAT       Sensitivity parameter for Foldseek  [default:
                                9.5]
  --keep-tmp-files              Keep temporary intermediate files (e.g. large
                                foldseek_results.tsv of all Foldseek hits)
  --max-seqs INTEGER            Maximum results per query sequence allowed to
                                pass the prefilter (saves disk space for
                                enormous datasets)  [default: 1000]
  --ultra-sensitive             Run with maximum sensitivity by skipping
                                Foldseek prefilter (not recommended for large
                                datasets)
  --extra-foldseek-params TEXT  Extra Foldseek search params
  --custom-db TEXT              Path to custom database
  --foldseek-gpu                Enable Foldseek-GPU search acceleration
  --custom-annotations PATH     Custom Foldseek DB annotations (2 column tsv:
                                Foldseek headers, description)
  --euk                         Eukaryotic input genome.
  --fast                        Skips Foldseek search against AFDB Clusters.
  --cds-program TEXT            CDS prediction tool string for compliant
                                outputs if using non-Bakta or Prokka input
  --trna-program TEXT           tRNA prediction tool string for compliant
                                outputs if using non-Bakta or Prokka input
  --tmrna-program TEXT          tmRNA profile prediction tool string for
                                compliant outputs if using non-Bakta or Prokka
                                input
  --rrna-program TEXT           rRNA profile prediction tool string for
                                compliant outputs if using non-Bakta or Prokka
                                input
  --ncrna-program TEXT          ncRNA profile prediction tool string for
                                compliant outputs if using non-Bakta or Prokka
                                input
```

### `baktfold proteins`

* `baktfold proteins` runs `baktfold proteins-predict` and `baktfold proteins-compare` together in one command. Recommended if you are running `baktfold` with an NVIDIA GPU, along with `--foldseek-gpu`.
* Also recommended if you are running `baktfold` in a lower-resource environment without a GPU - you will also need to specify `--cpu`.
* See `baktfold predict` and `baktfold compare` for more details on specific parameters

* Note: You can also use a `bakta_proteins` output JSON file as input if you have already annotated your protein FASTA with `bakta`


Example usage where NVIDIA GPU is available:

```bash
baktfold run -i bakta.json -o baktfold_output -d baktfold_db -t 8 --foldseek-gpu
```

```bash
Usage: baktfold proteins [OPTIONS]

  baktfold proteins-predict then comapare all in one - GPU recommended

Options:
  -h, --help                     Show this message and exit.
  -V, --version                  Show the version and exit.
  -i, --input PATH               Path to input file (amino acid FASTA or Bakta 
                                  proteins output JSON format) [required]
  -o, --output PATH              Output directory   [default: output_baktfold]
  -t, --threads INTEGER          Number of threads  [default: 1]
  -p, --prefix TEXT              Output files' prefix  [default: baktfold]
  -d, --database TEXT            Path to Baktfold's database
  -f, --force                    Force overwrites output directory
  --autotune                     Run autotuning to detect and set best batch
                                 size for local hardware (recommended only for
                                 large dataset, e.g. thousands of proteins)
  --batch-size INTEGER           Batch size for ProstT5 (1 is usually fastest)
                                 [default: 1]
  --cpu                          Use CPU only.
  --omit-probs                   Do not output per residue 3Di probabilities
                                 from ProstT5
  --save-per-residue-embeddings  Save ProstT5 embeddings per resuide in a H5
                                 file
  --save-per-protein-embeddings  Save ProstT5 embeddings as means per protein
                                 in a H5 file
  --mask-threshold FLOAT         Mask 3Di residues below this value of ProstT5
                                 confidence for Foldseek searches  [default:
                                 25]
  -e, --evalue FLOAT             Evalue threshold for Foldseek  [default:
                                 1e-3]
  -s, --sensitivity FLOAT        Sensitivity parameter for Foldseek  [default:
                                 9.5]
  --keep-tmp-files               Keep temporary intermediate files (e.g. large
                                 foldseek_results.tsv of all Foldseek hits)
  --max-seqs INTEGER             Maximum results per query sequence allowed to
                                 pass the prefilter (saves disk space for
                                 enormous datasets)  [default: 1000]
  --ultra-sensitive              Run with maximum sensitivity by skipping
                                 Foldseek prefilter (not recommended for large
                                 datasets)
  --extra-foldseek-params TEXT   Extra Foldseek search params
  --custom-db TEXT               Path to custom database
  --foldseek-gpu                 Enable Foldseek-GPU search acceleration
  --custom-annotations PATH      Custom Foldseek DB annotations (2 column tsv:
                                 Foldseek headers, description)
  --euk                          Eukaryotic input genome.
  --fast                         Skips Foldseek search against AFDB Clusters.
  --cds-program TEXT             CDS prediction tool string for compliant
                                 outputs if using non-Bakta or Prokka input
  --trna-program TEXT            tRNA prediction tool string for compliant
                                 outputs if using non-Bakta or Prokka input
  --tmrna-program TEXT           tmRNA profile prediction tool string for
                                 compliant outputs if using non-Bakta or
                                 Prokka input
  --rrna-program TEXT            rRNA profile prediction tool string for
                                 compliant outputs if using non-Bakta or
                                 Prokka input
  --ncrna-program TEXT           ncRNA profile prediction tool string for
                                 compliant outputs if using non-Bakta or
                                 Prokka input
```

### `baktfold createdb`

This in an auxillary command that allows you to create a Foldseek compatible database from AA and 3Di protein sequences (such as those created by `baktfold predict`). 

Example usage 

```bash
baktfold createdb --fasta_aa proteins_aa.fasta  --fasta_3di baktfold_3di.fasta -o my_foldseek_db  
```

```bash
Usage: baktfold createdb [OPTIONS]

  Creates foldseek DB from AA FASTA and 3Di FASTA input files

Options:
  -h, --help             Show this message and exit.
  -V, --version          Show the version and exit.
  --fasta-aa PATH        Path to input file (amino acid FASTA format)
                         [required]
  --fasta-3di PATH       Path to input file (3Di FASTA format)  [required]
  -o, --output PATH      Output directory  [default:
                         output_baktfold_foldseek_db]
  -t, --threads INTEGER  Number of threads  [default: 1]
  -p, --prefix TEXT      Foldseek database prefix  [default:
                         baktfold_foldseek_db]
  -f, --force            Force overwrites output directory
```

### `baktfold convert-prokka`

* This is an auxillary command that enables you to convert a Prokka GenBank format annotation output file into a Bakta JSON file compatible as input for `baktfold run`

* Example usage

```bash
baktfold convert-prokka -i prokka.gbk -o prokka_baktfold_input.json 
```

```bash
Usage: baktfold convert-prokka [OPTIONS]

  Converts Prokka GenBank to Bakta format json

Options:
  -h, --help          Show this message and exit.
  -V, --version       Show the version and exit.
  -i, --input PATH    Path to Prokka input (GenBank: .gbk)  [required]
  -o, --outfile PATH  Output file (Bakta: .json)  [default:
                      converted_bakta_formatted.json]
  -f, --force         Force overwrites the output file
  --verbose           Verbose output
```

### `baktfold convert-euk`

* This is an auxillary command that enables you to convert a eukaryote GenBank format annotation output file into a Bakta JSON file compatible as input for `baktfold run` with `--euk`

* Example usage

```bash
baktfold convert-euk -i eukaryote.gbk -o euk_baktfold_input.json 
```

```bash
Usage: baktfold convert-euk [OPTIONS]

  (Experimental) Converts eukaryotic GenBank to Bakta format json

Options:
  -h, --help          Show this message and exit.
  -V, --version       Show the version and exit.
  -i, --input PATH    Path to Prokka input (GenBank: .gbk)  [required]
  -o, --outfile PATH  Output file (Bakta: .json)  [default:
                      converted_bakta_formatted.json]
  -f, --force         Force overwrites the output file
  --verbose           Verbose output
```
