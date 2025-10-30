# baktfold
Rapid &amp; standardized annotation of bacterial genomes, MAGs &amp; plasmids using protein structural information

`baktfold` is a sensitive annotation tool for bacterial genomes, MAGs &amp; plasmids genomes using protein structural homology. 

`baktfold` is very similar to [phold](https://github.com/gbouras13/phold) but goes beyond phages to bacterial annotation. `baktfold` takes all "hypothetical proteins" from [Bakta's](https://github.com/oschwengers/bakta) output and uses the [ProstT5](https://github.com/mheinzinger/ProstT5) protein language model to rapidly translate protein amino acid sequences to the 3Di token alphabet used by [Foldseek](https://github.com/steineggerlab/foldseek). Foldseek is then used to search these against a series of databases (SwissProt, AlphaFold Database non-singleton clusters, PDB and CATH).

Additionally, instead of using ProstT5, you can specify protein structures that you have pre-computed for your hypothetical proteins.

You can also specify custom databases to search against using `--custom-db`

## Install

* We will be making baktfold available via pypi and bioconda shortly
* For now, please install from source
* The only non-python dependency is Foldseek

```bash
conda create -n baktfold foldseek
conda activate baktfold
git clone https://github.com/gbouras13/baktfold.git
cd baktfold
pip install .
baktfold --help
```

* If you have a Mac with M-series Apple Silicon, you may need to install a particular version of Pytorch to utilise GPU-acceleration
* The same is true if you use other non-NVIDIA e.g. AMD GPUs
* See [this link](https://phold.readthedocs.io/en/latest/install/) for some more detail and further links

* Database installation (use as many threads with `-t` as you can to speed up downloading)

```bash
baktfold install -d baktfold_db -t 8
```

* If you have an NVIDIA GPU, you will need to format the database to allow it to use Foldseek-GPU with `--foldseek-gpu`
    * Note: you can do this after downloading the database with the above command (it won't redownload the database, only do the relevant Foldseek database padding)

```bash
baktfold install -d baktfold_db --foldseek-gpu
```


## Example

* You will first need to run [bakta](https://github.com/oschwengers/bakta) and use the resulting `.json` file as input for `baktfold`
    * We will add other input formats eventually, put we would always recommend running Bakta first, as it is awesome and comprehensive
* To use `baktfold run` using a dummy test example

```bash
# with nvidia gpu 
baktfold run -i tests/tests_data/assembly_bakta_output/assembly.json  -o baktfold_output -f -t 8 -d baktfold_db   --foldseek-gpu
# without nvidia gpu available
baktfold run -i tests/tests_data/assembly_bakta_output/assembly.json  -o baktfold_output -f -t 8 -d baktfold_db   
```

* To use `baktfold proteins` using a dummy test example protein `.faa` file
* Note that this can be any `.faa` (It does not have to be the output of Bakta)

```bash
# with nvidia gpu 
baktfold proteins -i tests/tests_data/assembly.hypotheticals.faa  -o baktfold_proteins_output -f -t 8 -d baktfold_db   --foldseek-gpu
# without nvidia gpu available
baktfold proteins -i tests/tests_data/assembly.hypotheticals.faa  -o baktfold_proteins_output -f -t 8 -d baktfold_db   
```

## Usage

* The two most useful commands are `baktfold run` and `baktfold proteins`
* `baktfold run` accepts a __Bakta json file__ as its input, and by default, it will annotate all hypothetical CDS and return a variety of Bakta-like compliant output formats. All other annotations will be inherited from the Bakta output
* `baktfold proteins` accepts a protein FASTA `.faa` format file as input. It will annotate all protein sequences and return a variety of `bakta_proteins`-like output formats
* `baktfold predict` and `baktfold compare` split `baktfold run` into the ProstT5 and Foldseek modules, while `baktfold proteins-predict` and `baktfold proteins-compare` do the same for `baktfold proteins` (useful if you have non-NVIDIA GPUs)

* It is recommend you run baktfold with a GPU if you can.
* If you do not have a GPU, baktfold will still run, but the ProstT5 step will be fairly slow.
* If you have a NVIDIA GPU, you can also use the `--foldseek-gpu` parameter to accelerate Foldseek further

```bash
Usage: baktfold [OPTIONS] COMMAND [ARGS]...

Options:
  -h, --help     Show this message and exit.
  -V, --version  Show the version and exit.

Commands:
  citation          Print the citation(s) for this tool
  compare           Runs Foldseek vs baktfold db
  install           Installs ProstT5 model and baktfold database
  predict           Uses ProstT5 to predict 3Di tokens - GPU recommended
  proteins          baktfold protein-predict then comapare all in one - GPU...
  proteins-compare  Runs Foldseek vs baktfold db on proteins input
  proteins-predict  Runs ProstT5 on a multiFASTA input - GPU recommended
  run               baktfold predict then comapare all in one - GPU...
```

```bash
Usage: baktfold run [OPTIONS]

  baktfold predict then comapare all in one - GPU recommended

Options:
  -h, --help                     Show this message and exit.
  -V, --version                  Show the version and exit.
  -i, --input PATH               Path to input file in Bakta Genbank format or
                                 Bakta JSON format  [required]
  -o, --output PATH              Output directory   [default: output_baktfold]
  -t, --threads INTEGER          Number of threads  [default: 1]
  -p, --prefix TEXT              Prefix for output files  [default: baktfold]
  -d, --database TEXT            Specific path to installed baktfold database
  -f, --force                    Force overwrites the output directory
  --batch-size INTEGER           batch size for ProstT5. 1 is usually fastest.
                                 [default: 1]
  --cpu                          Use cpus only.
  --omit-probs                   Do not output per residue 3Di probabilities
                                 from ProstT5. Mean per protein 3Di
                                 probabilities will always be output.
  --save-per-residue-embeddings  Save the ProstT5 embeddings per resuide in a
                                 h5 file
  --save-per-protein-embeddings  Save the ProstT5 embeddings as means per
                                 protein in a h5 file
  --mask-threshold FLOAT         Masks 3Di residues below this value of
                                 ProstT5 confidence for Foldseek searches
                                 [default: 25]
  -e, --evalue FLOAT             Evalue threshold for Foldseek  [default:
                                 1e-3]
  -s, --sensitivity FLOAT        Sensitivity parameter for foldseek  [default:
                                 9.5]
  --keep-tmp-files               Keep temporary intermediate files,
                                 particularly the large foldseek_results.tsv
                                 of all Foldseek hits
  --max-seqs INTEGER             Maximum results per query sequence allowed to
                                 pass the prefilter. You may want to reduce
                                 this to save disk space for enormous datasets
                                 [default: 1000]
  --ultra-sensitive              Runs baktfold with maximum sensitivity by
                                 skipping Foldseek prefilter. Not recommended
                                 for large datasets.
  --extra-foldseek-params TEXT   Extra foldseek search params
  --custom-db TEXT               Path to custom database
  --foldseek-gpu                 Use this to enable compatibility with
                                 Foldseek-GPU search acceleration
  --custom-annotations PATH      Custom Foldseek DB annotations, 2 column tsv.
                                 Column 1 matches the Foldseek headers, column
                                 2 is the description.
```

## Citations

* A manuscript describing `baktfold` is in preparation.

* Please be sure to cite the following core dependencies - citing all bioinformatics tools that you use helps us, so helps you get better bioinformatics tools:

* Foldseek - (https://github.com/steineggerlab/foldseek) van Kempen M, Kim S, Tumescheit C, Mirdita M, Lee J, Gilchrist C, Söding J, and Steinegger M. Fast and accurate protein structure search with Foldseek. Nature Biotechnology (2023), [doi:10.1038/s41587-023-01773-0 ](https://www.nature.com/articles/s41587-023-01773-0)
* ProstT5 - (https://github.com/mheinzinger/ProstT5) Michael Heinzinger, Konstantin Weissenow, Joaquin Gomez Sanchez, Adrian Henkel, Martin Steinegger, Burkhard Rost. ProstT5: Bilingual language model for protein sequence and structure. NAR Genomics and Bioinformatics (2024) [doi:10.1101/2023.07.23.550085](https://doi.org/10.1093/nargab/lqae150) 

* Please also consider citing these databases where relevant:

* AFDB/SwissProt - Mihaly Varadi, Damian Bertoni, Paulyna Magana, Urmila Paramval, Ivanna Pidruchna, Malarvizhi Radhakrishnan, Maxim Tsenkov, Sreenath Nair, Milot Mirdita, Jingi Yeo, Oleg Kovalevskiy, Kathryn Tunyasuvunakool, Agata Laydon, Augustin Žídek, Hamish Tomlinson, Dhavanthi Hariharan, Josh Abrahamson, Tim Green, John Jumper, Ewan Birney, Martin Steinegger, Demis Hassabis, Sameer Velankar, AlphaFold Protein Structure Database in 2024: providing structure coverage for over 214 million protein sequences, Nucleic Acids Research, Volume 52, Issue D1, 5 January 2024, Pages D368–D375, [https://doi.org/10.1093/nar/gkad1011](https://doi.org/10.1093/nar/gkad1011)
* CATH - Orengo CA, Michie AD, Jones S, Jones DT, Swindells MB, Thornton JM. CATH--a hierarchic classification of protein domain structures. Structure. 1997 Aug 15;5(8):1093-108. doi: 10.1016/s0969-2126(97)00260-8. PMID: 9309224.
* PDB - H.M. Berman, J. Westbrook, Z. Feng, G. Gilliland, T.N. Bhat, H. Weissig, I.N. Shindyalov, P.E. Bourne, The Protein Data Bank (2000) Nucleic Acids Research 28: 235-242 [https://doi.org/10.1093/nar/28.1.235](https://doi.org/10.1093/nar/28.1.235)