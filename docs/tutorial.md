# `baktfold` Tutorial

* This tutorial assumes you have [conda](https://github.com/conda-forge/miniforge) installed and the correct channels available.
* For a quick start (this should detect your)

```bash
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
sh ./Miniforge3-$(uname)-$(uname -m).sh
```

* This tutorial uses the _Dolichospermum compactum NIES-806_ (GCF_002368115.1) genome available at https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_002368115.1/ specifically at https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/368/115/GCF_002368115.1_ASM236811v1/GCF_002368115.1_ASM236811v1_genomic.fna.gz
* We will use the already Bakta annotated genome available in the Baktfold repository (detailed below) - you can skip Step 0 in this case

## Step 0 Run Bakta

* Intallation instructions for this step is not included (but is easily done via conda using `conda create -n bakta bakta`)

```bash
DB_PATH="~/bakta_db/db"
THREADS="8"

bakta   --db "$DB_PATH" \
        --output bakta_out \
        --threads "$THREADS" \
        --prefix "GCF_002368115" --force \
        "GCF_002368115.1_ASM236811v1_genomic.fna.gz"
```

## Step 1 Get the test data from 

```bash
curl -L -o GCF_002368115.json https://raw.githubusercontent.com/gbouras13/baktfold/main/tests/test_data/GCF_002368115.json
```

## Step 2 Installing and Running `Baktfold`

* To install `baktfold` with conda from bioconda (this is assuming you have an NVIDIA GPU available):

```bash
conda create -n baktfoldENV -c conda-forge -c bioconda baktfold pytorch=*=cuda*
conda activate baktfoldENV
baktfold install -t 8 --foldseek-gpu
```
* If you do not have an NVIDIA GPU available, remove `pytorch=*=cuda*` and `--foldseek-gpu`
* For more installation options including non-NVIDIA GPUs including Mac Mseries laptops, see the [installation documentation](https://gbouras13.github.io/baktfold/install/).

## Step 3 Running `baktfold`

```bash
baktfold run -i GCF_002368115.json -o baktfold_output -t 8 -p GCF_002368115 --foldseek-gpu
```

* If you do not have a GPU available:

```bash
baktfold run -i GCF_002368115.json -o baktfold_output -t 8 -p GCF_002368115 --cpu
```

* If you have a non-NVIDIA GPU available (e.g. Mac Mseries) and have installed the correct compatible version of PyTorch:

```bash
baktfold run -i GCF_002368115.json -o baktfold_output -t 8 -p GCF_002368115
```

