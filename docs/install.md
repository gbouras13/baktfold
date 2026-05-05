# Installation

The best way to install `baktfold` is using conda, so you can install [Foldseek](https://github.com/steineggerlab/foldseek) (the only non-Python dependency) along with the Python dependencies.

We would highly recommend installing conda via [miniforge](https://github.com/conda-forge/miniforge). For more information, please follow the installation instructions  (https://github.com/conda-forge/miniforge). To install the latest version on your system (your operating system and architecture will automatically be detected):

```bash
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
sh ./Miniforge3-$(uname)-$(uname -m).sh
```

We then detail installing Baktfold in four different scenarios

1. Your machine does not have a GPU
2. Your machine has a modern NVIDIA GPU installed
3. Your machine is manufactured by Apple and has an Apple Silicon chip 
4. Your machine has a modern GPU installed from a different vendor (e.g. AMD)


1. If you have no GPU, the default bioconda install (CPU-only) will suffice:

```bash
conda create -n baktfoldENV -c conda-forge -c bioconda baktfold
```

2. If your machine has an NVIDIA GPU, you can create a conda environment with a CUDA-compatible PyTorch installation in one line automatically:

```bash
conda create -n baktfoldENV -c conda-forge -c bioconda baktfold pytorch=*=cuda*
```

3. If your Apple machine has an Apple Silicon chip with integrated GPU cores (e.g. M-series), you can create a conda environment with a GPU compatible PyTorch version as follows by installing from the PyTorch channel:
  
```bash
conda create -n baktfoldENV 
conda activate baktfoldENV
conda install -c pytorch pytorch torchvision torchaudio
conda install -c conda-forge -c bioconda baktfold
```

4. If you have a GPU from a different vendor (e.g. AMD) or have older hardware (e.g. CUDA drivers), we direct you to follow the instructions at https://pytorch.org and https://pytorch.org/get-started/previous-versions/ to install the appropriate version of PyTorch for your machine using pip . For example, to install  baktfold on a machine with an AMD GPU and ROCM 7.2

```bash
conda create -n baktfoldENV pip
conda activate baktfoldENV
pip install torch torchvision --index-url https://download.pytorch.org/whl/rocm7.2
conda install -c conda-forge -c bioconda baktfold
```



## Pip

You can also install `baktfold` using pip.

```bash
pip install baktfold
```

You will need to have [Foldseek](https://github.com/steineggerlab/foldseek) (ideally v10.941cd33) installed and available in the $PATH.

## Source

You can install the latest version of `baktfold` with potentially untested and unreleased changes into a conda environment using [conda](https://github.com/conda-forge/miniforge) as follows:

```bash
conda create -n baktfoldENV pip foldseek python=3.13
conda activate baktfoldENV
git clone https://github.com/gbouras13/baktfold.git
cd baktfold 
pip install -e .
```


# Database Installation

To download and install the `baktfold` database:

```bash
baktfold install -t <threads>
```

If you would like to specify a particular location for the database (e.g. if you use `baktfold` on a shared server), please use `-d`

```bash
baktfold install -d <path/to/databse_dir> -t <threads>
```

* Note: You will need at least 15GB of free space.

If you have an NVIDIA GPU available, you may wish to accelerate Foldseek using GPU. To do this, you will need to format the databases appropriately as follows

```bash
baktfold install -d <path/to/databse_dir> --foldseek-gpu -t <threads>
```

```bash
Usage: baktfold install [OPTIONS]

  Installs ProstT5 model and baktfold database

Options:
  -h, --help             Show this message and exit.
  -V, --version          Show the version and exit.
  -d, --database TEXT    Path to install Baktfold's database
  --foldseek-gpu         Enable Foldseek-GPU acceleration
  -t, --threads INTEGER  Number of threads  [default: 1]
```

