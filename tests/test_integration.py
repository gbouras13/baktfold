"""
Integration tests for baktfold

# to run pytest without remote and no gpu
pytest .


# to run with gpu
pytest  --gpu-available .

# to run with NVIDIA gpu available
pytest  --gpu-available --nvidia .

# to run with 8 threads 
pytest --gpu-available --nvidia --threads 8 .

"""

# import
import os
import shutil
# import functions
import subprocess
import sys
import unittest
from pathlib import Path
from unittest.mock import patch

import pytest
from loguru import logger

# import functions


# test data
test_data = Path("tests/test_data")
test_bakta_output = Path("tests/test_data/assembly_bakta_output")
test_bakta_proteins_output = Path("tests/test_data/assembly_bakta_proteins_output")
database_dir = Path(f"{test_data}/baktfold_db")

# inputs
input_json: Path = f"{test_bakta_output}/assembly.json"
input_fasta: Path = f"{test_data}/assembly.hypotheticals.faa"


pdb_dir = Path(f"{test_data}/pdbs")
cif_dir = Path(f"{test_data}/cifs")

output_dir = Path(f"{test_data}/test_outputs")
output_dir.mkdir(parents=True, exist_ok=True)

dummy_custom_db = Path(f"{test_data}/custom_db/dummy_custom_db")
dummy_custom_db_annotations = Path(f"{test_data}/custom_db/dummy_custom_db_annotations.tsv")

run_dir: Path = f"{output_dir}/run_json"
run_all_dir: Path = f"{output_dir}/run_json_all"
run_dir_extra: Path = f"{output_dir}/run_json_extra"
run_dir_custom_db: Path = f"{output_dir}/run_json_custom_db"
run_dir_custom_db_custom_annotations: Path = f"{output_dir}/run_json_custom_db_cuustom_annotations"

predict_dir: Path = f"{output_dir}/predict_json"
predict_embeddings_dir: Path = f"{output_dir}/predict_embeddings_json"

compare_dir: Path = f"{output_dir}/compare_json"
compare_pdb_dir: Path = f"{output_dir}/compare_pdb_json"
compare_cif_dir: Path = f"{output_dir}/compare_cif_json"

proteins_dir: Path = f"{output_dir}/proteins"

proteins_predict_dir: Path = f"{output_dir}/proteins_predict"

proteins_compare_dir: Path = f"{output_dir}/proteins_compare"
proteins_compare_pdb_dir: Path = f"{output_dir}/proteins_compare_pdb_json"
proteins_compare_cif_dir: Path = f"{output_dir}/proteins_compare_cif_json"


logger.add(lambda _: sys.exit(1), level="ERROR")
# threads = 1

def remove_directory(dir_path):
    if os.path.exists(dir_path):
        shutil.rmtree(dir_path)

@pytest.fixture(scope="session")
def gpu_available(pytestconfig):
    return pytestconfig.getoption("gpu_available")

@pytest.fixture(scope="session")
def nvidia(pytestconfig):
    return pytestconfig.getoption("nvidia")

@pytest.fixture(scope="session")
def threads(pytestconfig):
    return pytestconfig.getoption("threads")


def exec_command(cmnd, stdout=subprocess.PIPE, stderr=subprocess.PIPE):
    """executes shell command and returns stdout if completes exit code 0
    Parameters
    ----------
    cmnd : str
      shell command to be executed
    stdout, stderr : streams
      Default value (PIPE) intercepts process output, setting to None
      blocks this."""

    proc = subprocess.Popen(cmnd, shell=True, stdout=stdout, stderr=stderr)
    out, err = proc.communicate()
    if proc.returncode != 0:
        raise RuntimeError(f"FAILED: {cmnd}\n{err}")
    return out.decode("utf8") if out is not None else None

"""
install tests
"""

def test_install(threads, nvidia):
    """test baktfold install"""
    cmd = f"baktfold install -d {database_dir} -t {threads}"
    if nvidia:
       cmd = f"{cmd} --foldseek-gpu" 
    exec_command(cmd)

"""
run tests
"""

def test_run(gpu_available, threads, nvidia):
    """test baktfold run"""
    cmd = f"baktfold run -i {input_json} -o {run_dir} -t {threads} -d {database_dir} -f"
    if nvidia:
       cmd = f"{cmd} --foldseek-gpu" 
    if gpu_available is False:
        cmd = f"{cmd} --cpu"
    exec_command(cmd)

def test_run_autotune(gpu_available, threads):
    """test baktfold run with --autotune"""
    cmd = f"baktfold run -i {input_json} -o {run_dir} -t {threads} -d {database_dir} -f --autotune"
    if gpu_available is False:
        cmd = f"{cmd} --cpu"
    exec_command(cmd)

def test_run_all(gpu_available, threads, nvidia):
    """test baktfold run on all proteins not just hyps with -a"""
    cmd = f"baktfold run -i {input_json} -o {run_all_dir} -t {threads} -d {database_dir} -f -a"
    if nvidia:
       cmd = f"{cmd} --foldseek-gpu" 
    if gpu_available is False:
        cmd = f"{cmd} --cpu"
    exec_command(cmd)

def test_run_extra_foldseek_params(gpu_available, threads, nvidia):
    """test baktfold run on all proteins not just hyps with -a"""
    cmd = f"baktfold run -i {input_json} -o {run_dir_extra} -t {threads} -d {database_dir} -f --extra-foldseek-params \"--cov-mode 2\""
    if nvidia:
       cmd = f"{cmd} --foldseek-gpu" 
    if gpu_available is False:
        cmd = f"{cmd} --cpu"
    exec_command(cmd)


def test_run_custom_db(gpu_available, threads, nvidia):
    """test baktfold run with custom db"""
    cmd = f"baktfold run -i {input_json} -o {run_dir_custom_db} -t {threads} -d {database_dir} --custom-db {dummy_custom_db} -f "
    if nvidia:
       cmd = f"{cmd} --foldseek-gpu" 
    if gpu_available is False:
        cmd = f"{cmd} --cpu"
    exec_command(cmd)

def test_run_custom_db_custom_annotations(gpu_available, threads, nvidia):
    """test baktfold run with custom db and custom db annotation tsv"""
    cmd = f"baktfold run -i {input_json} -o {run_dir_custom_db_custom_annotations} -t {threads} -d {database_dir} --custom-db {dummy_custom_db} --custom-annotations {dummy_custom_db_annotations} -f "
    if nvidia:
       cmd = f"{cmd} --foldseek-gpu" 
    if gpu_available is False:
        cmd = f"{cmd} --cpu"
    exec_command(cmd)

"""
predict tests
"""

def test_predict(gpu_available, threads, nvidia):
    """test baktfold predict"""
    cmd = f"baktfold predict -i {input_json} -o {predict_dir} -t {threads}  -d {database_dir} -f "
    if gpu_available is False:
        cmd = f"{cmd} --cpu"
    exec_command(cmd)


def test_predict_save_embeddings(gpu_available, threads, nvidia):
    """test baktfold predict and save embeddings"""
    cmd = f"baktfold predict -i {input_json} -o {predict_embeddings_dir} -t {threads}  -d {database_dir} -f --save-per-residue-embeddings --save-per-protein-embeddings"
    if gpu_available is False:
        cmd = f"{cmd} --cpu"
    exec_command(cmd)



"""
compare tests
"""

def test_compare(gpu_available, threads, nvidia):
    """test baktfold compare """
    cmd = f"baktfold compare -i {input_json} -o {compare_dir} --predictions-dir {predict_dir} -t {threads} -d {database_dir} -f"
    if nvidia:
        cmd = f"{cmd} --foldseek-gpu" 
    exec_command(cmd)


def test_compare_pdb(gpu_available, threads, nvidia):
    """test baktfold compare with pdbs input"""
    cmd = f"baktfold compare -i {input_json} -o {compare_pdb_dir} -t {threads} -d {database_dir} --structure-dir {pdb_dir} -f"
    if nvidia:
        cmd = f"{cmd} --foldseek-gpu" 
    exec_command(cmd)

def test_compare_cif(gpu_available, threads, nvidia):
    """test baktfold compare with cifs input"""
    cmd = f"baktfold compare -i {input_json} -o {compare_cif_dir} -t {threads} -d {database_dir} --structure-dir {cif_dir} -f"
    if nvidia:
        cmd = f"{cmd} --foldseek-gpu" 
    exec_command(cmd)

"""
proteins 
"""


def test_proteins(gpu_available, threads, nvidia):
    """test baktfold proteins"""
    cmd = f"baktfold proteins -i {input_fasta} -o {proteins_dir} -t {threads} -d {database_dir} -f"
    if nvidia:
       cmd = f"{cmd} --foldseek-gpu" 
    if gpu_available is False:
        cmd = f"{cmd} --cpu"
    exec_command(cmd)


"""
proteins-predict
"""

def test_proteins_predict(gpu_available, threads, nvidia):
    """test baktfold proteins-predict"""
    cmd = f"baktfold proteins-predict -i {input_fasta} -o {proteins_predict_dir} -t {threads} -d {database_dir} -f"
    if gpu_available is False:
        cmd = f"{cmd} --cpu"
    exec_command(cmd)


"""
proteins-compare
"""

def test_proteins_compare(gpu_available, threads, nvidia):
    """test baktfold proteins-compare"""
    cmd = f"baktfold proteins-compare -i {input_fasta}  -o {proteins_compare_dir} --predictions-dir {proteins_predict_dir} -t {threads} -d {database_dir} -f"
    if nvidia:
       cmd = f"{cmd} --foldseek-gpu" 
    exec_command(cmd)

def test_proteins_compare_pdb(gpu_available, threads, nvidia):
    """test baktfold proteins-compare with pdbs input"""
    cmd = f"baktfold proteins-compare -i {input_fasta} -o {proteins_compare_pdb_dir} -t {threads} -d {database_dir} --structure-dir {pdb_dir} -f"
    if nvidia:
        cmd = f"{cmd} --foldseek-gpu" 
    exec_command(cmd)

def test_proteins_compare_cif(gpu_available, threads, nvidia):
    """test baktfold proteins-compare with cifs input"""
    cmd = f"baktfold proteins-compare -i {input_fasta} -o {proteins_compare_cif_dir} -t {threads} -d {database_dir} --structure-dir {cif_dir} -f"
    if nvidia:
        cmd = f"{cmd} --foldseek-gpu" 
    exec_command(cmd)


def test_autotune(gpu_available, threads, nvidia):
    """test autotune"""

    if gpu_available:
        min_batch = 1
        sample_seqs = 602
        max_batch = 301
        step = 20
    else:
        min_batch = 1
        sample_seqs = 10
        max_batch = 10
        step = 9

    cmd = f"baktfold autotune -t {threads} -d {database_dir}  --min_batch {min_batch}  --sample_seqs {sample_seqs} --max_batch {max_batch} --step {step}"
    if gpu_available is False:
        cmd = f"{cmd} --cpu"

    exec_command(cmd)




# class testFails(unittest.TestCase):
#     """Tests for fails"""
   


remove_directory(output_dir)
# remove_directory(database_dir)
