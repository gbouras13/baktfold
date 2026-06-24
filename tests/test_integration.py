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

# to run with euks
pytest --gpu-available --nvidia --euks --threads 8 .

"""

# import
import importlib.util
import json
import os
import shutil
import tempfile
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
test_prokka_output = Path("tests/test_data/assembly_prokka_output")
test_bakta_proteins_output = Path("tests/test_data/assembly_bakta_proteins_output")
database_dir = Path(f"{test_data}/baktfold_db")

# inputs
input_json: Path = f"{test_bakta_output}/assembly.json"
input_no_fs_hits_json: Path = f"{test_data}/SAMEA111266571.bakta.json"
input_fasta: Path = f"{test_data}/assembly.hypotheticals.faa"
input_proteins_json: Path = f"{test_data}/assembly_bakta_proteins_output_all/assembly.json"
input_pipe_fasta: Path = f"{test_data}/pipe.faa"
input_prok_gbk: Path = f"{test_prokka_output}/PROKKA_02192026.gbk"
input_prok_json: Path = f"{test_data}/assembly_prokka.json"
input_euk_gbk: Path = f"{test_data}/protist.gbk.gz"
input_ncbi_gbk: Path = f"{test_data}/clado.gbk.gz"
input_funannotate_gbk: Path = f"{test_data}/funannotate.gbk.gz"
input_fungi_gbk: Path = f"{test_data}/Aaosphaeria_arxii_cbs_175_79_gca_010015735.Aaoar1.62.nonchromosomal.gbk.gz"

pdb_dir = Path(f"{test_data}/pdbs")
cif_dir = Path(f"{test_data}/cifs")

output_dir = Path(f"{test_data}/test_outputs")
output_dir.mkdir(parents=True, exist_ok=True)

output_prok_json: Path = f"{output_dir}/assembly_prokka.json"
output_euk_json: Path = f"{output_dir}/protist.json"
output_funannotate_json: Path = f"{output_dir}/funannotate.json"
output_fungi_json: Path = f"{output_dir}/fungi.json"
output_ncbi_json: Path = f"{output_dir}/ncbi.json"

dummy_custom_db = Path(f"{test_data}/custom_db/dummy_custom_db")
dummy_custom_db_annotations = Path(f"{test_data}/custom_db/dummy_custom_db_annotations.tsv")

run_dir: Path = f"{output_dir}/run_json"
json_recon_dir: Path = f"{output_dir}/json_reconstruct"
json_recon_trna_dir: Path = f"{output_dir}/json_reconstruct_trna"
json_recon_proteins_dir: Path = f"{output_dir}/json_reconstruct_proteins"
run_prok_dir: Path = f"{output_dir}/run_prok_json"
run_euk_dir: Path = f"{output_dir}/run_protist_json"
run_funannotate_dir: Path = f"{output_dir}/run_funannotate_json"
run_fungi_dir: Path = f"{output_dir}/run_fungi_json"
run_ncbi_dir: Path = f"{output_dir}/run_ncbi_clado_json"


run_fast_dir: Path = f"{output_dir}/run_json_fast"
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
proteins_dir_from_json = f"{output_dir}/proteins_json"
proteins_pipe_dir: Path = f"{output_dir}/proteins_pipe"
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

@pytest.fixture(scope="session")
def euks(pytestconfig):
    return pytestconfig.getoption("euks")

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


# ── golden output regression ────────────────────────────────────────────────
# Each main test below calls assert_golden() on the output it just produced to
# compare it against a committed golden reference under
# tests/test_data/golden/<case>, using the same float-tolerant /
# timestamp-ignoring engine as run_comparison.sh (tests/compare_outputs.py).
# This adds golden regression coverage to every functionality with NO extra
# baktfold command runtime (it reuses the output the integration test already
# produced, after `baktfold install` has populated the database).
#
# Goldens are platform-sensitive (ProstT5 3Di prediction differs across
# hardware), so bootstrap / refresh them on the CI platform itself, then commit
# tests/test_data/golden/:
#
#     BAKTFOLD_UPDATE_GOLDEN=1 pytest tests/test_integration.py --gpu-available --threads 8
#
# Until a case's golden exists the comparison is skipped, so a not-yet-
# bootstrapped case never fails.
_co_spec = importlib.util.spec_from_file_location(
    "_baktfold_compare_outputs", Path(__file__).parent / "compare_outputs.py"
)
_compare_outputs = importlib.util.module_from_spec(_co_spec)
_co_spec.loader.exec_module(_compare_outputs)

GOLDEN_DIR = Path(test_data) / "golden"
UPDATE_GOLDEN = bool(os.environ.get("BAKTFOLD_UPDATE_GOLDEN"))


def assert_golden(produced, case, strict=False):
    """Compare a produced output dir (or single file) against the golden for ``case``.

    With BAKTFOLD_UPDATE_GOLDEN set, (re)writes the golden instead of comparing.
    Skips when the golden has not yet been bootstrapped.
    """
    produced = Path(produced)
    golden = GOLDEN_DIR / case

    if UPDATE_GOLDEN:
        if golden.exists():
            shutil.rmtree(golden)
        golden.mkdir(parents=True, exist_ok=True)
        if produced.is_dir():
            shutil.copytree(produced, golden, dirs_exist_ok=True)
        else:
            shutil.copy(produced, golden / produced.name)
        return

    if not golden.exists():
        pytest.skip(f"golden '{case}' not bootstrapped (set BAKTFOLD_UPDATE_GOLDEN=1 to create)")

    if produced.is_dir():
        diffs = _compare_outputs.compare_dirs(produced, golden, strict=strict)
    else:  # single-file output (e.g. convert-prokka): wrap for dir comparison
        with tempfile.TemporaryDirectory() as td:
            shutil.copy(produced, Path(td) / produced.name)
            diffs = _compare_outputs.compare_dirs(Path(td), golden, strict=strict)

    assert not diffs, f"golden mismatch for '{case}':\n" + "\n".join(diffs[:80])


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
    assert_golden(run_dir, "run")

def test_run_no_fs_hits(gpu_available, threads, nvidia):
    """test baktfold run on a genome with no foldseek hits for all dbs"""
    cmd = f"baktfold run -i {input_no_fs_hits_json} -o {run_dir} -t {threads} -d {database_dir} -e 1e-50 -f"
    if nvidia:
       cmd = f"{cmd} --foldseek-gpu" 
    if gpu_available is False:
        cmd = f"{cmd} --cpu"
    exec_command(cmd)
    assert_golden(run_dir, "run_no_fs_hits")

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
    assert_golden(run_all_dir, "run_all")

def test_run_fasta(gpu_available, threads, nvidia):
    """test baktfold run on all proteins just --fast"""
    cmd = f"baktfold run -i {input_json} -o {run_fast_dir} -t {threads} -d {database_dir} -f --fast"
    if nvidia:
       cmd = f"{cmd} --foldseek-gpu" 
    if gpu_available is False:
        cmd = f"{cmd} --cpu"
    exec_command(cmd)
    assert_golden(run_fast_dir, "run_fast")

def test_run_extra_foldseek_params(gpu_available, threads, nvidia):
    """test baktfold run on all proteins not just hyps with -a"""
    cmd = f"baktfold run -i {input_json} -o {run_dir_extra} -t {threads} -d {database_dir} -f --extra-foldseek-params \"--cov-mode 2\""
    if nvidia:
       cmd = f"{cmd} --foldseek-gpu" 
    if gpu_available is False:
        cmd = f"{cmd} --cpu"
    exec_command(cmd)
    assert_golden(run_dir_extra, "run_extra_foldseek_params")


def test_run_custom_db(gpu_available, threads, nvidia):
    """test baktfold run with custom db"""
    cmd = f"baktfold run -i {input_json} -o {run_dir_custom_db} -t {threads} -d {database_dir} --custom-db {dummy_custom_db} -f "
    if nvidia:
       cmd = f"{cmd} --foldseek-gpu" 
    if gpu_available is False:
        cmd = f"{cmd} --cpu"
    exec_command(cmd)
    assert_golden(run_dir_custom_db, "run_custom_db")

def test_run_custom_db_custom_annotations(gpu_available, threads, nvidia):
    """test baktfold run with custom db and custom db annotation tsv"""
    cmd = f"baktfold run -i {input_json} -o {run_dir_custom_db_custom_annotations} -t {threads} -d {database_dir} --custom-db {dummy_custom_db} --custom-annotations {dummy_custom_db_annotations} -f "
    if nvidia:
       cmd = f"{cmd} --foldseek-gpu" 
    if gpu_available is False:
        cmd = f"{cmd} --cpu"
    exec_command(cmd)
    assert_golden(run_dir_custom_db_custom_annotations, "run_custom_db_custom_annotations")



"""
predict tests
"""

def test_predict(gpu_available, threads, nvidia):
    """test baktfold predict"""
    cmd = f"baktfold predict -i {input_json} -o {predict_dir} -t {threads}  -d {database_dir} -f "
    if gpu_available is False:
        cmd = f"{cmd} --cpu"
    exec_command(cmd)
    assert_golden(predict_dir, "predict")


def test_predict_save_embeddings(gpu_available, threads, nvidia):
    """test baktfold predict and save embeddings"""
    cmd = f"baktfold predict -i {input_json} -o {predict_embeddings_dir} -t {threads}  -d {database_dir} -f --save-per-residue-embeddings --save-per-protein-embeddings"
    if gpu_available is False:
        cmd = f"{cmd} --cpu"
    exec_command(cmd)
    assert_golden(predict_embeddings_dir, "predict_embeddings")



"""
compare tests
"""

def test_compare(gpu_available, threads, nvidia):
    """test baktfold compare """
    cmd = f"baktfold compare -i {input_json} -o {compare_dir} --predictions-dir {predict_dir} -t {threads} -d {database_dir} -f"
    if nvidia:
        cmd = f"{cmd} --foldseek-gpu" 
    exec_command(cmd)
    assert_golden(compare_dir, "compare")


def test_compare_pdb(gpu_available, threads, nvidia):
    """test baktfold compare with pdbs input"""
    cmd = f"baktfold compare -i {input_json} -o {compare_pdb_dir} -t {threads} -d {database_dir} --structure-dir {pdb_dir} -f"
    if nvidia:
        cmd = f"{cmd} --foldseek-gpu" 
    exec_command(cmd)
    assert_golden(compare_pdb_dir, "compare_pdb")

def test_compare_cif(gpu_available, threads, nvidia):
    """test baktfold compare with cifs input"""
    cmd = f"baktfold compare -i {input_json} -o {compare_cif_dir} -t {threads} -d {database_dir} --structure-dir {cif_dir} -f"
    if nvidia:
        cmd = f"{cmd} --foldseek-gpu" 
    exec_command(cmd)
    assert_golden(compare_cif_dir, "compare_cif")

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
    assert_golden(proteins_dir, "proteins")

def test_proteins_json(gpu_available, threads, nvidia):
    """test baktfold proteins with json input"""
    cmd = f"baktfold proteins -i {input_proteins_json} -o {proteins_dir_from_json} -t {threads} -d {database_dir} -f"
    if nvidia:
       cmd = f"{cmd} --foldseek-gpu" 
    if gpu_available is False:
        cmd = f"{cmd} --cpu"
    exec_command(cmd)
    assert_golden(proteins_dir_from_json, "proteins_json")

def test_proteins_pipe(gpu_available, threads, nvidia):
    """test baktfold proteins where some inputs have | in header"""
    cmd = f"baktfold proteins -i {input_pipe_fasta} -o {proteins_pipe_dir} -t {threads} -d {database_dir} -f"
    if nvidia:
       cmd = f"{cmd} --foldseek-gpu" 
    if gpu_available is False:
        cmd = f"{cmd} --cpu"
    exec_command(cmd)
    assert_golden(proteins_pipe_dir, "proteins_pipe")

"""
proteins-predict
"""

def test_proteins_predict(gpu_available, threads, nvidia):
    """test baktfold proteins-predict"""
    cmd = f"baktfold proteins-predict -i {input_fasta} -o {proteins_predict_dir} -t {threads} -d {database_dir} -f"
    if gpu_available is False:
        cmd = f"{cmd} --cpu"
    exec_command(cmd)
    assert_golden(proteins_predict_dir, "proteins_predict")


"""
proteins-compare
"""

def test_proteins_compare(gpu_available, threads, nvidia):
    """test baktfold proteins-compare"""
    cmd = f"baktfold proteins-compare -i {input_fasta}  -o {proteins_compare_dir} --predictions-dir {proteins_predict_dir} -t {threads} -d {database_dir} -f"
    if nvidia:
       cmd = f"{cmd} --foldseek-gpu" 
    exec_command(cmd)
    assert_golden(proteins_compare_dir, "proteins_compare")

def test_proteins_compare_pdb(gpu_available, threads, nvidia):
    """test baktfold proteins-compare with pdbs input"""
    cmd = f"baktfold proteins-compare -i {input_fasta} -o {proteins_compare_pdb_dir} -t {threads} -d {database_dir} --structure-dir {pdb_dir} -f"
    if nvidia:
        cmd = f"{cmd} --foldseek-gpu" 
    exec_command(cmd)
    assert_golden(proteins_compare_pdb_dir, "proteins_compare_pdb")

def test_proteins_compare_cif(gpu_available, threads, nvidia):
    """test baktfold proteins-compare with cifs input"""
    cmd = f"baktfold proteins-compare -i {input_fasta} -o {proteins_compare_cif_dir} -t {threads} -d {database_dir} --structure-dir {cif_dir} -f"
    if nvidia:
        cmd = f"{cmd} --foldseek-gpu" 
    exec_command(cmd)
    assert_golden(proteins_compare_cif_dir, "proteins_compare_cif")

"""
autotune
"""


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

    cmd = f"baktfold autotune -t {threads} -d {database_dir}  --min-batch {min_batch}  --sample-seqs {sample_seqs} --max-batch {max_batch} --step {step}"
    if gpu_available is False:
        cmd = f"{cmd} --cpu"

    exec_command(cmd)

"""
prokka
"""

def test_convert_prokka(gpu_available, threads, nvidia):
    """test baktfold convert-prokka"""
    cmd = f"baktfold convert-prokka -i {input_prok_gbk} -o {output_prok_json} "
    exec_command(cmd)
    assert_golden(output_prok_json, "convert_prokka")

def test_json(gpu_available, threads, nvidia):
    """test baktfold json: reconstitute all genome outputs from a Bakta JSON (no database/foldseek required)"""
    cmd = f"baktfold json -i {input_json} -o {json_recon_dir} -f"
    exec_command(cmd)
    # every reconstitutable genome output is produced
    for ext in ("gff3", "gbff", "embl", "tsv", "inference.tsv", "faa", "ffn", "fna", "summary.txt", "json"):
        assert Path(f"{json_recon_dir}/baktfold.{ext}").exists()
    # the self-describing provenance block is persisted for faithful re-reconstruction
    data = json.loads(Path(f"{json_recon_dir}/baktfold.json").read_text())
    assert data.get("baktfold_run", {}).get("mode") == "genome"
    # INSDC genetic-code qualifier regression: must be /transl_table=11, not the
    # non-standard /translation_table or the boolean the old swapped args produced
    gbff = Path(f"{json_recon_dir}/baktfold.gbff").read_text()
    assert "/transl_table=11" in gbff
    assert "/translation_table" not in gbff

def test_json_trna_inference(gpu_available, threads, nvidia):
    """test baktfold json: bakta tRNA /inference is profile:tRNAscan:2.0 (matches the GFF source column)"""
    cmd = f"baktfold json -i {input_no_fs_hits_json} -o {json_recon_trna_dir} -f"
    exec_command(cmd)
    gbff = Path(f"{json_recon_trna_dir}/baktfold.gbff").read_text()
    assert "profile:tRNAscan:2.0" in gbff            # bakta tRNA inference
    assert "profile:tRNAscan-SE:2.0.12" not in gbff  # other_genbank program string must not leak

def test_json_proteins(gpu_available, threads, nvidia):
    """test baktfold json: reconstitute proteins outputs from a Bakta proteins JSON"""
    cmd = f"baktfold json -i {input_proteins_json} -o {json_recon_proteins_dir} -f"
    exec_command(cmd)
    # proteins mode emits exactly these
    for ext in ("tsv", "faa", "summary.txt", "json"):
        assert Path(f"{json_recon_proteins_dir}/baktfold.{ext}").exists()
    # and none of the genome-only formats
    for ext in ("gff3", "gbff", "embl", "ffn", "fna"):
        assert not Path(f"{json_recon_proteins_dir}/baktfold.{ext}").exists()
    data = json.loads(Path(f"{json_recon_proteins_dir}/baktfold.json").read_text())
    assert data.get("baktfold_run", {}).get("mode") == "proteins"

def test_run_prokka(gpu_available, threads, nvidia):
    """test baktfold run with prokka input"""
    cmd = f"baktfold run -i {output_prok_json} -o {run_prok_dir} -t {threads} -d {database_dir} -f"
    if nvidia:
       cmd = f"{cmd} --foldseek-gpu" 
    if gpu_available is False:
        cmd = f"{cmd} --cpu"
    exec_command(cmd)
    assert_golden(run_prok_dir, "run_prokka")


"""
euk
"""

def test_convert_euk(gpu_available, threads, nvidia, euks):
    """test baktfold convert-euk"""
    cmd = f"baktfold convert-euk -i {input_euk_gbk} -o {output_euk_json} -f"
    if euks:
        exec_command(cmd)
    else:
        pass

def test_run_euk(gpu_available, threads, nvidia, euks):
    """test baktfold run with euk input https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_000208925.2/"""
    cmd = f"baktfold run -i {output_euk_json} -o {run_euk_dir} -t {threads} -d {database_dir} -f --fast --euk"
    if nvidia:
       cmd = f"{cmd} --foldseek-gpu" 
    if gpu_available is False:
        cmd = f"{cmd} --cpu"
    if euks:
        exec_command(cmd)
    else:
        pass

"""
ensembl fungi 
https://ftp.ebi.ac.uk/ensemblgenomes/pub/release-62/fungi/genbank/fungi_ascomycota5_collection/aaosphaeria_arxii_cbs_175_79_gca_010015735/
wget https://ftp.ebi.ac.uk/ensemblgenomes/pub/release-62/fungi/genbank/fungi_ascomycota5_collection/aaosphaeria_arxii_cbs_175_79_gca_010015735/Aaosphaeria_arxii_cbs_175_79_gca_010015735.Aaoar1.62.nonchromosomal.dat.gz
"""

def test_convert_fungi(gpu_available, threads, nvidia, euks):
    """test baktfold convert-euk"""
    cmd = f"baktfold convert-euk -i {input_fungi_gbk} -o {output_fungi_json} -f"
    if euks:
        exec_command(cmd)
    else:
        pass

def test_run_fungi(gpu_available, threads, nvidia, euks):
    """test baktfold run with euk input https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_000208925.2/"""
    cmd = f"baktfold run -i {output_fungi_json} -o {run_fungi_dir} -t {threads} -d {database_dir} -f --euk"
    if nvidia:
       cmd = f"{cmd} --foldseek-gpu" 
    if gpu_available is False:
        cmd = f"{cmd} --cpu"
    if euks:
        exec_command(cmd)
    else:
        pass


"""
funannotate
https://ftp.ebi.ac.uk/ensemblgenomes/pub/release-62/fungi/genbank/fungi_ascomycota5_collection/aaosphaeria_arxii_cbs_175_79_gca_010015735/
"""

def test_convert_funannotate(gpu_available, threads, nvidia, euks):
    """test baktfold convert-euk"""
    cmd = f"baktfold convert-euk -i {input_funannotate_gbk} -o {output_funannotate_json} -f"
    if euks:
        exec_command(cmd)
    else:
        pass

def test_run_funannotate(gpu_available, threads, nvidia, euks):
    """test baktfold run with funannotate input"""
    cmd = f"baktfold run -i {output_funannotate_json} -o {run_funannotate_dir} -t {threads} -d {database_dir} -f --euk"
    if nvidia:
       cmd = f"{cmd} --foldseek-gpu" 
    if gpu_available is False:
        cmd = f"{cmd} --cpu"
    if euks:
        exec_command(cmd)
    else:
        pass

class testFails(unittest.TestCase):
    """Tests for fails"""

"""
NCBI Assembly
https://www.ncbi.nlm.nih.gov/datasets/genome/?bioproject=PRJEB55036
https://www.ncbi.nlm.nih.gov/bioproject/PRJEB55036
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/947/184/155/GCA_947184155.2_Cgoreaui_SCF055-01_v2.1/GCA_947184155.2_Cgoreaui_SCF055-01_v2.1_genomic.gbff.gz
"""

def test_download_genome(gpu_available, threads, nvidia, euks):
    """test baktfold run with euk input https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_000208925.2/"""
    url="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/947/184/155/GCA_947184155.2_Cgoreaui_SCF055-01_v2.1/GCA_947184155.2_Cgoreaui_SCF055-01_v2.1_genomic.gbff.gz"
    cmd = f"wget {url} -O {input_ncbi_gbk} "
    if euks:
        exec_command(cmd)
    else:
        pass

def test_convert_ncbi(gpu_available, threads, nvidia, euks):
    """test baktfold convert-euk"""
    cmd = f"baktfold convert-euk -i {input_ncbi_gbk} -o {output_ncbi_json} -f"
    if euks:
        exec_command(cmd)
    else:
        pass

def test_run_ncbi(gpu_available, threads, nvidia, euks):
    """test baktfold run with ncbi input"""
    cmd = f"baktfold run -i {output_ncbi_json} -o {run_ncbi_dir} -t {threads} -d {database_dir} -f --euk"
    if nvidia:
       cmd = f"{cmd} --foldseek-gpu" 
    if gpu_available is False:
        cmd = f"{cmd} --cpu"
    if euks:
        exec_command(cmd)
    else:
        pass

remove_directory(output_dir)
# remove_directory(database_dir)
