import hashlib
# import os
import shutil
import tarfile
from pathlib import Path

import requests
from alive_progress import alive_bar
from loguru import logger

from baktfold.utils.util import remove_directory
from baktfold.utils.external_tools import ExternalTool

# set this if changes
CURRENT_DB_VERSION: str = "0.0.1"

# to hold information about the different DBs
VERSION_DICTIONARY = {
    "0.0.1": {
        "md5": "b1eba2ac1a35e9c34b125887cb4aaf51",
        "major": 0,
        "minor": 0,
        "minorest": 1,
        "db_url": "https://zenodo.org/records/17347516/files/baktfold_db.tar.gz",
        "dir_name": "baktfold_db",
        "tarball": "baktfold_db.tar.gz",
        "prostt5_backup_url": "https://zenodo.org/records/11234657/files/models--Rostlab--ProstT5_fp16.tar.gz",
        "prostt5_backup_tarball": "models--Rostlab--ProstT5_fp16.tar.gz",
        "prostt5_backup_md5": "118c1997e6d2cb5025abda95d36681e0",
    },
}



BAKTFOLD_DB_NAMES = [
    "AFDBClusters",
    "AFDBClusters_ca",
    "AFDBClusters_ca.dbtype",
    "AFDBClusters_ca.index",
    "AFDBClusters.dbtype",
    "AFDBClusters_h",
    "AFDBClusters_h.index",
    "AFDBClusters_h.dbtype",
    "AFDBClusters.index",
    "AFDBClusters_ss",
    "AFDBClusters_ss.dbtype",
    "AFDBClusters_ss.index",
    "baktfold.db",
    "cath",
    "cath_ca",
    "cath_ca.dbtype",
    "cath_ca.index",
    "cath.dbtype",
    "cath_h",
    "cath_h.index",
    "cath_h.dbtype",
    "cath.index",
    "cath_ss",
    "cath_ss.dbtype",
    "cath_ss.index",
    "pdb",
    "pdb_ca",
    "pdb_ca.dbtype",
    "pdb_ca.index",
    "pdb.dbtype",
    "pdb_h",
    "pdb_h.index",
    "pdb_h.dbtype",
    "pdb.index",
    "pdb_ss",
    "pdb_ss.dbtype",
    "pdb_ss.index",
    "swissprot",
    "swissprot_ca",
    "swissprot_ca.dbtype",
    "swissprot_ca.index",
    "swissprot.dbtype",
    "swissprot_h",
    "swissprot_h.index",
    "swissprot_h.dbtype",
    "swissprot.index",
    "swissprot_ss",
    "swissprot_ss.dbtype",
    "swissprot_ss.index",
]


PROSTT5_MD5_DICTIONARY = {
    "refs": {"main": "962133e8e2bff04ec1768fa58dd788f3"},
    "blobs": {
        "2c19eb6e3b583f52d34b903b5978d3d30b6b7682": "8fd03e945174de0818746ecbde1aad8e",
        "60fe6bb247c90b8545d7b73820cd796ce6dcbd59": "ce32377619690072ced051ec7fc6c5f0",
        "6fc7be92c58e238f20a6cdea5a87b123a4ad35e2": "1deb27443d0d9b90b8e5704e752373e2",
        "74da7b4afcde53faa570114b530c726135bdfcdb813dec3abfb27f9d44db7324": "6ad28d9980aaec37d3935072d204e520",
        "b1a9ffcef73280cc57f090ad6446b4116b574b6c75d83ccc32778282f7f00855": "45f066fc3b0d87fb9f98bb0ddb97a3dc",
        "e9322396e6e75ecf8da41a9527e24dfa4eeea505": "b1cdd31ea50a37bf84cc0d7ef11820c8",
    },
}


baktfold_DB_FOLDSEEK_GPU_NAMES = [
    "AFDBClusters_gpu",
    "AFDBClusters_gpu_ca",
    "AFDBClusters_gpu_ca.dbtype",
    "AFDBClusters_gpu_ca.index",
    "AFDBClusters_gpu.dbtype",
    "AFDBClusters_gpu_h",
    "AFDBClusters_gpu_h.index",
    "AFDBClusters_gpu_h.dbtype",
    "AFDBClusters_gpu.index",
    "AFDBClusters_gpu_ss",
    "AFDBClusters_gpu_ss.dbtype",
    "AFDBClusters_gpu_ss.index",
    "cath_gpu",
    "cath_gpu_ca",
    "cath_gpu_ca.dbtype",
    "cath_gpu_ca.index",
    "cath_gpu.dbtype",
    "cath_gpu_h",
    "cath_gpu_h.index",
    "cath_gpu_h.dbtype",
    "cath_gpu.index",
    "cath_gpu_ss",
    "cath_gpu_ss.dbtype",
    "cath_gpu_ss.index",
    "pdb_gpu",
    "pdb_gpu_ca",
    "pdb_gpu_ca.dbtype",
    "pdb_gpu_ca.index",
    "pdb_gpu.dbtype",
    "pdb_gpu_h",
    "pdb_gpu_h.index",
    "pdb_gpu_h.dbtype",
    "pdb_gpu.index",
    "pdb_gpu_ss",
    "pdb_gpu_ss.dbtype",
    "pdb_gpu_ss.index",
    "swissprot_gpu",
    "swissprot_gpu_ca",
    "swissprot_gpu_ca.dbtype",
    "swissprot_gpu_ca.index",
    "swissprot_gpu.dbtype",
    "swissprot_gpu_h",
    "swissprot_gpu_h.index",
    "swissprot_gpu_h.dbtype",
    "swissprot_gpu.index",
    "swissprot_gpu_ss",
    "swissprot_gpu_ss.dbtype",
    "swissprot_gpu_ss.index",
]


def install_database(db_dir: Path, foldseek_gpu: bool, threads: int) -> None:
    """
    Install the baktfold database.

    Args:
        db_dir Path: The directory where the database should be installed.
        foldseek_gpu bool: Whether to install foldseek-gpu compatible baktfold db
        threads int: Number of threads available (makes downloading faster)
    """

    # check the database is installed
    logger.info(f"Checking baktfold database installation in {db_dir}.")
    downloaded_flag, gpu_flag = check_db_installation(db_dir, foldseek_gpu)
    if downloaded_flag:
        logger.info("All baktfold databases files are present")
    else:
        logger.info("Some baktfold databases files are missing")

        DICT = VERSION_DICTIONARY
        db_url = DICT[CURRENT_DB_VERSION]["db_url"]
        logger.info(f"Downloading baktfold DB from {db_url}")

        requiredmd5 = DICT[CURRENT_DB_VERSION]["md5"]
        tarball = DICT[CURRENT_DB_VERSION]["tarball"]
        
        tarball_path = Path(f"{db_dir}/{tarball}")
        logdir = Path(db_dir) / "logdir"

        try: 
            download(db_url, tarball_path, logdir, threads)
        except:
            logger.warning(
                f"Could not download file from Zenodo using aria2c! url={db_url}, path={tarball_path}"
            )
            logger.warning(f"Trying now with requests")
            download_requests(db_url, tarball_path)


        md5_sum = calc_md5_sum(tarball_path)

        if md5_sum == requiredmd5:
            logger.info(f"baktfold database file download OK: {md5_sum}")
        else:
            logger.error(
                f"Error: corrupt database file! MD5 should be '{requiredmd5}' but is '{md5_sum}'"
            )

        logger.info(
            f"Extracting baktfold database tarball: file={tarball_path}, output={db_dir}"
        )
        untar(tarball_path, db_dir, DICT)
        tarball_path.unlink()

    if foldseek_gpu:
        if gpu_flag:
            logger.info("All baktfold database files compatible with Foldseek-GPU are present")
        else:
            logger.info("Some baktfold database files compatible with Foldseek-GPU are missing")
            logger.info("Creating them")
            foldseek_makepaddedseqdb(db_dir)



"""
lots of this code from the marvellous bakta https://github.com/oschwengers/bakta, db.py specifically
"""


"""
aria2c bottlenecked by Zenodo but still faster than wget
dependency of Foldseek so it is always present
"""

def download(db_url: str, tarball_path: Path, logdir: Path, threads: int) -> None:
    """
    Download the database from the given URL using aria2c.

    Args:
        db_url (str): The URL of the database.
        tarball_path (Path): The path where the downloaded tarball should be saved.
        logdir (Path): The path to store logs
        threads (int): Number of threads for aria2c
    """

    cmd = (
        f'--dir {tarball_path.parent} '
        f'--out {tarball_path.name} '
        f'--max-connection-per-server={threads} '
        f'--allow-overwrite=true '
        f'{db_url}'
    )

    download_db = ExternalTool(
        tool="aria2c",
        input=f"",
        output=f"",
        params=f"{cmd}",
        logdir=logdir,
    )

    ExternalTool.run_download(download_db)

def download_requests(db_url: str, tarball_path: Path):
    """
    Downloads a file from a given URL using the requests library.

    Args:
      db_url (str): The URL of the file to download.
      tarball_path (Path): The path to save the downloaded file.

    Returns:
      None

    Examples:
      >>> download_requests("https://zenodo.org/records/17347516/files/baktfold_db.tar.gz", Path("baktfold_db.tar.gz"))
    """

    headers = {
        "User-Agent": f"baktfold/{CURRENT_DB_VERSION} (contact: george.bouras@adelaide.edu.au)"
    }

    try:
        with tarball_path.open("wb") as fh_out, requests.get(
            db_url, stream=True, headers=headers
        ) as resp:
            total_length = resp.headers.get("content-length")
            if total_length is not None:  # content length header is set
                total_length = int(total_length)
            with alive_bar(total=total_length, scale="SI") as bar:
                for data in resp.iter_content(chunk_size=1024 * 1024):
                    fh_out.write(data)
                    bar(count=len(data))
    except:
        logger.error(
            f"ERROR: Could not download file from Zenodo! url={db_url}, path={tarball_path}"
        )


def download_zenodo_prostT5(model_dir, logdir, threads):
    """
    Download the ProstT5 model from Zenodo

    Args:
        db_url (str): The URL of the database.
        tarball_path (Path): The path where the downloaded tarball should be saved.
    """

    db_url = VERSION_DICTIONARY[CURRENT_DB_VERSION]["prostt5_backup_url"]
    requiredmd5 = VERSION_DICTIONARY[CURRENT_DB_VERSION]["prostt5_backup_md5"]

    logger.info(f"Downloading ProstT5 model backup from {db_url}")

    tarball = VERSION_DICTIONARY[CURRENT_DB_VERSION]["prostt5_backup_tarball"]
    tarball_path = Path(f"{model_dir}/{tarball}")


    try: 
        download(db_url, tarball_path, logdir, threads)
    except IOError:
        logger.warning(
            f"Could not download file from Zenodo using aria2c! url={db_url}, path={tarball_path}"
        )
        logger.warning(f"Trying now with requests")
        download_requests(db_url, tarball_path)

    md5_sum = calc_md5_sum(tarball_path)

    if md5_sum == requiredmd5:
        logger.info(f"ProstT5 model backup file download OK: {md5_sum}")
    else:
        logger.error(
            f"Error: corrupt file! MD5 should be '{requiredmd5}' but is '{md5_sum}'"
        )

    logger.info(
        f"Extracting ProstT5 model backup tarball: file={tarball_path}, output={model_dir}"
    )

    try:
        with tarball_path.open("rb") as fh_in, tarfile.open(
            fileobj=fh_in, mode="r:gz"
        ) as tar_file:
            tar_file.extractall(path=str(model_dir))

    except OSError:
        logger.warning("Encountered OSError: {}".format(OSError))
        logger.error(f"Could not extract {tarball_path} to {model_dir}")

    tarball_path.unlink()


def check_prostT5_download(model_dir: Path, model_name: str) -> bool:
    """
     Args:
        model_dir (Path): Directory where the model and tokenizer is be stored.
        model_name (str): Name of the pre-trained T5 model.
    Returns:
        bool: bool to tell baktfold whether to download ProstT5
    """

    # assumes already has been downloaded
    download = False

    if model_name == "Rostlab/ProstT5_fp16":

        model_sub_dir = "models--Rostlab--ProstT5_fp16"
        DICT = PROSTT5_MD5_DICTIONARY


    for key in DICT:
        for nested_key in DICT[key]:
            file_path = Path(
                f"{model_dir}/{model_sub_dir}/{key}/{nested_key}"
            )

            # check file exists
            if file_path.exists():
                md5_sum = calc_md5_sum(file_path)
                if md5_sum != DICT[key][nested_key]:
                    logger.warning(
                        f"Corrupt model file {file_path}! MD5 should be '{DICT[key][nested_key]}' but is '{md5_sum}'"
                    )
                    download = True
            else:
                logger.warning(f"Model file {file_path} does not exist.")
                download = True
    
    return download


def calc_md5_sum(tarball_path: Path, buffer_size: int = 1024 * 1024) -> str:
    """
    Calculate the MD5 checksum of the given file.

    Args:
        tarball_path (Path): The path to the file for which the MD5 checksum needs to be calculated.
        buffer_size (int): The buffer size for reading the file.

    Returns:
        str: The MD5 checksum of the file.
    """

    md5 = hashlib.md5()
    with tarball_path.open("rb") as fh:
        data = fh.read(buffer_size)
        while data:
            md5.update(data)
            data = fh.read(buffer_size)
    return md5.hexdigest()


def untar(tarball_path: Path, output_path: Path, DICT: dict) -> None:
    """
    Extract the tarball to the output path.

    Args:
        tarball_path (Path): The path to the tarball file.
        output_path (Path): The path where the contents of the tarball should be extracted.
        DICT (dict): version dictionary
    """
    try:
        with tarball_path.open("rb") as fh_in, tarfile.open(
            fileobj=fh_in, mode="r:gz"
        ) as tar_file:
            tar_file.extractall(path=str(output_path))

        tarpath = Path(output_path) / DICT[CURRENT_DB_VERSION]["dir_name"]

        # Get a list of all files in the directory
        files_to_move = [f for f in tarpath.iterdir() if f.is_file()]

        # Move each file to the destination directory
        for file_name in files_to_move:
            destination_path = output_path / file_name.name
            shutil.move(file_name, destination_path)
        # remove the directory
        remove_directory(tarpath)

    except OSError:
        logger.warning("Encountered OSError: {}".format(OSError))
        logger.error(f"Could not extract {tarball_path} to {output_path}")


def check_db_installation(db_dir: Path, foldseek_gpu: bool) -> bool:
    """
    Check if the baktfold database is installed.

    Args:
        db_dir Path: The directory where the database is installed.
        foldseek_gpu bool: Whether to install foldseek-gpu compatible baktfold db

    Returns:
        bool: True if all required files are present, False otherwise.
    """
    downloaded_flag = True
    for file_name in BAKTFOLD_DB_NAMES:
        path = Path(db_dir) / file_name
        if not path.is_file():
            logger.warning(f"baktfold Database file {path} is missing")
            downloaded_flag = False
            break
    
    gpu_flag = True
    if foldseek_gpu:
        for file_name in baktfold_DB_FOLDSEEK_GPU_NAMES:
            path = Path(db_dir) / file_name
            if not path.is_file():
                logger.warning(f"baktfold Foldseek-GPU Database file {path} is missing")
                gpu_flag = False
                break 

    return downloaded_flag, gpu_flag


def validate_db(database: str, default_dir: str, foldseek_gpu: bool) -> Path:
    """
    Validates the baktfold database is installed.

    Args:
        database str: The directory where the database is installed.
        default_dir str: Default DB location
        foldseek_gpu bool: Whether to install foldseek-gpu compatible baktfold db

    Returns:
        bool: True if all required files are present, False otherwise.
    """
    # set default DB if not specified
    if database is not None:
        database: Path = Path(database)
    else:
        database = Path(default_dir)

    # check the database is installed
    logger.info(f"Checking baktfold database installation in {database}")
    downloaded_flag, gpu_flag = check_db_installation(database, foldseek_gpu)
    if downloaded_flag == True:
        logger.info("All baktfold databases files are present")
    else:
        if database == Path(default_dir):  # default
            logger.error(
                f"baktfold database not found. Please run baktfold install to download and install the baktfold database"
            )
        else:  # specific
            logger.error(
                f"baktfold database not found. Please run baktfold install -d {database} to download and install the baktfold database"
            )
    if foldseek_gpu:
        if gpu_flag:
            logger.info("All baktfold database files compatible with Foldseek-GPU are present")
        else:
            logger.error(
                f"baktfold database files compatible with Foldseek-GPU not found. Please run baktfold install -d {database} --foldseek_gpu"
            )


    return database

def foldseek_makepaddedseqdb(db_dir: Path) -> None:
    """
    Runs the Foldseek makepaddedseqdb command on a given database directory.

    Args:
      db_dir (Path): The path to the database directory.

    Returns:
      None

    Examples:
      >>> foldseek_makepaddedseqdb(Path("baktfold_db"))
    """

    dbs = ["AFDBClusters", "pdb", "cath", "swissprot"]
    logdir = Path(db_dir) / "logdir"

    for db_name in dbs:
        db_path = Path(db_dir) / db_name
        db_path_gpu = Path(db_dir) / f"{db_name}_gpu"

        foldseek_makepaddedseqdb = ExternalTool(
            tool="foldseek",
            input="",
            output="",
            params=f"makepaddedseqdb {db_path} {db_path_gpu}",
            logdir=logdir,
        )

        ExternalTool.run_tool(foldseek_makepaddedseqdb)



    
