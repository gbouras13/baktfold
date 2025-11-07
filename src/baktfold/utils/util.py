import os
import shutil
import sys
import time
from pathlib import Path
from typing import Any, Dict
import baktfold.bakta.config as cfg
import baktfold.bakta.constants as bc
import click
from loguru import logger
from datetime import datetime


class OrderedCommands(click.Group):
    """This class will preserve the order of subcommands, which is useful when printing --help"""

    def list_commands(self, ctx: click.Context):
        return list(self.commands)


def baktfold_base(rel_path):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), rel_path)


def get_version():
    with open(baktfold_base("VERSION"), "r") as f:
        version = f.readline()
    return version


def echo_click(msg, log=None):
    click.echo(msg, nl=False, err=True)
    if log:
        with open(log, "a") as lo:
            lo.write(msg)


def print_citation():
    with open(baktfold_base("CITATION"), "r") as f:
        for line in f:
            echo_click(line)


log_fmt = (
    "[<green>{time:YYYY-MM-DD HH:mm:ss}</green>] <level>{level: <8}</level> | "
    "<level>{message}</level>"
)

"""
begin and end functions
"""


def begin_baktfold(params: Dict[str, Any], subcommand: str) -> int:
    """
    Begin baktfold process.

    Parameters:
        params (Dict[str, Any]): A dictionary of parameters for baktfold.
        subcommand (str): Subcommand indicating the baktfold operation.

    Returns:
        int: Start time of the baktfold process.
    """
    # get start time
    start_time = time.time()

    cfg.run_start = datetime.now()

    # initial logging stuff
    log_file = os.path.join(params["--output"], f"baktfold_{subcommand}_{start_time}.log")
    # adds log file
    logger.add(log_file)
    logger.add(lambda _: sys.exit(1), level="ERROR")

    print_splash()
    logger.info("baktfold: rapid & standardized annotation of bacterial genomes, MAGs & plasmids using protein structural information")

    logger.info(f"You are using baktfold version {get_version()}")
    logger.info("Repository homepage is https://github.com/gbouras13/baktfold")
    logger.info(f"You are running baktfold {subcommand}")
    logger.info(f"Listing parameters")
    for key, value in params.items():
        logger.info(f"Parameter: {key} {value}")

    return start_time


def end_baktfold(start_time: float, subcommand: str) -> None:
    """
    Finish baktfold process and log elapsed time.

    Parameters:
        start_time (float): Start time of the process.
        subcommand (str): Subcommand name indicating the baktfold operation.

    Returns:
        None
    """

    # Determine elapsed time
    elapsed_time = time.time() - start_time
    elapsed_time = round(elapsed_time, 2)

    cfg.run_end = datetime.now()
    run_duration = (cfg.run_end - cfg.run_start).total_seconds()
    # logger.info(f'If you use these results please cite Baktfold: https://doi.org/{bc.BAKTA_DOI}')
    logger.info(f'If you use these results please cite Baktfold: https://github.com/gbouras13/baktfold')
    logger.info(f'baktfold {subcommand} successfully finished in {int(run_duration / 60):02}:{int(run_duration % 60):02} [mm:ss].')
   

    # Show elapsed time for the process
    logger.info(f"baktfold {subcommand} has finished")
    logger.info("Elapsed time: " + str(elapsed_time) + " seconds")


# need the logo here eventually
def print_splash():
    click.echo(
        """\b

  _           _    _    __      _     _ 
 | |         | |  | |  / _|    | |   | |
 | |__   __ _| | _| |_| |_ ___ | | __| |
 | '_ \ / _` | |/ / __|  _/ _ \| |/ _` |
 | |_) | (_| |   <| |_| || (_) | | (_| |
 |_.__/ \__,_|_|\_\\__|_| \___/|_|\__,_|
                                        
                                        
"""
    )


def remove_file(file_path: Path) -> None:
    """
    Remove a file if it exists.

    Parameters:
        file_path (Path): Path to the file to remove.

    Returns:
        None
    """
    if file_path.exists():
        file_path.unlink()  # Use unlink to remove the file


def remove_directory(dir_path: Path) -> None:
    """
    Remove a directory and all its contents if it exists.

    Parameters:
        dir_path (Path): Path to the directory to remove.

    Returns:
        None
    """
    if dir_path.exists():
        shutil.rmtree(dir_path)


def touch_file(path: Path) -> None:
    """
    Update the access and modification times of a file to the current time, creating the file if it does not exist.

    Parameters:
        path (Path): Path to the file.

    Returns:
        None
    """
    with open(path, "a"):
        os.utime(path, None)


def clean_up_temporary_files(output: Path, prefix: str) -> None:
    """
    Clean up temporary files generated during the baktfold process.

    Parameters:
        output (Path): Path to the output directory.
        prefix (str): prefix str


    Returns:
        None
    """
    
    baktfold_aa: Path = Path(output) / f"{prefix}_aa.fasta"
    result_tsv_swissprot: Path = Path(output) / "foldseek_results_swissprot.tsv"
    result_tsv_afdb: Path = Path(output) / "foldseek_results_afdb_clusters.tsv"
    result_tsv_pdb: Path = Path(output) / "foldseek_results_pdb.tsv"
    result_tsv_custom: Path = Path(output) / "foldseek_results_custom.tsv"
    foldseek_db: Path = Path(output) / "foldseek_db"
    result_db_base: Path = Path(output) / "result_db"
    temp_db: Path = Path(output) / "temp_db"
    
    remove_directory(result_db_base)
    remove_directory(temp_db)
    remove_directory(foldseek_db)

    remove_file(baktfold_aa)
    remove_file(result_tsv_swissprot)
    remove_file(result_tsv_afdb)
    remove_file(result_tsv_pdb)
    remove_file(result_tsv_custom)

