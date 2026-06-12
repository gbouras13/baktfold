import os
import shutil
import sys
import tempfile
import time
from contextlib import contextmanager
from pathlib import Path
from typing import Any, Dict, Iterator, List, Union

from loguru import logger
from datetime import datetime

import baktfold.bakta.config as cfg
import baktfold.bakta.constants as bc
import click

from Bio import SeqIO


@contextmanager
def atomic_write_path(target: Union[str, Path]) -> Iterator[Path]:
    """Yield a sibling temp path that is renamed over ``target`` on success.

    On any exception (including KeyboardInterrupt), the temp is removed and
    ``target`` is left exactly as it was before the with-block.
    """
    target = Path(target)
    target.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(
        prefix=f".{target.name}.",
        suffix=".tmp",
        dir=str(target.parent),
    )
    os.close(fd)
    tmp_path = Path(tmp_name)
    try:
        yield tmp_path
    except BaseException:
        try:
            tmp_path.unlink()
        except FileNotFoundError:
            pass
        raise
    else:
        os.replace(tmp_path, target)


class OrderedCommands(click.Group):
    """This class will preserve the order of subcommands, which is useful when printing --help"""

    def list_commands(self, ctx: click.Context):
        """
        Returns a list of subcommands in the order they were added.

        Args:
          ctx (click.Context): The click context.

        Returns:
          list: A list of subcommands in the order they were added.
        """
        return list(self.commands)


def baktfold_base(rel_path):
    """
    Returns the absolute path to the given relative path.

    Args:
      rel_path (str): The relative path to the file.

    Returns:
      str: The absolute path to the file.
    """
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), rel_path)


def get_version():
    """
    Returns the version number from the VERSION file.

    Returns:
      str: The version number.
    """
    with open(baktfold_base("VERSION"), "r") as f:
        version = f.readline()
    return version


def echo_click(msg, log=None):
    """
    Prints a message to stdout and optionally to a log file.

    Args:
      msg (str): The message to print.
      log (str): The path to the log file.

    Returns:
      None
    """
    click.echo(msg, nl=False, err=True)
    if log:
        with open(log, "a") as lo:
            lo.write(msg)


def print_citation():
    """
    Prints the contents of the CITATION file to stdout.

    Returns:
      None
    """
    with open(baktfold_base("CITATION"), "r") as f:
        for line in f:
            echo_click(line)


log_fmt = (
    "[<green>{time:YYYY-MM-DD HH:mm:ss}</green>] <level>{level: <8}</level> | "
    "<level>{message}</level>"
)

# Module-level register of every loguru sink that ``begin_baktfold`` installed.
# loguru's logger is a process-wide singleton — without tracking, every
# re-invocation stacks a new file handler and a new ``sys.exit``-on-error
# handler on top of the previous ones, multiplying log output unboundedly.
_BAKTFOLD_SINK_IDS: List[int] = []


def _remove_baktfold_sinks() -> None:
    """Idempotently remove every sink installed by a prior begin_baktfold."""
    while _BAKTFOLD_SINK_IDS:
        sink_id = _BAKTFOLD_SINK_IDS.pop()
        try:
            logger.remove(sink_id)
        except ValueError:
            pass  # already removed elsewhere


"""
begin and end functions
"""


def begin_baktfold(params: Dict[str, Any], subcommand: str, no_log: bool = False) -> int:
    """
    Begin baktfold process.

    Parameters:
        params (Dict[str, Any]): A dictionary of parameters for baktfold.
        subcommand (str): Subcommand indicating the baktfold operation.
        no_log (bool): No log file

    Returns:
        int: Start time of the baktfold process.
    """
    # Tear down any sinks from a prior call before installing fresh ones.
    _remove_baktfold_sinks()

    # get start time
    start_time = time.time()

    cfg.run_start = datetime.now()

    # initial logging stuff — track ids so they can be removed in end_baktfold.
    if not no_log:
        log_file = os.path.join(params["--output"], f"baktfold_{subcommand}_{start_time}.log")
        _BAKTFOLD_SINK_IDS.append(logger.add(log_file))
    _BAKTFOLD_SINK_IDS.append(logger.add(lambda _: sys.exit(1), level="ERROR"))

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

    # Clean up sinks so a subsequent call (or test) starts with a clean logger.
    _remove_baktfold_sinks()


# need the logo here eventually
def print_splash():
    """
    Prints the splash screen to stdout.

    Returns:
      None
    """
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
        shutil.rmtree(dir_path, ignore_errors=True)


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
    result_tsv_cath: Path = Path(output) / "foldseek_results_cath.tsv"
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
    remove_file(result_tsv_cath)

def get_type_rank(f):
    """
    ranks eukaryotic features 1) in order of gene -> mRNA -> CDS and gene -> tRNA
    dynamically adjusts if 5'UTR and 3'UTR is present
    """
    t = f['type']
    strand = f.get('strand', '+')  # default to + if missing

    # fixed ranks
    base_order = {
        'gene': 0,
        'mRNA': 1,
        'cds': 3,
        'tRNA': 6
    }

    # dynamic UTR ordering
    if t == bc.FEATURE_5UTR:
        return 2 if strand == '+' else 4
    if t == bc.FEATURE_3UTR:
        return 4 if strand == '+' else 2

    return base_order.get(t, 99)   # non-protein features become 99


def sort_euk_feature_key(f):
    """
    Sorts a feature dictionary by start, locus, type rank, and stop.

    Args:
      f (dict): The feature dictionary.

    Returns:
      tuple: A tuple of the sorted values.
    """
    start = f.get('start', float('inf'))
    stop = f.get('stop', float('inf'))
    locus = f.get('locus')
    type_rank = get_type_rank(f)

    if locus and type_rank != 99:
        # Within a locus → sort by type rank second and stop last (if multiple CDS e.g.)
        return (start, 0, locus, type_rank, stop)
    else:
        # Non-locus or non-gene features → sort only by start
        return (start, 1, '', 99, stop)

def replace_pipe_in_fasta(input_path):
    """Replace '~PIPE~' with '|' in FASTA headers, writing atomically.

    Streams line-by-line to a sibling temp file and renames it onto
    ``input_path`` on success.  A kill mid-write leaves the original intact.
    """
    with atomic_write_path(input_path) as tmp:
        with open(input_path, "r") as in_f, open(tmp, "w") as out_f:
            for line in in_f:
                if line.startswith(">") and "~PIPE~" in line:
                    line = line.replace("~PIPE~", "|")
                out_f.write(line)