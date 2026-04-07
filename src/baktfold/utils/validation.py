import shutil
import subprocess as sp
import sys
from pathlib import Path
from typing import Dict, Union

from Bio import SeqIO
from Bio.SeqUtils import IUPACData
from loguru import logger

from xopen import xopen

from baktfold.io.handle_genbank import  get_genbank



def instantiate_dirs(output_dir: Union[str, Path], force: bool) -> Path:
    """
    Checks and instantiates the output directory.

    Parameters:
        output_dir (Union[str, Path]): Path to the output directory.
        force (bool): Force flag indicating whether to overwrite existing directory.

    Returns:
        Path: Final output directory path.
    """

    # Checks the output directory
    # remove outdir on force
    logger.add(lambda _: sys.exit(1), level="ERROR")
    logger.info(f"Checking the output directory {output_dir}")
    if force is True:
        if Path(output_dir).exists():
            logger.info(f"Removing {output_dir} because --force was specified")
            shutil.rmtree(output_dir)
        else:
            logger.info(
                "--force was specified even though the output directory does not already exist. Continuing"
            )
    else:
        if Path(output_dir).exists():
            logger.error(
                "Output directory already exists and force was not specified. Please specify -f or --force to overwrite the output directory"
            )

    # instantiate outdir
    if Path(output_dir).exists() is False:
        Path(output_dir).mkdir(parents=True, exist_ok=True)

def validate_outfile(outfile: Union[str, Path], force: bool) -> Path:
    """
    Checks and instantiates the output file for baktfold convert-prokka

    Parameters:
        outfile (Union[str, Path]): Path to the output file.
        force (bool): Force flag indicating whether to overwrite existing outfile.

    Returns:
        Path: Final output file path.
    """

    # Checks the output directory
    # remove outdir on force
    logger.add(lambda _: sys.exit(1), level="ERROR")
    logger.info(f"Checking the output file {outfile}")
    if force is True:
        if Path(outfile).exists():
            logger.info(f"Removing {outfile} because --force was specified")
            Path(outfile).unlink()
        else:
            logger.info(
                f"--force was specified even though the output file {outfile} does not already exist. Continuing"
            )
    else:
        if Path(outfile).exists():
            logger.error(
                f"Output file {outfile} already exists and force was not specified. Please specify -f or --force to overwrite the output file"
            )


def check_dependencies() -> None:
    """
    Checks the dependencies and versions of non Python programs (i.e. Foldseek)

    Parameters:
        None

    Returns:
        None

    """

    #############
    # foldseek
    #############
    try:
        process = sp.Popen(["foldseek", "version"], stdout=sp.PIPE, stderr=sp.STDOUT)
    except:
        logger.error("Foldseek not found. Please reinstall baktfold.")

    foldseek_out, _ = process.communicate()
    foldseek_out = foldseek_out.decode()

    foldseek_version = foldseek_out.strip()

    if "941cd33" in foldseek_version:
        foldseek_major_version=10
        foldseek_minor_version="941cd33"
        logger.info(
        f"Foldseek version found is v{foldseek_major_version}.{foldseek_minor_version}"
    )
    else:
        logger.warning(f"Foldseek version found is v{foldseek_version}")
        logger.warning(f"baktfold is recommended to be run with Foldseek v10.941cd33")
        logger.warning(f"Using a different Foldseek version is likely to work without issue, but this cannot be guaranteed.")


    logger.info("Foldseek version is ok")

def check_genbank_and_prokka(filepath, euk):
    """
    Validate that an input file is a readable GenBank file and check whether it was
    annotated using Prokka. The function transparently supports compressed files
    (e.g., .gz, .bz2, .xz, .zst) via `xopen`.

    Validation steps:
      • Attempts to parse the file as GenBank using Biopython.
      • Logs an error and returns None if no GenBank records are found.
      • Checks the COMMENT field of each record for a Prokka signature
        ("Annotated using prokka", case-insensitive).
      • If no Prokka annotation is detected, a warning is logged but parsing continues as it is a valid genbank.

    Parameters
    ----------
    filepath : str
        Path to the GenBank or compressed GenBank file.
    euk: flag
        whether or not the input is eukaryotic (skips prokka)

    Returns
    -------
    list[SeqRecord] or None
        A list of Biopython SeqRecord objects if parsing succeeds.
        Returns None if the file is not valid GenBank or cannot be parsed.
    """

    logger.add(lambda _: sys.exit(1), level="ERROR")

    is_valid_genbank = False
    is_prokka = False

    try:
        # Use xopen so gzip/bz2/xz/zst work automatically
        with xopen(filepath, "rb") as handle:
            # SeqIO.parse expects text handle -> decode
            # Use .read() is too big; instead wrap in TextIOWrapper
            import io
            text_handle = io.TextIOWrapper(handle, encoding="utf-8", errors="replace")

            records = list(SeqIO.parse(text_handle, "genbank"))

        if not records:
            logger.error(f"Input file {filepath} is not GenBank format. Please check your input")
            return None
        else:
            is_valid_genbank = True


        # Scan comments for Prokka signature
        if not euk:
            for rec in records:
                comment = rec.annotations.get("comment", "") or ""
                if "annotated using prokka" in comment.lower():
                    is_prokka = True
                    break
                    
            
            if is_prokka is False:
                logger.warning(f"Input file {filepath} does not appear to come from Prokka.")
                logger.warning(f"Conversion will proceed but no guarantee of success.")

    except Exception:
        logger.error(f"There was an error parsing {filepath}. Please check your input")
        return None

    return records



def is_fasta(filepath: Path) -> bool:
    try:
        # Use xopen so gzip/bz2/xz/zst work automatically
        with xopen(filepath, "rb") as f:
            first_line = f.readline().strip()
            if not first_line.startswith(">"):
                return False

            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    continue
                if not all(c.isalpha() or c in "-*." for c in line):
                    return False

        return True
    except Exception:
        return False