"""
Module for manipulating genbank files
some taken from phynteny https://github.com/susiegriggo/Phynteny
"""

import binascii
import gzip
import multiprocessing.pool
from datetime import datetime
from pathlib import Path
from typing import IO, Dict, Union
from datetime import datetime

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord
from loguru import logger

from baktfold.utils.util import get_version


# imports


def open_protein_fasta_file(input_file: str) -> Union[IO[str], gzip.GzipFile]:
    """
    Open a fasta file, whether it is gzipped or plain text.

    Parameters:
    input_file (str): The path to the fasta file, either gzipped or plain.

    Returns:
    Union[IO[str], gzip.GzipFile]: A file handle to the opened fasta file.
    """
    input_file = Path(input_file)

    if input_file.suffix == ".gz":
        return gzip.open(input_file, "rt")
    else:
        return open(input_file, "r")


def is_gzip_file(f: Path) -> bool:
    """
    Method copied from Phispy see https://github.com/linsalrob/PhiSpy/blob/master/PhiSpyModules/helper_functions.py

    This is an elegant solution to test whether a file is gzipped by reading the first two characters.
    I also use a version of this in fastq_pair if you want a C version :)
    See https://stackoverflow.com/questions/3703276/how-to-tell-if-a-file-is-gzip-compressed for inspiration
    Args:
        f (Path): The file to test.

    Returns:
        bool: True if the file is gzip compressed, otherwise False.
    """
    with open(f, "rb") as i:
        return binascii.hexlify(i.read(2)) == b"1f8b"



def get_genbank(genbank: Path) -> dict:
    """
    Convert a GenBank file to a dictionary.

    This function reads a GenBank file and converts it into a dictionary.

    Args:
        genbank (Path): Path to the GenBank file.

    Returns:
        dict: A dictionary representation of the GenBank file.

    Raises:
        ValueError: If the provided file is not a GenBank file.
    """

    logger.info(f"Checking if input {genbank} is a Genbank format file")
    logger.info(f"If so, also detecting the likely input style out of Pharokka, Bakta and NCBI Refseq style.")
    def parse_records(handle):
        """
    Parses a genbank file and returns a list of SeqRecords.

    Args:
      file_path (str): The path to the genbank file to parse.
      file_format (str): The format of the genbank file. Defaults to 'genbank'.

    Returns:
      list: A list of SeqRecords parsed from the genbank file.

    Examples:
      >>> parse_records('example.gb')
      [SeqRecord(seq=Seq('ATGC'), id='example', name='example', description='example', dbxrefs=[]), ...]
    """
        try:
            records = list(SeqIO.parse(handle, "gb"))
            if not records:
                return {}, None
            gb_dict = {record.id: record for record in records}
            record = records[0]

            comment = record.annotations.get("comment", "")
            cds_feature = next((f for f in record.features if f.type == "CDS"), None)

            if cds_feature is None:
                logger.error(f"{genbank} appears to be a Genbank formatted file but no CDS was found. Please check your input.")
                return gb_dict, None

            # Check if 'Bakta' appears in the Comment - will appear there
            if "Bakta" in comment and "locus_tag" in cds_feature.qualifiers:
                logger.info(f"Detected Bakta style input Genbank. Using locus_tag qualifier from Bakta as the CDS IDs for Phold.")
                method = "Bakta"
            else:
                if "phrog" not in cds_feature.qualifiers and "protein_id" in cds_feature.qualifiers:
                    logger.info(f"Detected NCBI Refseq style input Genbank. Using protein_id qualifier as the CDS IDs for Phold.")
                    method = "NCBI"
                elif "phrog" in cds_feature.qualifiers and "ID" in cds_feature.qualifiers:
                    logger.info(f"Detected Pharokka style input Genbank. Using ID qualifier from Pharokka as the CDS IDs for Phold.")
                    method = "Pharokka"
                else:
                    logger.error(
                                f"Feature {cds_feature} could not be parsed. Therefore, the input style format for {genbank} could not be detected. Please check your input."
                            )
            return identify_long_ids(gb_dict), method
        except Exception as e:
            logger.warning(f"{genbank} is not a genbank file")
            return {}, None

    try:
        if is_gzip_file(genbank.strip()):
            with gzip.open(genbank.strip(), "rt") as handle:
                return parse_records(handle)
        else:
            with open(genbank.strip(), "rt") as handle:
                return parse_records(handle)
    except Exception as e:
        logger.warning(f"{genbank} is not a genbank file")
        return {}, None





def identify_long_ids(gb_dict: dict) -> dict:
    """

    Checks all feature IDs in gb_dict. If longer than 54 chars (line break from Pharokka/biopython reading GBK files), removes the space

    Args:
        dict: A dictionary representation of the GenBank file.

    Returns:
        dict: A dictionary representation of the GenBank file.
    """

    # remove spaces in ID/locus tag
    for record_id, record in gb_dict.items():
        for cds_feature in record.features:
            try:
                # if pharokka > 54 char IDs/locus tage, phold/biopython will parse with a space
                # no spaces in
                # for really long CDS IDs (over 54 chars), a space will be introduced
                # this is because the ID will go over a second line
                # weird bug noticed it on the Mgnify contigs annotated with Pharokka
                cds_id = cds_feature.qualifiers["ID"][0]
                if len(cds_id) >= 54:
                    logger.warning(
                        f"The CDS ID is {cds_id} is longer than 54 characters. It is recommended that you use short contig headers (which will therefore lead to shorter CDS ids)."
                    )
                    cds_feature.qualifiers["ID"][0] = cds_feature.qualifiers["ID"][
                        0
                    ].replace(" ", "")
            except:
                # will be GenBank/NCBI formatted
                # ID isn't a field and should be properly formatted - famous last words probably
                continue

    return gb_dict


def get_proteins(fasta: Path) -> dict:
    """
    Convert an Amino Acid FASTA file to a dictionary.

    This function reads a AA FASTA file and converts it into a dictionary.

    Args:
        fasta (Path): Path to the FASTA file.

    Returns:
        dict: A dictionary representation of the FASTA file.

    Raises:
        ValueError: If the provided file is not a FASTA file.
    """

    if is_gzip_file(fasta.strip()):
        try:
            fasta_dict = {}
            with gzip.open(fasta.strip(), "rt") as handle:
                sequence_id = ""
                sequence = ""
                for line in handle:
                    line = line.strip()
                    if line.startswith(">"):
                        if sequence_id:
                            fasta_dict[sequence_id] = sequence
                        sequence_id = line[1:]
                        sequence = ""
                    else:
                        sequence += line
                if sequence_id:
                    fasta_dict[sequence_id] = sequence
            handle.close()
        except ValueError:
            logger.error(f"{fasta.strip()} is not a FASTA file!")
            raise

    else:
        try:
            fasta_dict = {}
            with open(fasta.strip(), "rt", errors="ignore") as handle:
                sequence_id = ""
                sequence = ""
                for line in handle:
                    line = line.strip()
                    if line.startswith(">"):
                        if sequence_id:
                            fasta_dict[sequence_id] = sequence
                        sequence_id = line[1:]
                        sequence = ""
                    else:
                        sequence += line
                if sequence_id:
                    fasta_dict[sequence_id] = sequence
            handle.close()
        except ValueError:
            logger.error(f"{fasta.strip()} is not a FASTA file!")
            raise

    return fasta_dict
