#!/usr/bin/env python3
"""

Some code adapted from @mheinzinger 

https://github.com/mheinzinger/ProstT5/blob/main/scripts/generate_foldseek_db.py

"""


import os
import shutil
from pathlib import Path

from baktfold.utils.external_tools import ExternalTool
from baktfold.utils.util import remove_file

from Bio import SeqIO
from loguru import logger


def generate_foldseek_db_from_aa_3di(
    fasta_aa: Path, fasta_3di: Path, foldseek_db_path: Path, logdir: Path, prefix: str
) -> None:
    """
    Generate Foldseek database from amino-acid and 3Di sequences.

    Args:
        fasta_aa (Path): Path to the amino-acid FASTA file.
        fasta_3di (Path): Path to the 3Di FASTA file.
        foldseek_db_path (Path): Path to the directory where Foldseek database will be stored.
        logdir (Path): Path to the directory where logs will be stored.
        prefix (str): Prefix for the Foldseek database.

    Returns:
        None
    """
    # read in amino-acid sequences
    sequences_aa = {}
    for record in SeqIO.parse(fasta_aa, "fasta"):
        sequences_aa[record.id] = str(record.seq)

    # read in 3Di strings
    sequences_3di = {}
    for record in SeqIO.parse(fasta_3di, "fasta"):
        if not record.id in sequences_aa.keys():
            logger.warning(
                "Warning: ignoring 3Di entry {}, since it is not in the amino-acid FASTA file".format(
                    record.id
                )
            )
        else:
            sequences_3di[record.id] = str(record.seq)  #no upper if masked

    # assert that we parsed 3Di strings for all sequences in the amino-acid FASTA file
    for id in sequences_aa.keys():
        if not id in sequences_3di.keys():
            logger.warning(
                "Warning: entry {} in amino-acid FASTA file has no corresponding 3Di string".format(
                    id
                )
            )
            logger.warning("Removing: entry {} from the Foldseek database ".format(id))
            sequences_aa = {
                id: sequence
                for id, sequence in sequences_aa.items()
                if id in sequences_3di
            }

    # generate TSV file contents
    tsv_aa = ""
    tsv_3di = ""
    tsv_header = ""
    for i, id in enumerate(sequences_aa.keys()):
        tsv_aa += "{}\t{}\n".format(str(i + 1), sequences_aa[id])
        tsv_3di += "{}\t{}\n".format(str(i + 1), sequences_3di[id])
        tsv_header += "{}\t{}\n".format(str(i + 1), id)

    #### write temp tsv files

    # write TSV files
    temp_aa_tsv: Path = Path(foldseek_db_path) / "aa.tsv"
    temp_3di_tsv: Path = Path(foldseek_db_path) / "3di.tsv"
    temp_header_tsv: Path = Path(foldseek_db_path) / "header.tsv"
    with open(temp_aa_tsv, "w") as f:
        f.write(tsv_aa)
    with open(temp_3di_tsv, "w") as f:
        f.write(tsv_3di)
    with open(temp_header_tsv, "w") as f:
        f.write(tsv_header)

    # create foldseek db names

    short_db_name = f"{prefix}"
    aa_db_name: Path = Path(foldseek_db_path) / short_db_name
    tsv_db_name: Path = Path(foldseek_db_path) / f"{short_db_name}_ss"
    header_db_name: Path = Path(foldseek_db_path) / f"{short_db_name}_h"

    # create Foldseek database with foldseek tsv2db

    foldseek_tsv2db(temp_aa_tsv, aa_db_name, 0, logdir)
    foldseek_tsv2db(temp_3di_tsv, tsv_db_name, 0, logdir)
    foldseek_tsv2db(temp_header_tsv, header_db_name, 12, logdir)

    # clean up
    remove_file(temp_aa_tsv)
    remove_file(temp_3di_tsv)
    remove_file(temp_header_tsv)


def foldseek_tsv2db(
    in_tsv: Path, out_db_name: Path, db_type: int, logdir: Path
) -> None:
    """
    Convert a Foldseek TSV file to a Foldseek database.

    Args:
        in_tsv (Path): Path to the input TSV file.
        out_db_name (Path): Path for the output Foldseek database.
        db_type (int): Type of the output database.
        logdir (Path): Path to the directory where logs will be stored.

    Returns:
        None
    """
    foldseek_tsv2db = ExternalTool(
        tool="foldseek",
        input=f"",
        output=f"",
        params=f"tsv2db {in_tsv} {out_db_name}  --output-dbtype {str(db_type)} ",
        logdir=logdir,
    )

    ExternalTool.run_tool(foldseek_tsv2db)


def generate_foldseek_db_from_structures(
    fasta_aa: Path,
    foldseek_db_path: Path,
    structure_dir: Path,
    logdir: Path,
    prefix: str,
    proteins_flag: bool,
) -> None:
    """
    Generate Foldseek database from PDB files.

    Args:
        fasta_aa (Path): Path to the amino-acid FASTA file.
        foldseek_db_path (Path): Path to the directory where Foldseek database will be stored.
        structure_dir (Path): Path to the directory containing .pdb or .cif structure files.
        logdir (Path): Path to the directory where logs will be stored.
        prefix (str): Prefix for the Foldseek database.
        proteins_flag (bool): Flag - True if proteins-compare is run

    Returns:
        None
    """

    # read in amino-acid sequences
    sequences_aa = {}
    for record in SeqIO.parse(fasta_aa, "fasta"):
        sequences_aa[record.id] = str(record.seq)

    # lists all the pdb files

    structure_files = [
        file
        for file in os.listdir(structure_dir)
        if file.endswith(".pdb") or file.endswith(".cif")
    ]

    num_structures = len(structure_files)

    num_structures = 0

    # Checks that ID is in the pdbs

    no_structure_cds_ids = []

    for cds_id in sequences_aa.keys():

        matching_files = [
            file
            for file in structure_files
            if f"{cds_id}.pdb" == file or f"{cds_id}.cif" == file
        ]

        if len(matching_files) == 1:
            num_structures += 1

        # should neve happen but in case
        if len(matching_files) > 1:
            logger.warning(f"More than 1 structures found for {cds_id}")
            logger.warning("Taking the first one")
            num_structures += 1
        elif len(matching_files) == 0:
            logger.warning(f"No structure found for {cds_id}")
            logger.warning(f"{cds_id} will be ignored in annotation")
            no_structure_cds_ids.append(cds_id)

    if num_structures == 0:
        logger.error(
            f"No structures with matching CDS ids were found at all. Check the {structure_dir} directory"
        )

    # generate the db
    short_db_name = f"{prefix}"
    structure_db_name: Path = Path(foldseek_db_path) / short_db_name
    query_structure_dir = structure_dir


    foldseek_createdb_from_structures = ExternalTool(
        tool="foldseek",
        input=f"",
        output=f"",
        params=f"createdb {query_structure_dir} {structure_db_name} ",
        logdir=logdir,
    )

    ExternalTool.run_tool(foldseek_createdb_from_structures)


#### foldseek_gpu 


def create_foldseek_prostt5_gpu_db(
    fasta_aa: Path, foldseek_db_path: Path, db_dir: Path, logdir: Path
) -> None:
    """
    Convert a Foldseek DB with ProstT5 3Di predictions using Foldseek-GPU

    Args:
        fasta_aa (Path): Path to the amino-acid FASTA file.
        foldseek_db_path (Path): Path to the directory where Foldseek database will be stored.
        db_dir (Path): Path to the baktfold DB
        logdir (Path): Path to the directory where logs will be stored.
    Returns:
        None
    """

    prostt5_db_path = Path(db_dir) / "prostt5_weights"
    
    foldseek_createdb_prostt5 = ExternalTool(
        tool="foldseek",
        input=f"",
        output=f"",
        params=f"createdb {fasta_aa} {foldseek_db_path}  --prostt5-model {prostt5_db_path}  ",
        logdir=logdir,
    )

    ExternalTool.run_tool(foldseek_createdb_prostt5)
