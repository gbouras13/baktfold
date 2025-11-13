#!/usr/bin/env python3

import shutil
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import numpy as np
import pandas as pd
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord
from loguru import logger

from baktfold.features.create_foldseek_db import (
    generate_foldseek_db_from_aa_3di, generate_foldseek_db_from_structures)
from baktfold.features.run_foldseek import create_result_tsv, run_foldseek_search, summarise_hits
from baktfold.io.handle_genbank import write_genbank
import baktfold.io.io as io
from baktfold.results.tophit import get_tophit
import baktfold.bakta.pstc as pstc
from baktfold.utils.util import remove_file

def subcommand_compare(
    hypotheticals: Dict,
    output: Path,
    threads: int,
    evalue: float,
    sensitivity: float,
    database: Path,
    prefix: str,
    predictions_dir: Optional[Path],
    structures: bool,
    structure_dir: Optional[Path],
    logdir: Path,
    proteins_flag: bool,
    max_seqs: int,
    ultra_sensitive: bool,
    extra_foldseek_params: str,
    custom_db: str,
    foldseek_gpu: bool,
    custom_annotations: Optional[Path]
) -> bool:
    """
    Compare 3Di or PDB structures to the baktfold DB

    Parameters:
        hypotheticals (Dict):  hypothetical features dictionary
        output (Path): Path to the output directory.
        threads (int): Number of threads to use.
        evalue (float): E-value threshold.
        card_vfdb_evalue (float): E-value threshold for CARD and VFDB databases.
        sensitivity (float): Sensitivity threshold.
        database (Path): Path to the reference database.
        prefix (str): Prefix for output files.
        predictions_dir (Optional[Path]): Path to the directory containing predictions.
        structures (bool): Flag indicating whether structures files are used.
        structure_dir (Optional[Path]): Path to the directory containing structures (.pdb or .cif) files.
        logdir (Path): Path to the directory for log files.
        proteins_flag (bool): Flag indicating whether proteins are used.
        max_seqs (int): Maximum results per query sequence allowed to pass the prefilter for foldseek.
        ultra_sensitive (bool): Whether to skip foldseek prefilter for maximum sensitivity
        extra_foldseek_params (str): Extra foldseek search parameters
        custom_db (str): Custom foldseek database
        foldseek_gpu (bool): Use Foldseek-GPU acceleration and ungappedprefilter
        custom_annotations (Optional[Path]): Path to the tsv containing the custom_db annotations, 2 columns 

    Returns:
        bool: True if sub-databases are created successfully, False otherwise.
    """


    # input predictions or structures
    if structures is False:
        # prostT5
        fasta_aa_input: Path = Path(predictions_dir) / f"{prefix}_aa.fasta"
        fasta_3di_input: Path = Path(predictions_dir) / f"{prefix}_3di.fasta"

    fasta_aa: Path = Path(output) / f"{prefix}_aa.fasta"
    fasta_3di: Path = Path(output) / f"{prefix}_3di.fasta"

    ## copy the AA and 3Di from predictions directory 
    # if structures is false and baktfold compare is the command
    # Otherwise it will just copy itself

    if structures is False:
        if fasta_3di_input.exists():
            logger.info(
                f"Checked that the 3Di CDS file {fasta_3di_input} exists from baktfold predict"
            )
            if fasta_3di.exists() is False:
                shutil.copyfile(fasta_3di_input, fasta_3di)
        else:
            logger.error(
                f"The 3Di CDS file {fasta_3di_input} does not exist. Please run baktfold predict and/or check the prediction directory {predictions_dir}"
            )
        # copy the aa to file
        if fasta_aa_input.exists():
            logger.info(
                f"Checked that the AA CDS file {fasta_aa_input} exists from baktfold predict."
            )
            if fasta_aa.exists() is False:
                shutil.copyfile(fasta_aa_input, fasta_aa)
        else:
            logger.error(
                f"The AA CDS file {fasta_aa_input} does not exist. Please run baktfold predict and/or check the prediction directory {predictions_dir}"
                )

    ## write the AAs to file if structures is true because can't just copy from prediction_dir
    else:
        ## write the CDS to file
        logger.info(f"Writing the AAs to file {fasta_aa}.")

        with open(fasta_aa, "w+") as out_f:
            for entry in hypotheticals:
                header = f">{entry['locus']}\n"
                seq = f"{entry['aa']}\n"
                out_f.write(header)
                out_f.write(seq)


    ############
    # create foldseek db
    ############

    foldseek_query_db_path: Path = Path(output) / "foldseek_db"
    foldseek_query_db_path.mkdir(parents=True, exist_ok=True)

    if structures is True:
        logger.info("Creating a foldseek query database from structures.")

        generate_foldseek_db_from_structures(
            fasta_aa,
            foldseek_query_db_path,
            structure_dir,
            logdir,
            prefix,
            proteins_flag,
        )
    else:
        generate_foldseek_db_from_aa_3di(
            fasta_aa, fasta_3di, foldseek_query_db_path, logdir, prefix
        )

    short_db_name = prefix

    # db search 

    database_name = "swissprot"

    if short_db_name == database_name:
        logger.error(
            f"Please choose a different -p {prefix} as this conflicts with the {database_name}"
        )

    #####
    # foldseek search
    #####

    query_db: Path = Path(foldseek_query_db_path) / short_db_name
    target_db: Path = Path(database) / database_name

    # make result and temp dirs
    result_db_base: Path = Path(output) / "result_db"
    result_db_base.mkdir(parents=True, exist_ok=True)
    result_db: Path = Path(result_db_base) / "result_db"

    temp_db: Path = Path(output) / "temp_db"
    temp_db.mkdir(parents=True, exist_ok=True)

    # make result tsv
    result_tsv: Path = Path(output) / "foldseek_results_swissprot.tsv"

    # run foldseek search
    run_foldseek_search(
        query_db,
        target_db,
        result_db,
        temp_db,
        threads,
        logdir,
        evalue,
        sensitivity,
        max_seqs,
        ultra_sensitive,
        extra_foldseek_params,
        foldseek_gpu,
        structures
    )

       
    create_result_tsv(query_db, target_db, result_db, result_tsv, logdir, foldseek_gpu, structures, threads)

    swissprot_df = get_tophit(result_tsv, structures, cath=False)




    #####
    # foldseek search AFDB Clusters
    #####




    database_name = "AFDBClusters"

    if short_db_name == database_name:
        logger.error(
            f"Please choose a different -p {prefix} as this conflicts with the {database_name}"
        )

    query_db: Path = Path(foldseek_query_db_path) / short_db_name
    target_db: Path = Path(database) / database_name

    # make result and temp dirs
    result_db_base: Path = Path(output) / "result_db"
    result_db_base.mkdir(parents=True, exist_ok=True)
    result_db: Path = Path(result_db_base) / "result_afdb_db"

    temp_db: Path = Path(output) / "temp_db"
    temp_db.mkdir(parents=True, exist_ok=True)

    # make result tsv
    result_tsv: Path = Path(output) / "foldseek_results_afdb_clusters.tsv"

    # run foldseek search
    run_foldseek_search(
        query_db,
        target_db,
        result_db,
        temp_db,
        threads,
        logdir,
        evalue,
        sensitivity,
        max_seqs,
        ultra_sensitive,
        extra_foldseek_params,
        foldseek_gpu,
        structures
    )

       
    create_result_tsv(query_db, target_db, result_db, result_tsv, logdir, foldseek_gpu, structures, threads)

    afdbclusters_df = get_tophit(result_tsv,structures, cath=False)

    #####
    # foldseek search pdb
    #####


    database_name = "pdb"

    if short_db_name == database_name:
        logger.error(
            f"Please choose a different -p {prefix} as this conflicts with the {database_name}"
        )

    query_db: Path = Path(foldseek_query_db_path) / short_db_name
    target_db: Path = Path(database) / database_name

    # make result and temp dirs
    result_db_base: Path = Path(output) / "result_db"
    result_db_base.mkdir(parents=True, exist_ok=True)
    result_db: Path = Path(result_db_base) / "result_pdb_db"

    temp_db: Path = Path(output) / "temp_db"
    temp_db.mkdir(parents=True, exist_ok=True)

    # make result tsv
    result_tsv: Path = Path(output) / "foldseek_results_pdb.tsv"

    # run foldseek search
    run_foldseek_search(
        query_db,
        target_db,
        result_db,
        temp_db,
        threads,
        logdir,
        evalue,
        sensitivity,
        max_seqs,
        ultra_sensitive,
        extra_foldseek_params,
        foldseek_gpu,
        structures
    )

       
    create_result_tsv(query_db, target_db, result_db, result_tsv, logdir, foldseek_gpu, structures, threads)

    pdb_df = get_tophit(result_tsv,structures, cath=False)


    #####
    # foldseek search cath
    #####


    database_name = "cath"

    if short_db_name == database_name:
        logger.error(
            f"Please choose a different -p {prefix} as this conflicts with the {database_name}"
        )

    query_db: Path = Path(foldseek_query_db_path) / short_db_name
    target_db: Path = Path(database) / database_name

    # make result and temp dirs
    result_db_base: Path = Path(output) / "result_db"
    result_db_base.mkdir(parents=True, exist_ok=True)
    result_db: Path = Path(result_db_base) / "result_cath_db"
    result_db_greedy_best_hits: Path = Path(result_db_base) / "result_cath_db_greedy_best_hits"

    temp_db: Path = Path(output) / "temp_db"
    temp_db.mkdir(parents=True, exist_ok=True)

    # make result tsv
    result_tsv: Path = Path(output) / "foldseek_results_cath.tsv"
    result_greedy_tsv: Path = Path(output) /  "foldseek_results_cath_greedy_tophit"

    # run foldseek search
    run_foldseek_search(
        query_db,
        target_db,
        result_db,
        temp_db,
        threads,
        logdir,
        evalue,
        sensitivity,
        max_seqs,
        ultra_sensitive,
        extra_foldseek_params,
        foldseek_gpu,
        structures
    )

    # this keeps the greedy best hits for cath
    # we actually don't keep the single tophit - multidomain/fold proteins should have multiple non-overlapping CATH hits
    # this is equivalent to using --greedy-best-hits with foldseek easy-search
    summarise_hits(result_db, result_db_greedy_best_hits, logdir, threads)

    # saves all CATH hits first
    create_result_tsv(query_db, target_db, result_db, result_tsv, logdir, foldseek_gpu, structures, threads)
    # save greedy CATH tophits
    create_result_tsv(query_db, target_db, result_db_greedy_best_hits, result_greedy_tsv, logdir, foldseek_gpu, structures, threads)

    # this just reads it in with appropriate headers
    cath_df = get_tophit(result_greedy_tsv, structures, cath=True)

    # write tophits
    swissprot_tophit_path: Path = Path(output) / "baktfold_swissprot_tophit.tsv"
    io.write_foldseek_tophit(swissprot_df, swissprot_tophit_path)

    afdb_tophit_path: Path = Path(output) / "baktfold_afdbclusters_tophit.tsv"
    io.write_foldseek_tophit(afdbclusters_df, afdb_tophit_path)

    pdb_tophit_path: Path = Path(output) / "baktfold_pdb_tophit.tsv"
    io.write_foldseek_tophit(pdb_df, pdb_tophit_path)

    cath_tophit_path: Path = Path(output) / "baktfold_cath_tophit.tsv"
    io.write_foldseek_tophit(cath_df, cath_tophit_path)
    # remove result_greedy_tsv (identical to tophit, will make it confusing)
    remove_file(result_greedy_tsv) 


    # custom db output 

    #####
    # custom db
    #####


    if custom_db:

        try:

            logger.info(f"Foldseek will also be run against your custom database {custom_db}")
            # make result and temp dirs
            result_db_custom: Path = Path(result_db_base) / "result_db_custom"
            result_tsv_custom: Path = Path(output) / "foldseek_results_custom.tsv"

            run_foldseek_search(
            query_db,
            Path(custom_db),
            result_db_custom,
            temp_db,
            threads,
            logdir,
            evalue,
            sensitivity,
            max_seqs,
            ultra_sensitive,
            extra_foldseek_params,
            foldseek_gpu,
            structures
        )

            create_result_tsv(query_db, Path(custom_db),
                result_db_custom,
                result_tsv_custom, logdir, foldseek_gpu, structures, threads)

            custom_df = get_tophit(result_tsv_custom,structures, cath=False)

            custom_db_tophit_path: Path = Path(output) / "baktfold_custom_db_tophit.tsv"
            io.write_foldseek_tophit(custom_df, custom_db_tophit_path)
        
        except:
            logger.error(f"Foldseek failed to run against your custom database {custom_db}. Please check that it is formatted correctly as a Foldseek database")


    ####
    # lookup
    ####

    if proteins_flag: # baktfold proteins

        # note aas passed as hypotheticals to the overall function - so in and out as aas

        aas = pstc.parse(hypotheticals, swissprot_df, 'swissprot')
        aas = pstc.parse(aas, afdbclusters_df, 'afdb')
        aas = pstc.parse(aas, pdb_df, 'pdb')
        aas = pstc.parse(aas, cath_df, 'cath')
        if custom_db:
            aas = pstc.parse(aas, custom_df, 'custom_db')

        # get the lookup descriptions for each of them
        # this requires the DB

        #aas = pstc.lookup(aas, Path(database), custom_annotations)
        aas = pstc.lookup_sql(aas, Path(database), threads)
        # add the custom annotations if it is provided
        if custom_annotations:
            aas = pstc.lookup_custom(aas, Path(database), custom_annotations)

        return aas

    else: # baktfold run

        # add the Swissprot and AFDB and PDB tophits to the json
        hypotheticals = pstc.parse(hypotheticals, swissprot_df, 'swissprot')
        hypotheticals = pstc.parse(hypotheticals, afdbclusters_df, 'afdb')
        hypotheticals = pstc.parse(hypotheticals, pdb_df, 'pdb')
        hypotheticals = pstc.parse(hypotheticals, cath_df, 'cath')
        if custom_db:
            hypotheticals = pstc.parse(hypotheticals, custom_df, 'custom_db')

        # get the lookup descriptions for each of them
        # hypotheticals = pstc.lookup(hypotheticals, Path(database), custom_annotations)
        hypotheticals = pstc.lookup_sql(hypotheticals, Path(database), threads)
        if custom_annotations:
            hypotheticals = pstc.lookup_custom(hypotheticals, Path(database), custom_annotations)

        return hypotheticals


    





