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
from baktfold.features.run_foldseek import create_result_tsv, run_foldseek_search
from baktfold.io.handle_genbank import write_genbank
from baktfold.results.tophit import get_tophit
import baktfold.bakta.pstc as pstc

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
    foldseek_gpu: bool
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
        print('deal with this later should be trivial')



    ############
    # create foldseek db
    ############

    foldseek_query_db_path: Path = Path(output) / "foldseek_db"
    foldseek_query_db_path.mkdir(parents=True, exist_ok=True)

    if structures is True:
        logger.info("Creating a foldseek query database from structures.")
        logger.info("deal with this later also")

        generate_foldseek_db_from_structures(
            fasta_aa,
            foldseek_query_db_path,
            structure_dir,
            filtered_structures_path,
            logdir,
            prefix,
            filter_structures,
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
    result_tsv: Path = Path(output) / "foldseek_results.tsv"

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

    swissprot_df = get_tophit(result_tsv, structures)











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

    afdbclusters_df = get_tophit(result_tsv,structures)


    ####
    # looup
    ####
    # add the Swissprot and AFDB tophits to the json
    hypotheticals = pstc.parse(hypotheticals, swissprot_df, 'swissprot')
    hypotheticals = pstc.parse(hypotheticals, afdbclusters_df, 'afdb')

    # get the lookup descriptions for each of them
    hypotheticals = pstc.lookup(hypotheticals, Path(database))


    return hypotheticals


    

    # custom db output 

    #####
    # custom db
    #####

    # if custom_db:

    #     logger.info(f"Foldseek will also be run against your custom database {custom_db}")
    #     # make result and temp dirs
    #     result_db_custom: Path = Path(result_db_base) / "result_db_custom"
    #     result_tsv_custom: Path = Path(output) / "foldseek_results_custom.tsv"

    #     #try:
    #     run_foldseek_search(
    #     query_db,
    #     Path(custom_db),
    #     result_db_custom,
    #     temp_db,
    #     threads,
    #     logdir,
    #     evalue,
    #     sensitivity,
    #     max_seqs,
    #     ultra_sensitive,
    #     extra_foldseek_params,
    #     foldseek_gpu,
    #     structures,
    #     clustered_db=False # no custom db cluster searching
    # )

      
    
    #     # make result tsv
    #     create_result_tsv(query_db, Path(custom_db), result_db_custom,  result_tsv_custom, logdir, foldseek_gpu, structures, threads)

    #     tophit_custom_df = get_topcustom_hits(
    #     result_tsv_custom, structures, proteins_flag)



    #     #### merge 
    #     # left merge on the cds_id to get every cds/contig id (make it easier for downstream processing)

    #     if not proteins_flag: # if not proteins, need the contig_id
    #         all_cds_df = merged_df[['contig_id','cds_id']]
    #     else:
    #         all_cds_df = merged_df[['cds_id']]
    #     tophit_custom_df = all_cds_df.merge(tophit_custom_df, how="left", on='cds_id') # cds_id will always be unique

    #     # get the final column order required
    #     if proteins_flag: # no contig_id
    #         columns_order = ['cds_id'] + [col for col in tophit_custom_df.columns if col not in ['contig_id', 'cds_id']]
    #     else: # including structures, will have contig_id too in the merged_df
    #         columns_order = ['contig_id', 'cds_id'] + [col for col in tophit_custom_df.columns if col not in ['contig_id', 'cds_id']]
    #     tophit_custom_df = tophit_custom_df[columns_order]

    #     # get coverages
    #     tophit_custom_df = calculate_qcov_tcov(tophit_custom_df)
    #     custom_hits_path: Path = Path(output) / f"{prefix}_custom_database_hits.tsv"
    #     tophit_custom_df.to_csv(custom_hits_path, index=False, sep="\t")
        
    #     # except:
    #     #     logger.error(f"Foldseek failed to run against your custom database {custom_db}. Please check that it is formatted correctly as a Foldseek database")



    return True
