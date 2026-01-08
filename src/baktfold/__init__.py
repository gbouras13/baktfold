#!/usr/bin/env python3

__version__ = '0.1.0'

import gzip
from pathlib import Path

from datetime import datetime
import click
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation, SeqFeature
from loguru import logger

import baktfold.bakta.constants as bc
import baktfold.bakta.annotation as anno
from baktfold.io.json_in import parse_json_input
from baktfold.io.fasta_in import parse_protein_input
from baktfold.databases.db import install_database, validate_db
from baktfold.features.create_foldseek_db import generate_foldseek_db_from_aa_3di
from baktfold.features.predict_3Di import get_T5_model
from baktfold.io.handle_genbank import open_protein_fasta_file
from baktfold.subcommands.compare import subcommand_compare
from baktfold.subcommands.predict import subcommand_predict
from baktfold.utils.constants import DB_DIR, CNN_DIR
from baktfold.utils.util import (begin_baktfold, clean_up_temporary_files, end_baktfold,
                              get_version, print_citation, sort_euk_feature_key)
from baktfold.utils.validation import (check_dependencies, instantiate_dirs,validate_outfile,
                                    check_genbank_and_prokka)

from baktfold.io.prokka_gbk_to_json import prokka_gbk_to_json
from baktfold.io.eukaryotic_to_json import eukaryotic_gbk_to_json
import baktfold.bakta.config as cfg
import baktfold.io.io as io
from baktfold.features.autotune import run_autotune
from importlib.resources import files

log_fmt = (
    "[<green>{time:YYYY-MM-DD HH:mm:ss}</green>] <level>{level: <8}</level> | "
    "<level>{message}</level>"
)

"""
common options
"""


def common_options(func):
    """Common command line args
    Define common command line args here, and include them with the @common_options decorator below.
    """
    options = [
        click.option(
            "-o",
            "--output",
            default="output_baktfold",
            show_default=True,
            type=click.Path(),
            help="Output directory ",
        ),
        click.option(
            "-t",
            "--threads",
            help="Number of threads",
            default=1,
            type=int,
            show_default=True,
        ),
        click.option(
            "-p",
            "--prefix",
            default="baktfold",
            help="Prefix for output files",
            type=str,
            show_default=True,
        ),
        click.option(
            "-d",
            "--database",
            type=str,
            default=None,
            help="Specific path to installed baktfold database",
        ),
        click.option(
            "-f",
            "--force",
            is_flag=True,
            help="Force overwrites the output directory",
        ),
    ]
    for option in reversed(options):
        func = option(func)
    return func


"""
predict only options
"""


def predict_options(func):
    """predict command line args"""
    options = [
        click.option(
            "--autotune",
            is_flag=True,
            help="Run autotuning to detect and automatically use best batch size for your hardware. Recommended only if you have a large dataset (e.g. thousands of proteins), or else autotuning will add rather than save runtime.",
        ),
        click.option(
            "--batch-size",
            default=1,
            help="batch size for ProstT5. 1 is usually fastest.",
            show_default=True,
        ),
        click.option(
            "--cpu",
            is_flag=True,
            help="Use cpus only.",
        ),
        click.option(
            "--omit-probs",
            is_flag=True,
            help="Do not output per residue 3Di probabilities from ProstT5. Mean per protein 3Di probabilities will always be output.",
        ),
        click.option(
            "--save-per-residue-embeddings",
            is_flag=True,
            help="Save the ProstT5 embeddings per resuide in a h5 file ",
        ),
        click.option(
            "--save-per-protein-embeddings",
            is_flag=True,
            help="Save the ProstT5 embeddings as means per protein in a h5 file",
        ),
        click.option(
            "--mask-threshold",
            default=25,
            help="Masks 3Di residues below this value of ProstT5 confidence for Foldseek searches",
            type=float,
            show_default=True,
        )
    ]
    for option in reversed(options):
        func = option(func)
    return func


"""
compare only options
"""


def compare_options(func):
    """compare command line args"""
    options = [
        click.option(
            "-e",
            "--evalue",
            default="1e-3",
            type=float,
            help="Evalue threshold for Foldseek",
            show_default=True,
        ),
        click.option(
            "-s",
            "--sensitivity",
            default="9.5",
            help="Sensitivity parameter for foldseek",
            type=float,
            show_default=True,
        ),
        click.option(
            "--keep-tmp-files",
            is_flag=True,
            help="Keep temporary intermediate files, particularly the large foldseek_results.tsv of all Foldseek hits",
        ),
        click.option(
            "--max-seqs",
            type=int,
            default=1000,
            show_default=True,
            help="Maximum results per query sequence allowed to pass the prefilter. You may want to reduce this to save disk space for enormous datasets",
        ),
        click.option(
            "--ultra-sensitive",
            is_flag=True,
            help="Runs baktfold with maximum sensitivity by skipping Foldseek prefilter. Not recommended for large datasets.",
        ),
        click.option(
            "--extra-foldseek-params",
            type=str,
            help="Extra foldseek search params"
        ),
        click.option(
            "--custom-db",
            type=str,
            help="Path to custom database"
        ),
        click.option(
            "--foldseek-gpu",
            is_flag=True,
            help="Use this to enable compatibility with Foldseek-GPU search acceleration",
        ),
        click.option(
            "--custom-annotations",
            type=click.Path(),
            help="Custom Foldseek DB annotations, 2 column tsv. Column 1 matches the Foldseek headers, column 2 is the description.",
        ),
        click.option(
            "--euk",
            is_flag=True,
            help="Eukaryotic input genome.",
        )
    ]
    for option in reversed(options):
        func = option(func)
    return func

"""
bakta input options

Only for baktfold run predict and compate
"""

def bakta_options(func):
    """compare command line args"""
    options = [
        click.option(
            "-a",
            "--all-proteins",
            is_flag=True,
            help="annotate all proteins (not just hypotheticals)",
        ),
    ]
    for option in reversed(options):
        func = option(func)
    return func

        


@click.group()
@click.help_option("--help", "-h")
@click.version_option(get_version(), "--version", "-V")
def main_cli():
    1 + 1


"""
run command
"""


@main_cli.command()
@click.help_option("--help", "-h")
@click.version_option(get_version(), "--version", "-V")
@click.pass_context
@click.option(
    "-i",
    "--input",
    help="Path to input file in Bakta Genbank format or Bakta JSON format",
    type=click.Path(),
    required=True,
)
@common_options
@predict_options
@compare_options
@bakta_options
def run(
    ctx,
    input,
    output,
    threads,
    prefix,
    evalue,
    force,
    database,
    autotune,
    batch_size,
    sensitivity,
    cpu,
    omit_probs,
    keep_tmp_files,
    max_seqs,
    save_per_residue_embeddings,
    save_per_protein_embeddings,
    ultra_sensitive,
    mask_threshold,
    extra_foldseek_params,
    custom_db,
    custom_annotations,
    foldseek_gpu,
    all_proteins,
    euk,
    **kwargs,
):
    """baktfold predict then comapare all in one - GPU recommended"""

    # validates the directory  (need to before I start baktfold or else no log file is written)
    instantiate_dirs(output, force)

    output: Path = Path(output)
    logdir: Path = Path(output) / "logs"

    params = {
        "--input": input,
        "--output": output,
        "--threads": threads,
        "--force": force,
        "--prefix": prefix,
        "--evalue": evalue,
        "--database": database,
        "--autotune": autotune,
        "--batch-size": batch_size,
        "--sensitivity": sensitivity,
        "--keep-tmp-files": keep_tmp_files,
        "--cpu": cpu,
        "--omit_probs": omit_probs,
        "--max-seqs": max_seqs,
        "--save-per-residue-embeddings": save_per_residue_embeddings,
        "--save-per-protein-embeddings": save_per_protein_embeddings,
        "--ultra-sensitive": ultra_sensitive,
        "--mask-threshold": mask_threshold,
        "--extra-foldseek-params": extra_foldseek_params,
        "--custom-db": custom_db,
        "--custom-annotations": custom_annotations,
        "--foldseek-gpu": foldseek_gpu,
        "--all-proteins": all_proteins,
        "--euk": euk
    }

    # initial logging etc
    start_time = begin_baktfold(params, "run")

    # check foldseek is installed
    check_dependencies()

    # check the database is installed and return it
    database = validate_db(database, DB_DIR, foldseek_gpu)

    ###
    # parse the json output and save hypotheticals as AA FASTA
    ###

    fasta_aa: Path = Path(output) / f"{prefix}_aa.fasta"
    data, features, has_duplicate_locus = parse_json_input(input, fasta_aa, all_proteins)

    ###
    # split features in hypotheticals and non hypotheticals
    ###


    if all_proteins:
        hypotheticals = [feat for feat in features if feat['type'] == bc.FEATURE_CDS ]
        
        non_hypothetical_features = [
        feat for feat in features
        if (feat['type'] != bc.FEATURE_CDS) 
    ]
    else:

        hypotheticals = [feat for feat in features if feat['type'] == bc.FEATURE_CDS and 'hypothetical' in feat]
        non_hypothetical_features = [
        feat for feat in features
        if (feat['type'] != bc.FEATURE_CDS) or 
        (feat['type'] == bc.FEATURE_CDS and 'hypothetical' not in feat)
    ]

    # put the CDS AA in a simple dictionary for ProstT5 code
    cds_dict = {}
    for feat in hypotheticals:
        if has_duplicate_locus:
            cds_dict[feat['id']] = feat['aa']
        else:
            cds_dict[feat['locus']] = feat['aa']


    # add a function to add 3Di to cds_dict

    # baktfold predict
    model_dir = database
    model_name = "Rostlab/ProstT5_fp16"
    checkpoint_path = Path(CNN_DIR) / "cnn_chkpnt" / "model.pt"


    if autotune:

        input_path = files("baktfold.features.autotune_data").joinpath("swissprot_5000.fasta.gz")

        step = 20
        min_batch = 1
        max_batch = 301
        sample_seqs = 602

        batch_size = run_autotune(
            input_path,
            model_dir,
            model_name,
            cpu,
            threads,
            step, 
            min_batch,
            max_batch, 
            sample_seqs)


    # hypotheticals is input to the function as it updates the 3Di feature

    hypotheticals = subcommand_predict(
        hypotheticals,
        cds_dict,
        output,
        prefix,
        cpu,
        omit_probs,
        model_dir,
        model_name,
        checkpoint_path,
        batch_size,
        save_per_residue_embeddings=save_per_residue_embeddings,
        save_per_protein_embeddings=save_per_protein_embeddings,
        threads=threads,
        mask_threshold=mask_threshold,
        has_duplicate_locus=has_duplicate_locus
    )

    # baktfold compare

    # predictions_dir is output as this will be where it lives for 'run'

    hypotheticals = subcommand_compare(
        hypotheticals,
        output,
        threads,
        evalue,
        sensitivity,
        database,
        prefix,
        predictions_dir=output,
        structures=False,
        structure_dir=None,
        logdir=logdir,
        proteins_flag=False,
        max_seqs=max_seqs,
        ultra_sensitive=ultra_sensitive,
        extra_foldseek_params=extra_foldseek_params,
        custom_db=custom_db,
        foldseek_gpu=foldseek_gpu,
        custom_annotations=custom_annotations,
        has_duplicate_locus=has_duplicate_locus
    )

    #####
    # update the hypotheticals 
    #####

    for cds in hypotheticals:
        anno.combine_annotation(cds)  # add on PSTC annotations and mark hypotheticals

    # recombine updated and existing features
    combined_features = non_hypothetical_features + hypotheticals  # recombine
    
    # Sort by ascending 'id'
    combined_features_sorted = sorted(combined_features, key=lambda x: x.get('id', ''))

    # put back in dictionary
    data['features'] = combined_features_sorted

    # map features by sequence for io
    features_by_sequence = {seq['id']: [] for seq in data['sequences']}

    for feature in data['features']:
        if 'discarded' not in feature:
            seq_features = features_by_sequence.get(feature['sequence'])
            if seq_features is not None:
                seq_features.append(feature)

    # flatten sorted features
    features = []
    for seq in data['sequences']:
        seq_features = features_by_sequence[seq['id']]
        if euk: # ensure gene -> mRNA -> CDS ordering for each locus tag
            seq_features.sort(key=sort_euk_feature_key)
        else:
            seq_features.sort(key=lambda k: k['start'])  # sort features by start position
        features.extend(seq_features)

    # overwrite feature list with sorted features
    data['features'] = features

    #####
    # don't include this for now as no gene symbols
    #####

    # logger.info('improve annotations...')
    # genes_with_improved_symbols = anno.select_gene_symbols([feature for feature in features if feature['type'] in [bc.FEATURE_CDS, bc.FEATURE_SORF]])
    # print(f'\trevised gene symbols: {len(genes_with_improved_symbols)}')

    ####
    # bakta output module
    ####



    logger.info('writing baktfold outputs')
    io.write_bakta_outputs(data, features, features_by_sequence, output, prefix, custom_db, euk, has_duplicate_locus)

    # cleanup the temp files
    if not keep_tmp_files:
        clean_up_temporary_files(output, prefix)

   
    # end baktfold
    end_baktfold(start_time, "run")




"""
proteins command
"""


@main_cli.command()
@click.help_option("--help", "-h")
@click.version_option(get_version(), "--version", "-V")
@click.pass_context
@click.option(
    "-i",
    "--input",
    help="Path to input file in amino acid FASTA format",
    type=click.Path(),
    required=True,
)
@common_options
@predict_options
@compare_options
def proteins(
    ctx,
    input,
    output,
    threads,
    prefix,
    evalue,
    force,
    database,
    autotune,
    batch_size,
    sensitivity,
    cpu,
    omit_probs,
    keep_tmp_files,
    max_seqs,
    save_per_residue_embeddings,
    save_per_protein_embeddings,
    ultra_sensitive,
    mask_threshold,
    extra_foldseek_params,
    custom_db,
    foldseek_gpu,
    custom_annotations,
    **kwargs,
):
    """baktfold proteins-predict then comapare all in one - GPU recommended"""

    # validates the directory  (need to before I start baktfold or else no log file is written)
    instantiate_dirs(output, force)

    output: Path = Path(output)
    logdir: Path = Path(output) / "logs"

    params = {
        "--input": input,
        "--output": output,
        "--threads": threads,
        "--force": force,
        "--prefix": prefix,
        "--evalue": evalue,
        "--database": database,
        "--autotune": autotune,
        "--batch_size": batch_size,
        "--sensitivity": sensitivity,
        "--keep-tmp-files": keep_tmp_files,
        "--cpu": cpu,
        "--omit-probs": omit_probs,
        "--max-seqs": max_seqs,
        "--save-per-residue-embeddings": save_per_residue_embeddings,
        "--save-per-protein-embeddings": save_per_protein_embeddings,
        "--ultra-sensitive": ultra_sensitive,
        "--mask-threshold": mask_threshold,
        "--extra-foldseek_params": extra_foldseek_params,
        "--custom-db": custom_db,
        "--foldseek-gpu": foldseek_gpu,
        "--custom-annotations": custom_annotations
    }

    # initial logging etc
    start_time = begin_baktfold(params, "proteins")

    # check foldseek is installed
    check_dependencies()

    # check the database is installed and return it
    database = validate_db(database, DB_DIR, foldseek_gpu)


    ###
    # parse the input and save hypotheticals as AA FASTA
    ###

    fasta_aa: Path = Path(output) / f"{prefix}_aa.fasta"
    # puts the CDS AA in a simple dictionary for ProstT5 code
    aas = parse_protein_input(input, fasta_aa)

    # put the CDS AA in a simple dictionary for ProstT5 code
    cds_dict = {}
    for feat in aas:
        cds_dict[feat['locus']] = feat['aa']

    # baktfold predict
    model_dir = database
    model_name = "Rostlab/ProstT5_fp16"
    checkpoint_path = Path(CNN_DIR) / "cnn_chkpnt" / "model.pt"


    if autotune:

        input_path = files("baktfold.features.autotune_data").joinpath("swissprot_5000.fasta.gz")

        step = 20
        min_batch = 1
        max_batch = 301
        sample_seqs = 602

        batch_size = run_autotune(
            input_path,
            model_dir,
            model_name,
            cpu,
            threads,
            step, 
            min_batch,
            max_batch, 
            sample_seqs)

    aas = subcommand_predict(
        aas,
        cds_dict,
        output,
        prefix,
        cpu,
        omit_probs,
        model_dir,
        model_name,
        checkpoint_path,
        batch_size,
        save_per_residue_embeddings=save_per_residue_embeddings,
        save_per_protein_embeddings=save_per_protein_embeddings,
        threads=threads,
        mask_threshold=mask_threshold,
        has_duplicate_locus=False
    )

    # baktfold compare


    # predictions_dir is output as this will be where it lives for 'run'

    aas = subcommand_compare(
        aas, # this is dummy, no operations happen here
        output,
        threads,
        evalue,
        sensitivity,
        database,
        prefix,
        predictions_dir=output,
        structures=False,
        structure_dir=None,
        logdir=logdir,
        proteins_flag=True,
        max_seqs=max_seqs,
        ultra_sensitive=ultra_sensitive,
        extra_foldseek_params=extra_foldseek_params,
        custom_db=custom_db,
        foldseek_gpu=foldseek_gpu,
        custom_annotations=custom_annotations,
        has_duplicate_locus=False
    )

    #####
    # update the hypotheticals 
    #####

    for aa in aas:
        anno.combine_annotation(aa)  # add on PSTC annotations and mark hypotheticals


    ####
    # bakta output module
    ####
    logger.info('writing baktfold outputs')

    cfg.run_end = datetime.now()
    run_duration = (cfg.run_end - cfg.run_start).total_seconds()

    for aa in aas:  # reset mock attributes
        aa['start'] = -1
        aa['stop'] = -1

    ############################################################################
    # Write output files
    # - write comprehensive annotation results as JSON
    # - write optional output files in TSV, FAA formats
    # - remove temp directory
    ############################################################################
    
    io.write_bakta_proteins_outputs(aas, output, prefix, custom_db)

    # cleanup the temp files
    if not keep_tmp_files:
        clean_up_temporary_files(output, prefix)

    # end baktfold
    end_baktfold(start_time, "proteins")





"""
predict command
Uses ProstT5 to predict 3Di sequences from bakta json input 
"""


@main_cli.command()
@click.help_option("--help", "-h")
@click.version_option(get_version(), "--version", "-V")
@click.pass_context
@click.option(
    "-i",
    "--input",
    help="Path to input file in Genbank format or nucleotide FASTA format",
    type=click.Path(),
    required=True,
)
@common_options
@predict_options
@bakta_options



def predict(
    ctx,
    input,
    output,
    threads,
    prefix,
    force,
    database,
    autotune,
    batch_size,
    cpu,
    omit_probs,
    save_per_residue_embeddings,
    save_per_protein_embeddings,
    mask_threshold,
    all_proteins,
    **kwargs,
):

    """Uses ProstT5 to predict 3Di tokens - GPU recommended"""


    # validates the directory  (need to before I start baktfold or else no log file is written)
    instantiate_dirs(output, force)

    output: Path = Path(output)
    logdir: Path = Path(output) / "logs"

    params = {
        "--input": input,
        "--output": output,
        "--threads": threads,
        "--force": force,
        "--prefix": prefix,
        "--database": database,
        "--autotune": autotune,
        "--batch-size": batch_size,
        "--cpu": cpu,
        "--omit-probs": omit_probs,
        "--save-per-residue-embeddings": save_per_residue_embeddings,
        "--save-per-protein-embeddings": save_per_protein_embeddings,
        "--mask-threshold": mask_threshold,
        "--all-proteins": all_proteins,

    }

    # initial logging etc
    start_time = begin_baktfold(params, "predict")

    # check foldseek is installed
    check_dependencies()

    # check the database is installed and return it
    database = validate_db(database, DB_DIR, foldseek_gpu=False) # dont need the foldseek gpu here as not an inout option


    ###
    # parse the json output and save hypotheticals as AA FASTA
    ###

    fasta_aa: Path = Path(output) / f"{prefix}_aa.fasta"
    data, features, has_duplicate_locus = parse_json_input(input, fasta_aa, all_proteins)

    ###
    # split features in hypotheticals and non hypotheticals
    ###

    if all_proteins:
        hypotheticals = [
            feat for feat in features
            if feat['type'] == bc.FEATURE_CDS
        ]
    else:
        hypotheticals = [
            feat for feat in features
            if feat['type'] == bc.FEATURE_CDS and 'hypothetical' in feat
        ]

    non_hypothetical_features = [
        feat for feat in features
        if feat not in hypotheticals
    ]



    # put the CDS AA in a simple dictionary for ProstT5 code
    cds_dict = {}
    if has_duplicate_locus:
        for feat in hypotheticals:
            cds_dict[feat['id']] = feat['aa']
    
    else:
        for feat in hypotheticals:
            cds_dict[feat['locus']] = feat['aa']


    # add a function to add 3Di to cds_dict

    # baktfold predict
    model_dir = database
    model_name = "Rostlab/ProstT5_fp16"
    checkpoint_path = Path(CNN_DIR) / "cnn_chkpnt" / "model.pt"

    if autotune:

        input_path = files("baktfold.features.autotune_data").joinpath("swissprot_5000.fasta.gz")

        step = 20
        min_batch = 1
        max_batch = 301
        sample_seqs = 602

        batch_size = run_autotune(
            input_path,
            model_dir,
            model_name,
            cpu,
            threads,
            step, 
            min_batch,
            max_batch, 
            sample_seqs)

    hypotheticals = subcommand_predict(
        hypotheticals,
        cds_dict,
        output,
        prefix,
        cpu,
        omit_probs,
        model_dir,
        model_name,
        checkpoint_path,
        batch_size,
        save_per_residue_embeddings=save_per_residue_embeddings,
        save_per_protein_embeddings=save_per_protein_embeddings,
        threads=threads,
        mask_threshold=mask_threshold,
        has_duplicate_locus=has_duplicate_locus
    )

    # end baktfold
    end_baktfold(start_time, "predict")


"""
compare command

runs Foldseek using either 1) output of baktfold predict or 2) user defined protein structures and generates compliant outputs

"""


@main_cli.command()
@click.help_option("--help", "-h")
@click.version_option(get_version(), "--version", "-V")
@click.pass_context
@click.option(
    "-i",
    "--input",
    help="Path to input file in Genbank format or nucleotide FASTA format",
    type=click.Path(),
    required=True,
)
@click.option(
    "--predictions-dir",
    help="Path to output directory from baktfold predict",
    type=click.Path(),
    default=None,
)
@click.option(
    "--structure-dir",
    help="Path to directory with .pdb or .cif file structures. The IDs need to be in the name of the file i.e id.pdb or id.cif",
    type=click.Path(),
    default=None,
)
@common_options
@compare_options
@bakta_options
def compare(
    ctx,
    input,
    output,
    threads,
    prefix,
    evalue,
    force,
    database,
    sensitivity,
    keep_tmp_files,
    predictions_dir,
    structure_dir,
    max_seqs,
    ultra_sensitive,
    extra_foldseek_params,
    custom_db,
    custom_annotations,
    foldseek_gpu,
    all_proteins,
    euk,
    **kwargs,
):
    """Runs Foldseek vs baktfold db"""

    # validates the directory  (need to before I start baktfold or else no log file is written)

    instantiate_dirs(output, force)

    output: Path = Path(output)
    logdir: Path = Path(output) / "logs"

    params = {
        "--input": input,
        "--output": output,
        "--threads": threads,
        "--force": force,
        "--prefix": prefix,
        "--evalue": evalue,
        "--database": database,
        "--sensitivity": sensitivity,
        "--predictions-dir": predictions_dir,
        "--structure-dir": structure_dir,
        "--keep-tmp-files": keep_tmp_files,
        "--max-seqs": max_seqs,
        "--ultra-sensitive": ultra_sensitive,
        "--extra-foldseek-params": extra_foldseek_params,
        "--custom-db": custom_db,
        "--custom-annotations": custom_annotations,
        "--foldseek-gpu": foldseek_gpu,
        "--all-proteins": all_proteins,
        "--euk": euk
    }

    # initial logging etc
    start_time = begin_baktfold(params, "compare")

    # check foldseek is installed
    check_dependencies()

    # check the database is installed and return it
    database = validate_db(database, DB_DIR, foldseek_gpu)

    # bool for the subcommand
    if (structure_dir):
        structures = True
        if predictions_dir:
            logger.warning(f"Both --predictions-dir {predictions_dir} and --structure-dir {structure_dir} detected")
            logger.warning(f"Proceeding with --predictions-dir {predictions_dir}")
            structures = False
    else:
        structures = False
        if not predictions_dir:
            logger.error(f"neither --predictions_dir or --structure-dir was specified. Please specify one.")


    ###
    # parse the json output and save hypotheticals as AA FASTA
    ###

    fasta_aa: Path = Path(output) / f"{prefix}_aa.fasta"
    data, features, has_duplicate_locus = parse_json_input(input, fasta_aa, all_proteins)

    ###
    # split features in hypotheticals and non hypotheticals
    ###

    if all_proteins:
        hypotheticals = [
            feat for feat in features
            if feat['type'] == bc.FEATURE_CDS
        ]
    else:
        hypotheticals = [
            feat for feat in features
            if feat['type'] == bc.FEATURE_CDS and 'hypothetical' in feat
        ]

    non_hypothetical_features = [
        feat for feat in features
        if feat not in hypotheticals
    ]

    # code to read in and append 3Di from ProstT5 to the dictionary for the json output

    if not structures:
        threedi_aa = Path(predictions_dir) / f"{prefix}_3di.fasta"
        predictions = {record.id: str(record.seq) for record in SeqIO.parse(threedi_aa, "fasta")}
        
        for feat in hypotheticals:
            if has_duplicate_locus:
                seq_id = feat["id"]
            else:
                seq_id = feat["locus"]
            threedi_seq = predictions.get(seq_id)
            feat["3di"] = threedi_seq if threedi_seq else ""

    hypotheticals = subcommand_compare(
        hypotheticals,
        output,
        threads,
        evalue,
        sensitivity,
        database,
        prefix,
        predictions_dir=predictions_dir,
        structures=structures,
        structure_dir=structure_dir,
        logdir=logdir,
        proteins_flag=False,
        max_seqs=max_seqs,
        ultra_sensitive=ultra_sensitive,
        extra_foldseek_params=extra_foldseek_params,
        custom_db=custom_db,
        foldseek_gpu=foldseek_gpu,
        custom_annotations=custom_annotations,
        has_duplicate_locus=has_duplicate_locus
    )


    for cds in hypotheticals:
        anno.combine_annotation(cds)  # add on PSTC annotations and mark hypotheticals

    # recombine updated and existing features
    combined_features = non_hypothetical_features + hypotheticals  # recombine
    
    # Sort by ascending 'id'
    combined_features_sorted = sorted(combined_features, key=lambda x: x.get('id', ''))

    # put back in dictionary
    data['features'] = combined_features_sorted


    # map features by sequence for io
    features_by_sequence = {seq['id']: [] for seq in data['sequences']}

    for feature in data['features']:
        if 'discarded' not in feature:
            seq_features = features_by_sequence.get(feature['sequence'])
            if seq_features is not None:
                seq_features.append(feature)

    # flatten sorted features
    features = []

    for seq in data['sequences']:
        seq_features = features_by_sequence[seq['id']]
        if euk: # ensure gene -> mRNA -> CDS ordering for each locus tag but overall keep the start as the prinary sort

            seq_features.sort(key=sort_euk_feature_key)

        else:
            seq_features.sort(key=lambda k: k['start'])  # sort features by start position
        features.extend(seq_features)

    # overwrite feature list with sorted features
    data['features'] = features

    #####
    # don't include this for now as no gene symbols
    #####

    # logger.info('improve annotations...')
    # genes_with_improved_symbols = anno.select_gene_symbols([feature for feature in features if feature['type'] in [bc.FEATURE_CDS, bc.FEATURE_SORF]])
    # print(f'\trevised gene symbols: {len(genes_with_improved_symbols)}')


    ####
    # bakta output module
    ####
    logger.info('writing baktfold outputs')
    io.write_bakta_outputs(data,features, features_by_sequence, output, prefix, custom_db, euk, has_duplicate_locus)

    # cleanup the temp files
    if not keep_tmp_files:
        clean_up_temporary_files(output, prefix)

    # end baktfold
    end_baktfold(start_time, "compare")




""" 
proteins-predict command
Uses ProstT5 to predict 3Di from a multiFASTA of proteins as input
"""


@main_cli.command()
@click.help_option("--help", "-h")
@click.version_option(get_version(), "--version", "-V")
@click.pass_context
@click.option(
    "-i",
    "--input",
    help="Path to input multiFASTA file",
    type=click.Path(),
    required=True,
)
@common_options
@predict_options
def proteins_predict(
    ctx,
    input,
    output,
    threads,
    prefix,
    force,
    database,
    autotune,
    batch_size,
    cpu,
    omit_probs,
    save_per_residue_embeddings,
    save_per_protein_embeddings,
    mask_threshold,
    **kwargs,
):

    """Runs ProstT5 on a multiFASTA input - GPU recommended"""

    # validates the directory  (need to before I start baktfold or else no log file is written)
    instantiate_dirs(output, force)

    output: Path = Path(output)
    logdir: Path = Path(output) / "logs"

    params = {
        "--input": input,
        "--output": output,
        "--threads": threads,
        "--force": force,
        "--prefix": prefix,
        "--database": database,
        "--autotune": autotune,
        "--batch-size": batch_size,
        "--cpu": cpu,
        "--omit-probs": omit_probs,
        "--save-per-residue-embeddings": save_per_residue_embeddings,
        "--save-per-protein-embeddings": save_per_protein_embeddings,
        "--mask-threshold": mask_threshold,

    }

    # initial logging etc
    start_time = begin_baktfold(params, "proteins-predict")

    # check foldseek is installed
    check_dependencies()

    # check the database is installed and return it
    database = validate_db(database, DB_DIR, foldseek_gpu=False) # dont need foldseek_gpu


    ###
    # parse the json output and save hypotheticals as AA FASTA
    ###

    fasta_aa: Path = Path(output) / f"{prefix}_aa.fasta"
    # puts the CDS AA in a simple dictionary for ProstT5 code
    aas = parse_protein_input(input, fasta_aa)

    # put the CDS AA in a simple dictionary for ProstT5 code
    cds_dict = {}
    for feat in aas:
        cds_dict[feat['locus']] = feat['aa']

    # baktfold predict
    model_dir = database
    model_name = "Rostlab/ProstT5_fp16"
    checkpoint_path = Path(CNN_DIR) / "cnn_chkpnt" / "model.pt"

    if autotune:

        input_path = files("baktfold.features.autotune_data").joinpath("swissprot_5000.fasta.gz")

        step = 20
        min_batch = 1
        max_batch = 301
        sample_seqs = 602

        batch_size = run_autotune(
            input_path,
            model_dir,
            model_name,
            cpu,
            threads,
            step, 
            min_batch,
            max_batch, 
            sample_seqs)

    aas = subcommand_predict(
        aas,
        cds_dict,
        output,
        prefix,
        cpu,
        omit_probs,
        model_dir,
        model_name,
        checkpoint_path,
        batch_size,
        save_per_residue_embeddings=save_per_residue_embeddings,
        save_per_protein_embeddings=save_per_protein_embeddings,
        threads=threads,
        mask_threshold=mask_threshold,
        has_duplicate_locus=False
    )

    # end baktfold
    end_baktfold(start_time, "proteins-predict")


""" 
proteins compare command
Runs Foldseek vs baktfold DBs for multiFASTA 3Di sequences (predicted with proteins-predict)
"""


@main_cli.command()
@click.help_option("--help", "-h")
@click.version_option(get_version(), "--version", "-V")
@click.pass_context
@click.option(
    "-i",
    "--input",
    help="Path to input file in multiFASTA format",
    type=click.Path(),
    required=True,
)
@click.option(
    "--predictions-dir",
    help="Path to output directory from baktfold proteins-predict",
    type=click.Path(),
)
@click.option(
    "--structure-dir",
    help="Path to directory with .pdb or .cif file structures. The CDS IDs need to be in the name of the file",
    type=click.Path(),
)
@common_options
@compare_options
def proteins_compare(
    ctx,
    input,
    output,
    threads,
    prefix,
    evalue,
    force,
    database,
    sensitivity,
    keep_tmp_files,
    predictions_dir,
    structure_dir,
    max_seqs,
    ultra_sensitive,
    extra_foldseek_params,
    custom_db,
    custom_annotations,
    foldseek_gpu,
    **kwargs,
):
    """Runs Foldseek vs baktfold db on proteins input"""

    # validates the directory  (need to before I start baktfold or else no log file is written)

    instantiate_dirs(output, force)

    output: Path = Path(output)
    logdir: Path = Path(output) / "logs"

    params = {
        "--input": input,
        "--output": output,
        "--threads": threads,
        "--force": force,
        "--prefix": prefix,
        "--evalue": evalue,
        "--database": database,
        "--sensitivity": sensitivity,
        "--predictions-dir": predictions_dir,
        "--structure-dir": structure_dir,
        "--keep-tmp-files": keep_tmp_files,
        "--max-seqs": max_seqs,
        "--ultra-sensitive": ultra_sensitive,
        "--extra-foldseek-params": extra_foldseek_params,
        "--custom-db": custom_db,
        "--custom-annotations": custom_annotations,
        "--foldseek-gpu": foldseek_gpu,
    }


    # initial logging etc
    start_time = begin_baktfold(params, "proteins-compare")

    # check foldseek is installed
    check_dependencies()

    # bool for the subcommand
    if (structure_dir):
        structures = True
        if predictions_dir:
            logger.warning(f"Both --predictions-dir {predictions_dir} and --structure-dir {structure_dir} detected")
            logger.warning(f"Proceeding with --predictions-dir {predictions_dir}")
            structures = False
    else:
        structures = False
        if not predictions_dir:
            logger.error(f"neither --predictions-dir or --structure-dir was specified. Please specify one.")

    # check if predictions_dir and structures


    ###
    # parse the json output and save hypotheticals as AA FASTA
    ###

    fasta_aa: Path = Path(output) / f"{prefix}_aa.fasta"
    # puts the CDS AA in a simple dictionary for ProstT5 code
    aas = parse_protein_input(input, fasta_aa)

    # for adding the 3Di to the dictionary
    if predictions_dir:
        threedi_aa : Path = Path(predictions_dir) / f"{prefix}_3di.fasta"
        predictions = {record.id: str(record.seq) for record in SeqIO.parse(threedi_aa, "fasta")}
        for aa in aas:
            seq_id = aa["locus"]
            threedi_seq = predictions.get(seq_id)
            aa["3di"] = threedi_seq


    aas = subcommand_compare(
        aas, 
        output,
        threads,
        evalue,
        sensitivity,
        database,
        prefix,
        predictions_dir,
        structures=structures,
        structure_dir=structure_dir,
        logdir=logdir,
        proteins_flag=True,
        max_seqs=max_seqs,
        ultra_sensitive=ultra_sensitive,
        extra_foldseek_params=extra_foldseek_params,
        custom_db=custom_db,
        foldseek_gpu=foldseek_gpu,
        custom_annotations=custom_annotations,
        has_duplicate_locus=False
    )

    #####
    # update the hypotheticals 
    #####

    for aa in aas:
        anno.combine_annotation(aa)  # add on PSTC annotations and mark hypotheticals



    ####
    # bakta output module
    ####
    logger.info('writing baktfold outputs')

    cfg.run_end = datetime.now()
    run_duration = (cfg.run_end - cfg.run_start).total_seconds()

    for aa in aas:  # reset mock attributes
        aa['start'] = -1
        aa['stop'] = -1

    ############################################################################
    # Write output files
    # - write comprehensive annotation results as JSON
    # - write optional output files in TSV, FAA formats
    # - remove temp directory
    ############################################################################
    
    io.write_bakta_proteins_outputs(aas, output, prefix, custom_db)


    # cleanup the temp files
    if not keep_tmp_files:
        clean_up_temporary_files(output, prefix)

    # end baktfold
    end_baktfold(start_time, "protein-compare")

"""
createdb command
"""


@main_cli.command()
@click.help_option("--help", "-h")
@click.version_option(get_version(), "--version", "-V")
@click.pass_context
@click.option(
    "--fasta-aa",
    help="Path to input Amino Acid FASTA file of proteins",
    type=click.Path(),
    required=True,
)
@click.option(
    "--fasta-3di",
    help="Path to input 3Di FASTA file of proteins",
    type=click.Path(),
    required=True,
)
@click.option(
    "-o",
    "--output",
    default="output_baktfold_foldseek_db",
    show_default=True,
    type=click.Path(),
    help="Output directory ",
)
@click.option(
    "-t",
    "--threads",
    help="Number of threads to use with Foldseek",
    default=1,
    type=int,
    show_default=True,
)
@click.option(
    "-p",
    "--prefix",
    default="baktfold_foldseek_db",
    help="Prefix for Foldseek database",
    type=str,
    show_default=True,
)
@click.option(
    "-f",
    "--force",
    is_flag=True,
    help="Force overwrites the output directory",
)
def createdb(
    ctx,
    fasta_aa,
    fasta_3di,
    output,
    threads,
    prefix,
    force,
    **kwargs,
):
    """Creates foldseek DB from AA FASTA and 3Di FASTA input files"""

    # validates the directory  (need to before I start phold or else no log file is written)
    instantiate_dirs(output, force)

    output: Path = Path(output)
    logdir: Path = Path(output) / "logs"

    params = {
        "--fasta-aa": fasta_aa,
        "--fasta-3di": fasta_3di,
        "--output": output,
        "--threads": threads,
        "--force": force,
        "--prefix": prefix,
    }

    # initial logging etc
    start_time = begin_baktfold(params, "createdb")

    # check foldseek is installed
    check_dependencies()

    logger.info(f"Creating the Foldseek database using {fasta_aa} and {fasta_3di}.")
    logger.info(
        f"The database will be saved in the {output} directory and be called {prefix}."
    )

    ############
    # create foldseek db
    ############

    foldseek_query_db_path: Path = Path(output)
    foldseek_query_db_path.mkdir(parents=True, exist_ok=True)

    generate_foldseek_db_from_aa_3di(
        fasta_aa, fasta_3di, foldseek_query_db_path, logdir, prefix
    )

    # end
    end_baktfold(start_time, "createdb")



"""
convert  options
"""

def convert_options(func):
    """compare command line args"""
    options = [
        click.option(
            "-i",
            "--input",
            help="Path to input Prokka GenBank (.gbk) output",
            type=click.Path(),
            required=True,
        ),
        click.option(
            "-o",
            "--outfile",
            default="converted_bakta_formatted.json",
            show_default=True,
            type=click.Path(),
            help="Output bakta format .json output",
        ),
        click.option(
            "-f",
            "--force",
            is_flag=True,
            help="Force overwrites the output file",
        )
    ]
    for option in reversed(options):
        func = option(func)
    return func
    


"""
convert prokka GenBank to Bakta formatted json
"""

@main_cli.command()
@click.help_option("--help", "-h")
@click.version_option(get_version(), "--version", "-V")
@click.pass_context
@convert_options
def convert_prokka(
    ctx,
    input,
    outfile,
    force,
    **kwargs,
):
    """Converts Prokka GenBank to Bakta format json"""

    # validates the output file - check it doesnt exist, if it does overwrite it
    validate_outfile(outfile, force)

    # validates input genbank and returns records
    records = check_genbank_and_prokka(input, euk=False)


    params = {
        "--input": input,
        "--outfile": outfile,
        "--force": force,
    }

    # initial logging etc
    start_time = begin_baktfold(params, "convert-prokka", no_log=True)

    # check foldseek is installed
    # check_dependencies()

    logger.info(f"Converting Prokka input GenBank file {input} to Bakta formatted .json file.")
    logger.info(
        f"This will be saved as {outfile}."
    )

    prokka_gbk_to_json(records, outfile)

    logger.info(f"Conversion successful.")
    logger.info(f"Bakta format JSON → {outfile}")
 
    # end
    end_baktfold(start_time, "convert-prokka")


"""
convert eukaryotic GenBank to Bakta formatted json
"""

@main_cli.command()
@click.help_option("--help", "-h")
@click.version_option(get_version(), "--version", "-V")
@click.pass_context
@convert_options
def convert_euk(
    ctx,
    input,
    outfile,
    force,
    **kwargs,
):
    """(Experimental) Converts eukaryotic GenBank to Bakta format json"""

    # validates the output file - check it doesnt exist, if it does overwrite it
    validate_outfile(outfile, force)

    # validates input genbank and returns records
    records = check_genbank_and_prokka(input, euk=True)


    params = {
        "--input": input,
        "--outfile": outfile,
        "--force": force,
    }

    # initial logging etc
    start_time = begin_baktfold(params, "convert-euk", no_log=True)

    # check foldseek is installed
    # check_dependencies()

    logger.info(f"Converting eukaryotic input GenBank file {input} to Bakta formatted .json file.")
    logger.info(
        f"This will be saved as {outfile}."
    )

    eukaryotic_gbk_to_json(records, outfile)

    logger.info(f"Conversion successful.")
    logger.info(f"Bakta format JSON → {outfile}")
 
    # end
    end_baktfold(start_time, "convert-euk")





"""
install command
"""


@main_cli.command()
@click.help_option("--help", "-h")
@click.version_option(get_version(), "--version", "-V")
@click.pass_context
@click.option(
    "-d",
    "--database",
    type=str,
    default=None,
    help="Specific path to install the baktfold database",
)
@click.option(
    "--foldseek-gpu",
    is_flag=True,
    help="Use this to enable compatibility with Foldseek-GPU acceleration",
)
@click.option(
            "-t",
            "--threads",
            help="Number of threads",
            default=1,
            type=int,
            show_default=True,
)
def install(
    ctx,
    database,
    foldseek_gpu,
    threads,
    **kwargs,
):
    """Installs ProstT5 model and baktfold database"""

    if database:
        logger.info(
            f"You have specified the {database} directory to store the baktfold database and ProstT5 model"
        )
        database: Path = Path(database)
    else:
        logger.info(
            f"Downloading the baktfold database into the default directory {DB_DIR}"
        )
        database = Path(DB_DIR)

    model_name = "Rostlab/ProstT5_fp16"

    logger.info(
        f"Checking that the {model_name} ProstT5 model is available in {database}"
    )

    # always install with cpu mode as guaranteed to be present
    cpu = True

    # load model (will be downloaded if not present)
    model, vocab = get_T5_model(database, model_name, cpu, threads=1)
    del model
    del vocab
    logger.info(f"ProstT5 model downloaded")

    # will check if db is present, and if not, download it
    install_database(database, foldseek_gpu, threads)

@main_cli.command()
@click.help_option("--help", "-h")
@click.version_option(get_version(), "--version", "-V")
@click.pass_context
@click.option(
    "-i",
    "--input",
    help="Optional path to input file of proteins if you do not want to use the default sample of 5000 Phold DB proteins",
    type=click.Path()
)
@click.option(
    "--cpu",
    is_flag=True,
    help="Use cpus only.",
)
@click.option(
    "-t",
    "--threads",
    help="Number of threads",
    default=1,
    type=int,
    show_default=True,
)
@click.option(
    "-d",
    "--database",
    type=str,
    default=None,
    help="Specific path to installed phold database",
)
@click.option(
    "--min_batch",
    show_default=True,
    type=int,
    default=1,
    help="Minimum batch size to test",
)
@click.option(
    "--step",
    show_default=True,
    type=int,
    default=10,
    help="Controls batch size step increment",
)
@click.option(
    "--max_batch",
    default=251,
    show_default=True,
    type=int,
    help="Maximum batch size to test",
)
@click.option(
    "--sample_seqs",
    default=500,
    show_default=True,
    type=int,
    help="Number of proteins to subsample from input.",
)
def autotune(
    ctx,
    input,
    cpu,
    threads,
    database,
    step,
    min_batch,
    max_batch,
    sample_seqs,
    **kwargs,
):
    """Determines optimal batch size for 3Di prediction with your hardware"""

    params = {
        "--input": input,
        "--threads": threads,
        "--cpu": cpu,
        "--database": database,
        "--step": step,
        "--min_batch": min_batch,
        "--max_batch": max_batch,
        "--sample_seqs": sample_seqs,
    }

    # initial logging etc
    start_time = begin_baktfold(params, "autotune", no_log=True)

    # check the database is installed
    database = validate_db(database, DB_DIR, foldseek_gpu=False)

    if input:
        input_path = input
    else:
        input_path = files("baktfold.features.autotune_data").joinpath("swissprot_5000.fasta.gz")

    model_dir = database
    model_name = "Rostlab/ProstT5_fp16"

    batch_size = run_autotune(
        input_path,
        model_dir,
        model_name,
        cpu,
        threads,
        step, 
        min_batch,
        max_batch, 
        sample_seqs)


@click.command()
def citation(**kwargs):
    """Print the citation(s) for this tool"""
    print_citation()


# main_cli.add_command(run)
main_cli.add_command(citation)


def main():
    main_cli()


if __name__ == "__main__":
    main()
