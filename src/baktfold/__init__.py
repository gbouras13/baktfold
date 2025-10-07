#!/usr/bin/env python3

__version__ = '0.0.1'

import gzip
from pathlib import Path

import click
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation, SeqFeature
from loguru import logger

import baktfold.bakta.constants as bc
from baktfold.bakta.json_io import parse_json_input
from baktfold.databases.db import install_database, validate_db
from baktfold.features.create_foldseek_db import generate_foldseek_db_from_aa_3di
from baktfold.features.predict_3Di import get_T5_model
from baktfold.io.handle_genbank import open_protein_fasta_file
from baktfold.subcommands.compare import subcommand_compare
from baktfold.subcommands.predict import subcommand_predict
from baktfold.utils.constants import DB_DIR, CNN_DIR
from baktfold.utils.util import (begin_baktfold, clean_up_temporary_files, end_baktfold,
                              get_version, print_citation)
from baktfold.utils.validation import (check_dependencies, instantiate_dirs,
                                    validate_input)
import baktfold.io.gff as gff

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
            "--batch_size",
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
            "--omit_probs",
            is_flag=True,
            help="Do not output per residue 3Di probabilities from ProstT5. Mean per protein 3Di probabilities will always be output.",
        ),
        click.option(
            "--save_per_residue_embeddings",
            is_flag=True,
            help="Save the ProstT5 embeddings per resuide in a h5 file ",
        ),
        click.option(
            "--save_per_protein_embeddings",
            is_flag=True,
            help="Save the ProstT5 embeddings as means per protein in a h5 file",
        ),
        click.option(
            "--mask_threshold",
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
            "--keep_tmp_files",
            is_flag=True,
            help="Keep temporary intermediate files, particularly the large foldseek_results.tsv of all Foldseek hits",
        ),
        click.option(
            "--max_seqs",
            type=int,
            default=1000,
            show_default=True,
            help="Maximum results per query sequence allowed to pass the prefilter. You may want to reduce this to save disk space for enormous datasets",
        ),
        click.option(
            "--ultra_sensitive",
            is_flag=True,
            help="Runs baktfold with maximum sensitivity by skipping Foldseek prefilter. Not recommended for large datasets.",
        ),
        click.option(
            "--extra_foldseek_params",
            type=str,
            help="Extra foldseek search params"
        ),
        click.option(
            "--custom_db",
            type=str,
            help="Path to custom database"
        ),
        click.option(
            "--foldseek_gpu",
            is_flag=True,
            help="Use this to enable compatibility with Foldseek-GPU search acceleration",
        )
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
def run(
    ctx,
    input,
    output,
    threads,
    prefix,
    evalue,
    force,
    database,
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
        "--batch_size": batch_size,
        "--sensitivity": sensitivity,
        "--keep_tmp_files": keep_tmp_files,
        "--cpu": cpu,
        "--omit_probs": omit_probs,
        "--max_seqs": max_seqs,
        "--save_per_residue_embeddings": save_per_residue_embeddings,
        "--save_per_protein_embeddings": save_per_protein_embeddings,
        "--ultra_sensitive": ultra_sensitive,
        "--mask_threshold": mask_threshold,
        "--extra_foldseek_params": extra_foldseek_params,
        "--custom_db": custom_db,
        "--foldseek_gpu": foldseek_gpu,
    }

    # initial logging etc
    start_time = begin_baktfold(params, "run")

    # check foldseek is installed
    check_dependencies()

    # check the database is installed and return it
    #database = validate_db(database, DB_DIR, foldseek_gpu)

    # validate input
    #fasta_flag, gb_dict, method = validate_input(input, threads)

    # baktfold predict
    model_dir = database
    model_name = "Rostlab/ProstT5_fp16"
    checkpoint_path = Path(CNN_DIR) / "cnn_chkpnt" / "model.pt"

    fasta_aa: Path = Path(output) / f"{prefix}_aa.fasta"

    data, features = parse_json_input(input, fasta_aa)

    hypotheticals = [feat for feat in features if feat['type'] == bc.FEATURE_CDS and 'hypothetical' in feat]
    non_hypothetical_features = [
    feat for feat in features
    if (feat['type'] != bc.FEATURE_CDS) or 
       (feat['type'] == bc.FEATURE_CDS and 'hypothetical' not in (feat.get('product') or '').lower())
]
    
    # for prostT5 code
    cds_dict = {}
    for feat in hypotheticals:
        cds_dict[feat['locus']] = feat['aa']


    # add a function to add 3Di to cds_dict

    subcommand_predict(
        cds_dict,
        output,
        prefix,
        cpu,
        omit_probs,
        model_dir,
        model_name,
        checkpoint_path,
        batch_size,
        proteins_flag=False,
        save_per_residue_embeddings=save_per_residue_embeddings,
        save_per_protein_embeddings=save_per_protein_embeddings,
        threads=threads,
        mask_threshold=mask_threshold,
    )

    # baktfold compare

    # predictions_dir is output as this will be where it lives

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
    )

    combined_features = non_hypothetical_features + hypotheticals  # recombine

    # Sort by ascending 'id'
    combined_features_sorted = sorted(combined_features, key=lambda x: x.get('id', ''))

    data['features'] = combined_features_sorted

    print('select features and preserve existing IDs...')
    logger.info('start feature selection and preserving existing feature IDs')

    # map features by sequence
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
        seq_features.sort(key=lambda k: k['start'])  # sort features by start position
        features.extend(seq_features)

    # overwrite feature list with sorted features
    data['features'] = features

    logger.info('selected features=%i', len(features))
    print(f'\tselected: {len(features)}')


    print('\tGFF3...')
    gff3_path: Path = Path(output) / f"{prefix}.gff3"
    gff.write_features(data, features_by_sequence, gff3_path)


    # end baktfold
    end_baktfold(start_time, "run")






"""
predict command
Uses ProstT5 to predict 3Di sequences 
"""


# @main_cli.command()
# @click.help_option("--help", "-h")
# @click.version_option(get_version(), "--version", "-V")
# @click.pass_context
# @click.option(
#     "-i",
#     "--input",
#     help="Path to input file in Genbank format or nucleotide FASTA format",
#     type=click.Path(),
#     required=True,
# )
# @common_options
# @predict_options
# def predict(
#     ctx,
#     input,
#     output,
#     threads,
#     prefix,
#     force,
#     database,
#     batch_size,
#     cpu,
#     omit_probs,
#     save_per_residue_embeddings,
#     save_per_protein_embeddings,
#     mask_threshold,
#     finetune,
#     vanilla,
#     hyps,
#     **kwargs,
# ):
#     """Uses ProstT5 to predict 3Di tokens - GPU recommended"""

#     # validates the directory  (need to before I start baktfold or else no log file is written)
#     instantiate_dirs(output, force)

#     output: Path = Path(output)
#     logdir: Path = Path(output) / "logs"

#     params = {
#         "--input": input,
#         "--output": output,
#         "--threads": threads,
#         "--force": force,
#         "--prefix": prefix,
#         "--database": database,
#         "--batch_size": batch_size,
#         "--cpu": cpu,
#         "--omit_probs": omit_probs,
#         "--save_per_residue_embeddings": save_per_residue_embeddings,
#         "--save_per_protein_embeddings": save_per_protein_embeddings,
#         "--mask_threshold": mask_threshold,
#         "--finetune": finetune,
#         "--vanilla": vanilla,
#         "--hyps": hyps
#     }

#     # initial logging etc
#     start_time = begin_baktfold(params, "predict")

#     # check the database is installed
#     database = validate_db(database, DB_DIR, foldseek_gpu=False)

#     # validate input
#     fasta_flag, gb_dict, method = validate_input(input, threads)

#     # runs baktfold predict subcommand
#     model_dir = database
#     model_name = "Rostlab/ProstT5_fp16"
#     checkpoint_path = Path(CNN_DIR) / "cnn_chkpnt" / "model.pt"

#     if finetune:
#         model_name = "gbouras13/ProstT5baktfold"
#         checkpoint_path = Path(CNN_DIR) / "cnn_chkpnt_finetune" / "baktfold_db_model.pth"
#         if vanilla:
#             checkpoint_path = Path(CNN_DIR) / "cnn_chkpnt_finetune" / "vanilla_model.pth"


#     subcommand_predict(
#         gb_dict,
#         method,
#         output,
#         prefix,
#         cpu,
#         omit_probs,
#         model_dir,
#         model_name,
#         checkpoint_path,
#         batch_size,
#         proteins_flag=False,
#         fasta_flag=fasta_flag,
#         save_per_residue_embeddings=save_per_residue_embeddings,
#         save_per_protein_embeddings=save_per_protein_embeddings,
#         threads=threads,
#         mask_threshold=mask_threshold,
#         hyps=hyps
#     )

#     # end baktfold
#     end_baktfold(start_time, "predict")


"""
compare command
"""


# @main_cli.command()
# @click.help_option("--help", "-h")
# @click.version_option(get_version(), "--version", "-V")
# @click.pass_context
# @click.option(
#     "-i",
#     "--input",
#     help="Path to input file in Genbank format or nucleotide FASTA format",
#     type=click.Path(),
#     required=True,
# )
# @click.option(
#     "--predictions_dir",
#     help="Path to output directory from baktfold predict",
#     type=click.Path(),
# )
# @click.option(
#     "--structures",
#     is_flag=True,
#     help="Use if you have .pdb or .cif file structures for the input proteins (e.g. with AF2/Colabfold .pdb or AF3 for .cif) in a directory that you specify with --structure_dir",
# )
# @click.option(
#     "--structure_dir",
#     help="Path to directory with .pdb or .cif file structures. The CDS IDs need to be in the name of the file",
#     type=click.Path(),
# )
# @click.option(
#     "--filter_structures",
#     is_flag=True,
#     help="Flag that creates a copy of the .pdb or .cif files structures with matching record IDs found in the input GenBank file. Helpful if you have a directory with lots of .pdb files and want to annotate only e.g. 1 phage.",
# )
# @common_options
# @compare_options
# def compare(
#     ctx,
#     input,
#     output,
#     threads,
#     prefix,
#     evalue,
#     force,
#     database,
#     sensitivity,
#     predictions_dir,
#     structures,
#     structure_dir,
#     filter_structures,
#     keep_tmp_files,
#     card_vfdb_evalue,
#     separate,
#     max_seqs,
#     ultra_sensitive,
#     extra_foldseek_params,
#     custom_db,
#     foldseek_gpu,
#     **kwargs,
# ):
#     """Runs Foldseek vs baktfold db"""

#     # validates the directory  (need to before I start baktfold or else no log file is written)

#     instantiate_dirs(output, force)

#     output: Path = Path(output)
#     logdir: Path = Path(output) / "logs"

#     params = {
#         "--input": input,
#         "--output": output,
#         "--threads": threads,
#         "--force": force,
#         "--prefix": prefix,
#         "--evalue": evalue,
#         "--database": database,
#         "--sensitivity": sensitivity,
#         "--predictions_dir": predictions_dir,
#         "--structures": structures,
#         "--structure_dir": structure_dir,
#         "--filter_structures": filter_structures,
#         "--keep_tmp_files": keep_tmp_files,
#         "--card_vfdb_evalue": card_vfdb_evalue,
#         "--separate": separate,
#         "--max_seqs": max_seqs,
#         "--ultra_sensitive": ultra_sensitive,
#         "--extra_foldseek_params": extra_foldseek_params,
#         "--custom_db": custom_db,
#         "--foldseek_gpu": foldseek_gpu,
#     }

#     # initial logging etc
#     start_time = begin_baktfold(params, "compare")

#     # check foldseek is installed
#     check_dependencies()

#     # check the database is installed
#     database = validate_db(database, DB_DIR, foldseek_gpu)

#     # validate fasta
#     fasta_flag, gb_dict, method = validate_input(input, threads)

#     subcommand_compare(
#         gb_dict,
#         output,
#         threads,
#         evalue,
#         card_vfdb_evalue,
#         sensitivity,
#         database,
#         prefix,
#         predictions_dir,
#         structures,
#         structure_dir,
#         logdir,
#         filter_structures,
#         remote_flag=False,
#         proteins_flag=False,
#         fasta_flag=fasta_flag,
#         separate=separate,
#         max_seqs=max_seqs,
#         ultra_sensitive=ultra_sensitive,
#         extra_foldseek_params=extra_foldseek_params,
#         custom_db=custom_db,
#         foldseek_gpu=foldseek_gpu,
#     )

#     # cleanup the temp files
#     if keep_tmp_files is False:
#         clean_up_temporary_files(output)

#     # end baktfold
#     end_baktfold(start_time, "compare")






# """ 
# proteins-predict command
# Uses ProstT5 to predict 3Di from a multiFASTA of proteins as input
# """


# @main_cli.command()
# @click.help_option("--help", "-h")
# @click.version_option(get_version(), "--version", "-V")
# @click.pass_context
# @click.option(
#     "-i",
#     "--input",
#     help="Path to input multiFASTA file",
#     type=click.Path(),
#     required=True,
# )
# @common_options
# @predict_options
# def proteins_predict(
#     ctx,
#     input,
#     output,
#     threads,
#     prefix,
#     force,
#     database,
#     batch_size,
#     cpu,
#     omit_probs,
#     save_per_residue_embeddings,
#     save_per_protein_embeddings,
#     mask_threshold,
#     finetune,
#     vanilla,
#     **kwargs,
# ):
#     """Runs ProstT5 on a multiFASTA input - GPU recommended"""

#     # validates the directory  (need to before baktfold starts or else no log file is written)
#     instantiate_dirs(output, force)

#     output: Path = Path(output)
#     logdir: Path = Path(output) / "logs"

#     params = {
#         "--input": input,
#         "--output": output,
#         "--threads": threads,
#         "--force": force,
#         "--prefix": prefix,
#         "--database": database,
#         "--batch_size": batch_size,
#         "--cpu": cpu,
#         "--omit_probs": omit_probs,
#         "--save_per_residue_embeddings": save_per_residue_embeddings,
#         "--save_per_protein_embeddings": save_per_protein_embeddings,
#         "--mask_threshold": mask_threshold,
#         "--finetune": finetune,
#         "--vanilla": vanilla
#     }

#     # initial logging etc
#     start_time = begin_baktfold(params, "proteins-predict")

#     # check the database is installed
#     database = validate_db(database, DB_DIR, foldseek_gpu=False)

#     # Dictionary to store the records
#     cds_dict = {}
#     # need a dummmy nested dict
#     cds_dict["proteins"] = {}

#     # Iterate through the multifasta file and save each Seqfeature to the dictionary
#     # 1 dummy record = proteins

#     with open_protein_fasta_file(input) as handle:  # handles gzip too
#         records = list(SeqIO.parse(handle, "fasta"))
#         if not records:
#             logger.warning(f"No proteins were found in your input file {input}.")
#             logger.error(
#                 f"Your input file {input} is likely not a amino acid FASTA file. Please check this."
#             )
#         for record in records:
#             prot_id = record.id
#             feature_location = FeatureLocation(0, len(record.seq))
#             # Seq needs to be saved as the first element in list hence the closed brackets [str(record.seq)]
#             seq_feature = SeqFeature(
#                 feature_location,
#                 type="CDS",
#                 qualifiers={
#                     "ID": record.id,
#                     "description": record.description,
#                     "translation": str(record.seq),
#                 },
#             )

#             cds_dict["proteins"][prot_id] = seq_feature

#     if not cds_dict:
#         logger.error(f"Error: no AA protein sequences found in {input} file")

#     # runs baktfold predict subcommand
#     model_dir = database
#     model_name = "Rostlab/ProstT5_fp16"
#     checkpoint_path = Path(CNN_DIR) / "cnn_chkpnt" / "model.pt"

#     if finetune:
#         model_name = "gbouras13/ProstT5baktfold"
#         checkpoint_path = Path(CNN_DIR) / "cnn_chkpnt_finetune" / "baktfold_db_model.pth"
#         if vanilla:
#             checkpoint_path = Path(CNN_DIR) / "cnn_chkpnt_finetune" / "vanilla_model.pth"

#     method = "pharokka" # this can be whatever for proteins, it wont matter - it is for genbank input

#     subcommand_predict(
#         cds_dict,
#         method,
#         output,
#         prefix,
#         cpu,
#         omit_probs,
#         model_dir,
#         model_name,
#         checkpoint_path,
#         batch_size,
#         proteins_flag=True,
#         fasta_flag=False,
#         save_per_residue_embeddings=save_per_residue_embeddings,
#         save_per_protein_embeddings=save_per_protein_embeddings,
#         threads=threads,
#         mask_threshold=mask_threshold,
#         hyps=False # always False for this as no Pharokka genbank to parse on input
#     )

#     # end baktfold
#     end_baktfold(start_time, "proteins-predict")


# """ 
# proteins compare command
# Runs Foldseek vs baktfold DB for multiFASTA 3Di sequences (predicted with proteins-predict)
# """


# @main_cli.command()
# @click.help_option("--help", "-h")
# @click.version_option(get_version(), "--version", "-V")
# @click.pass_context
# @click.option(
#     "-i",
#     "--input",
#     help="Path to input file in multiFASTA format",
#     type=click.Path(),
#     required=True,
# )
# @click.option(
#     "--predictions_dir",
#     help="Path to output directory from baktfold proteins-predict",
#     type=click.Path(),
# )
# @click.option(
#     "--structures",
#     is_flag=True,
#     help="Use if you have .pdb or .cif file structures for the input proteins (e.g. with AF2/Colabfold) in a directory that you specify with --structure_dir",
# )
# @click.option(
#     "--structure_dir",
#     help="Path to directory with .pdb or .cif file structures. The CDS IDs need to be in the name of the file",
#     type=click.Path(),
# )
# @click.option(
#     "--filter_structures",
#     is_flag=True,
#     help="Flag that creates a copy of the .pdb or .cif files structures with matching record IDs found in the input GenBank file. Helpful if you have a directory with lots of .pdb files and want to annotate only e.g. 1 phage.",
# )
# @common_options
# @compare_options
# def proteins_compare(
#     ctx,
#     input,
#     output,
#     threads,
#     prefix,
#     evalue,
#     force,
#     database,
#     sensitivity,
#     predictions_dir,
#     structures,
#     structure_dir,
#     filter_structures,
#     keep_tmp_files,
#     card_vfdb_evalue,
#     max_seqs,
#     ultra_sensitive,
#     extra_foldseek_params,
#     custom_db,
#     foldseek_gpu,
#     **kwargs
# ):
#     """Runs Foldseek vs baktfold db on proteins input"""

#     # validates the directory  (need to before I start baktfold or else no log file is written)

#     instantiate_dirs(output, force)

#     output: Path = Path(output)
#     logdir: Path = Path(output) / "logs"

#     params = {
#         "--input": input,
#         "--output": output,
#         "--threads": threads,
#         "--force": force,
#         "--prefix": prefix,
#         "--evalue": evalue,
#         "--database": database,
#         "--sensitivity": sensitivity,
#         "--predictions_dir": predictions_dir,
#         "--structures": structures,
#         "--structure_dir": structure_dir,
#         "--filter_structures": filter_structures,
#         "--keep_tmp_files": keep_tmp_files,
#         "--card_vfdb_evalue": card_vfdb_evalue,
#         "--max_seqs": max_seqs,
#         "--ultra_sensitive": ultra_sensitive,
#         "--extra_foldseek_params": extra_foldseek_params,
#         "--custom_db": custom_db,
#         "--foldseek_gpu": foldseek_gpu,
#     }

#     # initial logging etc
#     start_time = begin_baktfold(params, "proteins-compare")

#     # check foldseek is installed
#     check_dependencies()

#     # check the database is installed
#     database = validate_db(database, DB_DIR, foldseek_gpu)

#     # Dictionary to store the records
#     cds_dict = {}
#     # need a dummmy nested dict
#     cds_dict["proteins"] = {}

#     # Iterate through the multifasta file and save each Seqfeature to the dictionary
#     # 1 dummy record = proteins
#     with open_protein_fasta_file(input) as handle:  # handles gzip too
#         records = list(SeqIO.parse(handle, "fasta"))
#         if not records:
#             logger.warning(f"No proteins were found in your input file {input}.")
#             logger.error(
#                 f"Your input file {input} is likely not a amino acid FASTA file. Please check this."
#             )
#         for record in records:
#             prot_id = record.id
#             feature_location = FeatureLocation(0, len(record.seq))
#             # Seq needs to be saved as the first element in list hence the closed brackets [str(record.seq)]
#             seq_feature = SeqFeature(
#                 feature_location,
#                 type="CDS",
#                 qualifiers={
#                     "ID": record.id,
#                     "description": record.description,
#                     "translation": str(record.seq),
#                 },
#             )

#             cds_dict["proteins"][prot_id] = seq_feature

#     if not cds_dict:
#         logger.error(f"Error: no AA protein sequences found in {input} file")

#     success = subcommand_compare(
#         cds_dict,
#         output,
#         threads,
#         evalue,
#         card_vfdb_evalue,
#         sensitivity,
#         database,
#         prefix,
#         predictions_dir,
#         structures,
#         structure_dir,
#         logdir,
#         filter_structures,
#         remote_flag=False,
#         proteins_flag=True,
#         fasta_flag=False,
#         separate=False,
#         max_seqs=max_seqs,
#         ultra_sensitive=ultra_sensitive,
#         extra_foldseek_params=extra_foldseek_params,
#         custom_db=custom_db,
#         foldseek_gpu=foldseek_gpu,
#     )

#     # cleanup the temp files
#     if keep_tmp_files is False:
#         clean_up_temporary_files(output)

#     # end baktfold
#     end_baktfold(start_time, "proteins-compare")

"""
install command
"""


# @main_cli.command()
# @click.help_option("--help", "-h")
# @click.version_option(get_version(), "--version", "-V")
# @click.pass_context
# @click.option(
#     "-d",
#     "--database",
#     type=str,
#     default=None,
#     help="Specific path to install the baktfold database",
# )
# @click.option(
#     "--foldseek_gpu",
#     is_flag=True,
#     help="Use this to enable compatibility with Foldseek-GPU acceleration",
# )
# @click.option(
#     "--extended_db",
#     is_flag=True,
#     help=(
#         "Download the extended baktfold DB 3.16M including 1.8M efam and enVhog proteins without functional labels\n"
#         "instead of the default baktfold Search 1.36M. Using the extended database will likely marginally reduce\n"
#         "functional annotation sensitivity and increase runtime, but may find more hits overall\n"
#         "i.e. including to efam and enVhog proteins that have no functional labels."
#     )
# )
# @click.option(
#             "-t",
#             "--threads",
#             help="Number of threads",
#             default=1,
#             type=int,
#             show_default=True,
# )

# def install(
#     ctx,
#     database,
#     foldseek_gpu,
#     extended_db,
#     threads,
#     **kwargs,
# ):
#     """Installs ProstT5 model and baktfold database"""

#     if database is not None:
#         logger.info(
#             f"You have specified the {database} directory to store the baktfold database and ProstT5 model"
#         )
#         database: Path = Path(database)
#     else:
#         logger.info(
#             f"Downloading the baktfold database into the default directory {DB_DIR}"
#         )
#         database = Path(DB_DIR)

#     model_name = "Rostlab/ProstT5_fp16"

#     logger.info(
#         f"Checking that the {model_name} ProstT5 model is available in {database}"
#     )

#     # always install with cpu mode as guarantee to be present
#     cpu = True

#     # load model (will be downloaded if not present)
#     model, vocab = get_T5_model(database, model_name, cpu, threads=1)
#     del model
#     del vocab
#     logger.info(f"ProstT5 model downloaded")

#     # will check if db is present, and if not, download it
#     install_database(database, foldseek_gpu, extended_db, threads)


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
