#!/usr/bin/env python3

from pathlib import Path

from loguru import logger

from baktfold.features.predict_3Di import get_embeddings

def mask_low_confidence_aa(sequence, scores, threshold=0.5):
    """
    Masks all low confidence AA to X if their corresponding ProstT5 confidence score is below the given threshold.

    Parameters:
    sequence (str): The amino acid sequence.
    scores (List[float]): A list of confidence scores for each amino acid.
    threshold (float, optional): The confidence threshold below which amino acids are converted to lowercase. Default is 0.5.

    Returns:
    str: The modified amino acid sequence with low-confidence residues in lowercase.
    """
    return "".join('X' if float(score) < threshold else aa 
                   for aa, score in zip(sequence, *scores))

def subcommand_predict(
    cds_dict: dict,
    output: Path,
    prefix: str,
    cpu: bool,
    omit_probs: bool,
    model_dir: Path,
    model_name: str,
    checkpoint_path: Path,
    batch_size: int,
    proteins_flag: bool,
    save_per_residue_embeddings: bool,
    save_per_protein_embeddings: bool,
    threads: int,
    mask_threshold: float,
) -> bool:
    """
    Wrapper command for baktfold predict. Predicts embeddings using ProstT5 encoder + CNN prediction head.

    Args:
        cds_dict (Dict[str, any]): id:aa dictionary
        output (str): Output directory path.
        prefix (str): Prefix for output file names.
        cpu (bool): Flag indicating whether to use CPU for prediction.
        omit_probs (bool): Flag indicating whether to omit prediction probabilities from ProstT5.
        model_dir (str): Directory containing the ProstT5 model.
        model_name (str): Name of the ProstT5 model.
        checkpoint_path (Path): Path to ProstT5 CNN checkpoint.
        batch_size (int): Batch size for prediction.
        proteins_flag (bool): True if baktfold proteins-predict, false otherwise
        save_per_residue_embeddings (bool, optional): Whether to save per residue embeddings to h5 file. Defaults to False.
        save_per_protein_embeddings (bool, optional): Whether to save mean per protein embeddings to h5 file. Defaults to False.

    Returns:
        bool: True if prediction succeeds, False otherwise.
    """

    fasta_aa: Path = Path(output) / f"{prefix}_aa.fasta"

    ############
    # prostt5
    ############
   
    fasta_3di: Path = Path(output) / f"{prefix}_3di.fasta"
    # embeddings h5 - will only be generated if flag is true
    output_h5_per_residue: Path = Path(output) / f"{prefix}_embeddings_per_residue.h5"
    output_h5_per_protein: Path = Path(output) / f"{prefix}_embeddings_per_protein.h5"

    if cpu is True:
        half_precision = False
    else:
        half_precision = True

    if omit_probs:
        output_probs = False
    else:
        output_probs = True

    predictions = get_embeddings(
        cds_dict,
        output,
        prefix,
        model_dir,
        model_name,
        checkpoint_path,
        fasta_3di,
        output_h5_per_residue,
        output_h5_per_protein,
        half_precision=half_precision,
        max_residues=5000,
        max_seq_len=1000,
        max_batch=batch_size,
        cpu=cpu,
        output_probs=output_probs,
        proteins_flag=proteins_flag,
        save_per_residue_embeddings=save_per_residue_embeddings,
        save_per_protein_embeddings=save_per_protein_embeddings,
        threads=threads,
        mask_threshold=mask_threshold
    )

    mask_prop_threshold = mask_threshold/100


    ########
    ## write the AA CDS to file
    ######

    # with open(fasta_aa, "w+") as out_f:
    #     for contig_id, rest in cds_dict.items():
    #         aa_contig_dict = cds_dict[contig_id]
    #         prediction_contig_dict = predictions[contig_id]
    #         prediction_contig_dict = {
    #             k: v for k, v in prediction_contig_dict.items() if len(v[0]) > 0
    #         }
    #         for seq_id, cds_feature in aa_contig_dict.items():
    #             if proteins_flag is True:
    #                 out_f.write(f">{seq_id}\n")
    #             else:
    #                 out_f.write(f">{contig_id}:{seq_id}\n")

    #             prot_seq = cds_feature.qualifiers['translation']
    #             # prediction_contig_dict[seq_id][2] these are teh ProstT5 confidence scores from 0-1 - need to convert to list

    #             try:
    #                 # this will fail if ProstT5 OOM fails (or fails for some other reason)
    #                 prot_seq = mask_low_confidence_aa(prot_seq, prediction_contig_dict[seq_id][2].tolist(), threshold=mask_prop_threshold)
    #             except (KeyError, IndexError):
    #                 # in that case, just return 'X' aka masked proteins
    #                 prot_seq = "X" * len(prot_seq)

    #             out_f.write(f"{prot_seq}\n")

    return True
