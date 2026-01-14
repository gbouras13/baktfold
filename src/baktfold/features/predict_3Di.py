#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Code adapted from @mheinzinger 

https://github.com/mheinzinger/ProstT5/blob/main/scripts/predict_3Di_encoderOnly.py

"""

import csv
import json
import h5py
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

from baktfold.databases.db import check_prostT5_download, download_zenodo_prostT5
from baktfold.utils.constants import CNN_DIR

import numpy as np
import torch
import torch.nn.functional as F
from loguru import logger
from torch import nn
from transformers import T5EncoderModel, T5Tokenizer
from tqdm import tqdm

# Convolutional neural network (two convolutional layers)
class CNN(nn.Module):
    """
    Convolutional neural network (two convolutional layers).

    Args:
      nn.Module: The base class for all neural network modules.

    Returns:
      None.

    Examples:
      >>> CNN()
      None
    """
    def __init__(self):
        """
        Initialize the Convolutional Neural Network (CNN) model.
        """
        super(CNN, self).__init__()

        self.classifier = nn.Sequential(
            nn.Conv2d(1024, 32, kernel_size=(7, 1), padding=(3, 0)),  # 7x32
            nn.ReLU(),
            nn.Dropout(0.0),
            nn.Conv2d(32, 20, kernel_size=(7, 1), padding=(3, 0)),
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """
        Perform forward pass through the CNN.

        Args:
            x (torch.Tensor): Input tensor of shape (batch_size, sequence_length, embedding_size).

        Returns:
            torch.Tensor: Output tensor of shape (batch_size, sequence_length, num_classes).

        L = protein length
        B = batch-size
        F = number of features (1024 for embeddings)
        N = number of classes (20 for 3Di)
        """

        # Permute input tensor to match expected shape
        # Input shape: (batch_size, sequence_length, embedding_size)
        # Output shape: (batch_size, embedding_size, sequence_length, 1)

        x = x.permute(0, 2, 1).unsqueeze(
            dim=-1
        )  # IN: X = (B x L x F); OUT: (B x F x L, 1)

        # Pass the input through the classifier
        # Output shape: (batch_size, num_classes, sequence_length, 1)
        Yhat = self.classifier(x)  # OUT: Yhat_consurf = (B x N x L x 1)

        # Remove the singleton dimension from the output tensor
        # Output shape: (batch_size, num_classes, sequence_length)
        Yhat = Yhat.squeeze(dim=-1)  # IN: (B x N x L x 1); OUT: ( B x L x N )
        return Yhat


def get_T5_model(
    model_dir: Path, model_name: str, cpu: bool, threads: int
) -> (T5EncoderModel, T5Tokenizer):
    """
    Loads a T5 model and tokenizer.

    Args:
        model_dir (Path): Directory where the model and tokenizer is be stored.
        model_name (str): Name of the pre-trained T5 model.
        cpu (bool): Whether to use CPU only.
        threads (int): Number of cpu threads.

    Returns:
        Tuple[T5EncoderModel, T5Tokenizer]: Tuple containing the loaded T5 model and tokenizer.
    """

    # sets the device

    # Torch load will map back to device from state, which often is GPU:0.
    # to overcome, need to explicitly map to active device

    global device

    torch.set_num_threads(threads)

    if cpu is True:
        device = torch.device("cpu")
        dev_name = "cpu"
    else:
        # check for NVIDIA/cuda
        if torch.cuda.is_available():
            device = torch.device("cuda:0")
            dev_name = "cuda:0"
        # check for apple silicon/metal
        elif torch.backends.mps.is_available():
            device = torch.device("mps")
            dev_name = "mps"
        else:
            device = torch.device("cpu")
            dev_name = "cpu"
            if cpu is not True:
                logger.warning(
                    "No available GPU was found, but --cpu was not specified"
                )
                logger.warning("ProstT5 will be run with CPU only")

    # logger device only if the function is called
    logger.info("Using device: {}".format(dev_name))

    # make dir if doesnt exist
    Path(model_dir).mkdir(parents=True, exist_ok=True)

    # load
    logger.info(f"Loading T5 from: {model_dir}/{model_name}")
    logger.info(f"If {model_dir}/{model_name} is not found, it will be downloaded")

    # check ProstT5 is downloaded
    # flag assumes transformers takes from local file (see #44)
    localfile = True
    download = False

    download = check_prostT5_download(model_dir, model_name)
    if download:
        localfile = False
        logger.info("ProstT5 not found. Downloading ProstT5 from Hugging Face")
    try:
        model = T5EncoderModel.from_pretrained(
            model_name,
            cache_dir=f"{model_dir}/",
            force_download=download,
            local_files_only=localfile,
        ).to(device)
        
    except:
        logger.warning("Download from Hugging Face failed. Trying backup from Zenodo.")
        logdir = f"{model_dir}/logdir"
        download_zenodo_prostT5(model_dir, logdir, threads )

        model = T5EncoderModel.from_pretrained(
            model_name,
            cache_dir=f"{model_dir}/",
            force_download=False,
            local_files_only=True,
        ).to(device)

    model = model.eval()
    vocab = T5Tokenizer.from_pretrained(
        model_name, cache_dir=f"{model_dir}/", do_lower_case=False
    )

    logger.info(f"{model_name} loaded")

    return model, vocab


def write_embeddings(
    embeddings: Dict[str, Dict[str, Tuple[List[str], Any, Any]]],
    out_path: Path,
) -> None:
    """
    Write embeddings to an output file.

    Args:
        embeddings (Dict[str, Dict[str, Tuple[List[str], Any, Any]]]): Predictions dictionary containing contig IDs, sequence IDs, predictions, and additional information.
        out_path (Path): Path to the output file.

    Returns:
        None
    """

    with h5py.File(str(out_path), "w") as hf:
        for sequence_id, embedding in embeddings.items():
            hf.create_dataset(sequence_id, data=embedding)


def write_predictions(
    hypotheticals: Dict,
    predictions: Dict[str, Tuple[List[str], Any, Any]],
    out_path: Path,
    mask_threshold: float,
    has_duplicate_locus: bool
) -> None:
    """
    Write predictions to an output file.

    Args:
        hypotheticals Dict: Hypothetical protein feature dictionary from Bakta
        predictions (Dict[str,  Tuple[List[str], Any, Any]]): Predictions dictionary containing sequence IDs, predictions, and additional information.
        out_path (Path): Path to the output file.
        mask_threshold (float): between 0 and 100 - below this ProstT5 confidence, 3Di predictions are masked


    Returns:
        None
    """
    ss_mapping = {
        0: "A",
        1: "C",
        2: "D",
        3: "E",
        4: "F",
        5: "G",
        6: "H",
        7: "I",
        8: "K",
        9: "L",
        10: "M",
        11: "N",
        12: "P",
        13: "Q",
        14: "R",
        15: "S",
        16: "T",
        17: "V",
        18: "W",
        19: "Y",
        20: "X" # fully mask the low confidence 3Di residues with X not lower case (not working for Foldseek v10, but X does)
    }

    mask_prop_threshold = mask_threshold/100

    # Filter out entries where the length of the value is 0
    # Issue #47
    predictions = {
        k: v for k, v in predictions.items() if len(v[0]) > 0
    }


    # mask low confidence residues
    for seq_id, (pred, mean_prob, all_prob) in predictions.items():

        # mask to X
        for i in range(len(pred)):
            if all_prob[0][i] < mask_prop_threshold:
                pred[i] = 20

    with open(out_path, "w+") as out_f:
        for feat in hypotheticals:
            if has_duplicate_locus:
                seq_id = feat["id"]
            else:
                seq_id = feat["locus"]
            pred = predictions.get(seq_id)  # predictions = {seq_id: (yhats, _, _)} or None

            if pred is not None:
                yhats = pred[0]  
                threedi_seq = "".join(ss_mapping[int(yhat)] for yhat in yhats)
                feat["3di"] = threedi_seq  # update the feature dictionary
                out_f.write(f">{seq_id}\n{threedi_seq}\n")
            else:
                feat["3di"] = None  # missing prediction

    logger.info(f"Finished writing results to {out_path}")
    return None


def write_probs(
    predictions: Dict[str, Tuple[int, float, Union[int, np.ndarray]]],
    output_path_mean: Path,
    output_path_all: Path,
    original_keys: list[str],
) -> None:
    """
    Write all ProstT5 encoder + CNN probabilities and mean probabilities to output files.

    Args:
        predictions (Dict[str, Tuple[int, float, Union[int, np.ndarray]]]):
        Predictions dictionary containing  sequence IDs, probabilities, and additional information.
        output_path_mean (str): Path to the output file for mean probabilities.
        output_path_all (str): Path to the output file for all probabilities.

    Returns:
        None
    """

    with open(output_path_mean, "w+") as out_f:
        for seq_id in original_keys:
            if seq_id not in predictions:
                logger.warning(f"Missing ProstT5 mean confidence for {seq_id}")
                continue
            _, mean_prob, _ = predictions[seq_id]
            out_f.write(f"{seq_id},{mean_prob}\n")

    if output_path_all is not None:
        with open(output_path_all, "w+") as out_f:
            for seq_id in original_keys:
                if seq_id not in predictions:
                    logger.warning(f"Missing ProstT5 confidence for {seq_id}")
                    continue

                _, _, all_probs = predictions[seq_id]

                # convert to percentage
                all_probs = all_probs * 100

                # flatten
                if isinstance(all_probs, np.ndarray):
                    all_probs_list = all_probs.flatten().tolist()
                else:
                    all_probs_list = all_probs

                rounded_list = [round(num, 2) for num in all_probs_list]

                json_data = json.dumps(
                    {"seq_id": seq_id, "probability": rounded_list}
                )
                out_f.write(json_data + "\n")


def toCPU(tensor: torch.Tensor) -> np.ndarray:
    """
    Move a tensor to CPU and convert it to a NumPy array.

    Args:
        tensor (torch.Tensor): Input tensor.

    Returns:
        np.ndarray: NumPy array.
    """
    if len(tensor.shape) > 1:
        return tensor.detach().cpu().squeeze(dim=-1).numpy()
    else:
        return tensor.detach().cpu().numpy()


def load_predictor(checkpoint_path: Union[str, Path]) -> CNN:
    """
    Load a pre-trained CNN ProstT5 prediction head weights from a checkpoint file.

    Args:
        checkpoint_path (Union[str, Path]): Path to the checkpoint file.

    Returns:
        CNN: Loaded CNN model.
    """

    model = CNN()

    # checkpoint_path = Path(CNN_DIR) / "cnn_chkpnt" / "model.pt"

    state = torch.load(checkpoint_path, map_location=device)

    # regular ProstT5 CNN 
    if checkpoint_path.suffix == '.pt':
        model.load_state_dict(state["state_dict"])
    # finetuned
    else:
        model.load_state_dict(state)


    model = model.eval()
    model = model.to(device)

    return model


def get_embeddings(
    hypotheticals: Dict[str, Tuple[str, ...]],
    cds_dict: Dict[str, Tuple[str, ...]],
    out_path: Path,
    prefix: str,
    model_dir: Path,
    model_name: str,
    checkpoint_path: Path,
    output_3di: Path,
    output_h5_per_residue: Path,
    output_h5_per_protein: Path,
    half_precision: bool,
    max_residues: int = 100000,
    max_seq_len: int = 30000,
    max_batch: int = 10000,
    cpu: bool = False,
    output_probs: bool = True,
    save_per_residue_embeddings: bool = False,
    save_per_protein_embeddings: bool = False,
    threads: int = 1,
    mask_threshold: float = 0,
    has_duplicate_locus: bool = False
) -> bool:
    """
    Generate embeddings and predictions for protein sequences using ProstT5 encoder & CNN prediction head.

    Args:
        hypotheticals ( Dict[str, Tuple[str, ...]]):  dictionary containing CDS IDs feature information
        cds_dict ( Dict[str, Tuple[str, ...]]):  dictionary containing CDS IDs and corresponding protein sequences.
        out_path (Path): Path to the output directory.
        prefix (str): Prefix for the output files.
        model_dir (Path): Directory containing the pre-trained model.
        model_name (str): Name of the pre-trained model.
        output_3di (Path): Path to the output 3Di file.
        output_h5_per_residue (Path): Path to the output h5 per residue embeddings file.
        output_h5_per_protein (Path): Path to the output h5 per proteins embeddings file.
        half_precision (bool): Whether to use half precision for the models.
        max_residues (int, optional): Maximum number of residues allowed in a batch. Defaults to 3000.
        max_seq_len (int, optional): Maximum sequence length allowed. Defaults to 1000.
        max_batch (int, optional): Maximum batch size. Defaults to 100.
        cpu (bool, optional): Whether to use CPU for processing. Defaults to False.
        output_probs (bool, optional): Whether to output probabilities. Defaults to True.
        save_embeddings (bool, optional): Whether to save embeddings to h5 file. Defaults to False. Will  save per residue embeddings
        per_protein_embeddings (bool, optional): Whether to save per protein mean embeddings to h5 file. Defaults to False.
        threads (int): number of cpu threads
        mask_threshold (float) : 0-100 - below this ProstT5 confidence threshold, these residues are masked
        has_duplicate_locus (bool) : some euks have dupe locus tags

    Returns:
        bool: True if embeddings and predictions are generated successfully.
    """

    predictions = {}
    batch_predictions = {}

    if save_per_residue_embeddings:
        embeddings_per_residue = {}
        batch_embeddings_per_residue = {}
    if save_per_protein_embeddings:
        embeddings_per_protein = {}
        batch_embeddings_per_protein = {}


    prostt5_prefix = "<AA2fold>"


    model, vocab = get_T5_model(model_dir, model_name, cpu, threads)
    predictor = load_predictor(checkpoint_path)

    logger.info("Beginning ProstT5 predictions")

    if half_precision:
        model = model.half()
        predictor = predictor.half()
        logger.info("Using models in half-precision")
    else:
        logger.info("Using models in full-precision")

    # 3Di predictions

    fail_ids = []

    for k, v in list(cds_dict.items()):
        if len(v) == 0:
            logger.info(f"Skipping empty CDS entry as it has no amino acid string associated (likely pseudo): key={k}")
            del cds_dict[k]

    # sort sequences by length
    seq_dict = []
    fail_ids = []

    for k, seq in cds_dict.items():
        if isinstance(seq, str) and seq:
            clean_seq = (
                seq.replace("U", "X")
                .replace("Z", "X")
                .replace("O", "X")
            )
            seq_dict.append((k, clean_seq, len(clean_seq))) # in the correct format for the remainder of the code
        else:
            logger.warning(
                f"Protein header {k} is corrupt. It will be saved in fails.tsv"
            )
            fail_ids.append(k)


    original_keys = list(cds_dict.keys())
        # --- sort once ---
    seq_dict.sort(key=lambda x: x[2], reverse=True)

    batch = list()
    for seq_idx, (pdb_id, seq, slen) in enumerate(tqdm(seq_dict, desc=f"Predicting 3Di"), 1):

        # replace non-standard AAs
        seq = seq.replace("U", "X").replace("Z", "X").replace("O", "X")
        seq_len = len(seq)
        seq = prostt5_prefix + " " + " ".join(list(seq))
        batch.append((pdb_id, seq, seq_len))

        # count residues in current batch and add the last sequence length to
        # avoid that batches with (n_res_batch > max_residues) get processed
        n_res_batch = sum([s_len for _, _, s_len in batch]) + seq_len
        if (
            len(batch) >= max_batch
            or n_res_batch >= max_residues
            or seq_idx == len(seq_dict)
            or seq_len > max_seq_len
        ):
            pdb_ids, seqs, seq_lens = zip(*batch)
            batch = list()

            token_encoding = vocab.batch_encode_plus(
                seqs,
                add_special_tokens=True,
                padding="longest",
                return_tensors="pt",
            ).to(device)
            try:
                with torch.no_grad():
                    embedding_repr = model(
                        token_encoding.input_ids,
                        attention_mask=token_encoding.attention_mask,
                    )
            except RuntimeError:
                logger.warning(f" number of residues in batch {n_res_batch}")
                logger.warning(f" seq length is {seq_len}")
                logger.warning(f" ids are {pdb_ids}")
                logger.warning(
                    "RuntimeError during embedding for {} (L={})".format(
                        pdb_id, seq_len
                    )
                )
                for id in pdb_ids:
                    fail_ids.append(id)
                continue


            # ProtT5 appends a special tokens at the end of each sequence
            # Mask this also out during inference while taking into account the prostt5 prefix
            try:
                for idx, s_len in enumerate(seq_lens):
                    token_encoding.attention_mask[idx, s_len + 1] = 0

                # extract last hidden states (=embeddings)
                residue_embedding = embedding_repr.last_hidden_state.detach()
                # mask out padded elements in the attention output (can be non-zero) for further processing/prediction
                residue_embedding = (
                    residue_embedding
                    * token_encoding.attention_mask.unsqueeze(dim=-1)
                )
                # slice off embedding of special token prepended before to each sequence
                residue_embedding = residue_embedding[:, 1:]
                prediction = predictor(residue_embedding)


                # compute max probabilities per token/residue
                probabilities = toCPU(
                    torch.max(
                        F.softmax(prediction, dim=1), dim=1, keepdim=True
                    )[0]
                )

                prediction = toCPU(
                    torch.max(prediction, dim=1, keepdim=True)[1]
                ).astype(np.byte)


                    # to get logits
                # t = prediction.transpose(1, 2)  # changes ( B x L x N ) to ( B x N x L ) 
                # logits = t.detach().cpu().numpy()
                # print(logits[0])
                # print(logits[0].shape)

                    # prob_tensor = F.softmax(prediction, dim=1).cpu()
                    # #print(prob_tensor.shape)
                    # prob_matrix = prob_tensor.squeeze(0).transpose(0, 1)
                    # #print(prob_matrix.shape)
                    # all_sampled_states = {}
                    # i = 1
                    # while i < 11:
                    #     all_sampled_states[i] = torch.multinomial(prob_matrix, 1).squeeze(1)
                    #     i += 1

                    #sampled_states = torch.multinomial(prob_matrix, 1).squeeze(1)
                    #sampled_states = torch.multinomial(prob_matrix, 100, replacement=True)
                    #print(sampled_states)
                    #print(sampled_states.shape)

                # batch-size x seq_len x embedding_dim
                # extra token is added at the end of the seq
                for batch_idx, identifier in enumerate(pdb_ids):
                    s_len = seq_lens[batch_idx]

                    # save embeddings
                    if save_per_residue_embeddings or save_per_protein_embeddings:
                        try:
                            # account for prefix in offset
                            emb = embedding_repr.last_hidden_state[
                                batch_idx, 1 : s_len + 1
                            ]

                            if save_per_residue_embeddings:
                                batch_embeddings_per_residue[identifier] = (
                                    emb.detach().cpu().numpy().squeeze()
                                )

                            if save_per_protein_embeddings:
                                batch_embeddings_per_protein[identifier] = (
                                    emb.mean(dim=0).detach().cpu().numpy().squeeze()
                                )

                        except:
                            logger.warning(
                                f"Saving embeddings failed for {identifier}"
                            )

                    # slice off padding and special token appended to the end of the sequence
                    pred = prediction[batch_idx, :, 0:s_len].squeeze()

                    # always return the mean probs
                    mean_prob = round(
                            100 * np.mean(probabilities[batch_idx, :, 0:s_len]), 2
                        )

                    if output_probs:  # if you want the per-residue probs
                        all_prob = probabilities[batch_idx, :, 0:s_len]
                        batch_predictions[identifier] = (
                            pred,
                            mean_prob,
                            all_prob,
                        )
                    else:
                        batch_predictions[identifier] = (pred, mean_prob, None)

                    try:
                        len(batch_predictions[identifier][0])
                    except:
                        logger.warning(
                            f"{identifier} prediction has length 0"
                        )
                        fail_ids.append(identifier)
                        continue

                    if s_len != len(batch_predictions[identifier][0]):
                        logger.warning(
                            f"Length mismatch for {identifier}: is:{len(batch_predictions[identifier][0])} vs should:{s_len}"
                        )

                    for k in original_keys:
                        if k in batch_predictions:
                            predictions[k] = batch_predictions[k]
                
                    if save_per_residue_embeddings:
                        for k in original_keys:
                            if k in batch_predictions:
                                embeddings_per_residue[k] = batch_embeddings_per_residue[k]

                    if save_per_protein_embeddings:
                        for k in original_keys:
                            if k in batch_predictions:
                                embeddings_per_protein[k] = batch_embeddings_per_protein[k]

            except IndexError:
                logger.warning(
                    "Index error during prediction for {} (L={})".format(
                        pdb_id, seq_len
                    )
                )
                for id in pdb_ids:
                    fail_ids.append(id)
                continue

    # write list of fails if length > 0
    if len(fail_ids) > 0:
        fail_tsv: Path = Path(out_path) / "fails.tsv"

        # Convert the list to a list of lists
        data_as_list_of_lists = [[str(item)] for item in fail_ids]

        # Write the list to a TSV file
        with open(fail_tsv, "w", newline="") as file:
            tsv_writer = csv.writer(file, delimiter="\t")
            tsv_writer.writerows(data_as_list_of_lists)

    write_predictions(hypotheticals, predictions, output_3di,  mask_threshold, has_duplicate_locus)

    if save_per_residue_embeddings:
        write_embeddings(embeddings_per_residue, output_h5_per_residue)

    if save_per_protein_embeddings:
        write_embeddings(embeddings_per_protein, output_h5_per_protein)

    # always write the mean embeddings
    mean_probs_out_path: Path = (
        Path(out_path) / f"{prefix}_prostT5_3di_mean_probabilities.csv"
    )

    # output per residue probs
    if output_probs:
        all_probs_out_path: Path = (
            Path(out_path) / f"{prefix}_prostT5_3di_all_probabilities.json"
        )
    else:
        all_probs_out_path = None

    write_probs(predictions, mean_probs_out_path, all_probs_out_path, original_keys)

    return predictions
