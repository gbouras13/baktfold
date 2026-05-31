#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
3Di prediction for baktfold — wraps pholdlib's shared inference engine.

Baktfold-specific: flat cds_dict (no contig nesting), Bakta hypotheticals
format with in-place annotation updates, has_duplicate_locus support.

Code adapted from @mheinzinger
https://github.com/mheinzinger/ProstT5/blob/main/scripts/predict_3Di_encoderOnly.py
"""

from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import h5py
from loguru import logger

# ── pholdlib shared components ────────────────────────────────────────────────
from pholdlib.prostt5.model import CNN, get_T5_model, load_predictor, toCPU  # noqa: F401
from pholdlib.prostt5.inference import run_prostt5_inference
from pholdlib.prostt5.output import (
    SS_MAPPING,
    write_fail_ids,
    write_probs,
)

# ── baktfold-specific DB helpers ──────────────────────────────────────────────
from baktfold.databases.db import check_prostT5_download, download_zenodo_prostT5


# ─────────────────────────────────────────────────────────────────────────────
# HDF5 embedding writer (flat seq_id structure — no contig nesting in baktfold)
# ─────────────────────────────────────────────────────────────────────────────

def write_embeddings(embeddings: Dict[str, Any], out_path: Path) -> None:
    """Write per-residue or per-protein embeddings to HDF5 (flat key structure)."""
    with h5py.File(str(out_path), "w") as hf:
        for sequence_id, embedding in embeddings.items():
            hf.create_dataset(sequence_id, data=embedding)


# ─────────────────────────────────────────────────────────────────────────────
# 3Di FASTA + Bakta annotation writer (baktfold-specific)
# ─────────────────────────────────────────────────────────────────────────────

def write_predictions(
    hypotheticals: List[Dict],
    predictions: Dict[str, Tuple],
    out_path: Path,
    mask_threshold: float,
    has_duplicate_locus: bool = False,
) -> None:
    """Write 3Di predictions to FASTA and update Bakta hypotheticals in-place.

    Args:
        hypotheticals: List of Bakta feature dicts. Each is mutated in-place
                       with a ``"3di"`` key set to the predicted 3Di string
                       (or None if prediction failed / was skipped).
        predictions: Flat ``{seq_id: (pred, mean_prob, all_prob)}`` dict.
        out_path: Output FASTA path.
        mask_threshold: Residues with max softmax prob (0–100) below this
                        threshold are replaced with 'X'.
        has_duplicate_locus: If True, use ``feat["id"]`` as seq_id (needed for
                             eukaryotic inputs that may have duplicate locus tags).
                             Otherwise use ``feat["locus"]``.
    """
    mask_prop = mask_threshold / 100

    # drop zero-length predictions (issue #47)
    predictions = {k: v for k, v in predictions.items() if len(v[0]) > 0}

    # apply confidence masking in-place on pred index arrays
    for seq_id, (pred, mean_prob, all_prob) in predictions.items():
        for i in range(len(pred)):
            if all_prob[0][i] < mask_prop:
                pred[i] = 20  # 'X'

    with open(out_path, "w+") as out_f:
        for feat in hypotheticals:
            seq_id = feat["id"] if has_duplicate_locus else feat["locus"]
            pred_tuple = predictions.get(seq_id)
            if pred_tuple is not None:
                yhats = pred_tuple[0]
                threedi_seq = "".join(SS_MAPPING[int(y)] for y in yhats)
                feat["3di"] = threedi_seq  # mutate Bakta feature dict in-place
                out_f.write(f">{seq_id}\n{threedi_seq}\n")
            else:
                feat["3di"] = None  # no prediction (OOM / corrupt entry)

    logger.info(f"Finished writing 3Di FASTA to {out_path}")


# ─────────────────────────────────────────────────────────────────────────────
# Main entry point
# ─────────────────────────────────────────────────────────────────────────────

def get_embeddings(
    hypotheticals: List[Dict],
    cds_dict: Dict[str, str],
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
    has_duplicate_locus: bool = False,
) -> Dict:
    """Run ProstT5 + CNN 3Di prediction for all sequences in *cds_dict*.

    Args:
        hypotheticals: List of Bakta feature dicts (mutated in-place with "3di").
        cds_dict: Flat ``{seq_id: amino_acid_str}`` dict.
        out_path: Directory for output files.
        prefix: Filename prefix for CSV / JSONL outputs.
        model_dir: Directory where ProstT5 is cached.
        model_name: HuggingFace model identifier.
        checkpoint_path: Path to the CNN ``.pt`` checkpoint.
        output_3di: Output FASTA path for 3Di sequences.
        output_h5_per_residue: HDF5 path for per-residue embeddings.
        output_h5_per_protein: HDF5 path for per-protein embeddings.
        half_precision: If True, cast model + predictor to fp16 after loading.
        max_residues: Max total residues per inference batch.
        max_seq_len: Sequences longer than this flush a batch immediately.
        max_batch: Max sequences per batch.
        cpu: Force CPU inference.
        output_probs: Whether to write per-residue probability JSONL.
        save_per_residue_embeddings: Save per-residue HDF5.
        save_per_protein_embeddings: Save per-protein HDF5.
        threads: Number of CPU threads for torch.
        mask_threshold: Residues with max softmax prob < threshold/100 → 'X'.
        has_duplicate_locus: If True use feat["id"] rather than feat["locus"].

    Returns:
        predictions: Flat ``{seq_id: (pred, mean_prob, all_prob)}`` dict,
                     in original cds_dict key order.
    """
    # ── load model ──────────────────────────────────────────────────────────
    model, vocab, device = get_T5_model(
        model_dir, model_name, cpu, threads,
        check_fn=check_prostT5_download,
        zenodo_fn=download_zenodo_prostT5,
    )
    predictor = load_predictor(checkpoint_path, device)

    logger.info("Beginning ProstT5 predictions")

    if half_precision:
        model = model.half()
        predictor = predictor.half()
        logger.info("Using models in half-precision")
    else:
        logger.info("Using models in full-precision")

    # ── build seq_dict (skip empty / non-string entries) ────────────────────
    original_keys = list(cds_dict.keys())
    seq_dict: List[Tuple] = []
    fail_ids: List[str] = []

    for k, seq in cds_dict.items():
        if isinstance(seq, str) and seq:
            seq_dict.append((k, seq, len(seq)))
        else:
            logger.warning(
                f"Protein header {k} is corrupt or empty — will be saved in fails.tsv"
            )
            fail_ids.append(k)

    # sort descending by length (minimises padding in each batch)
    seq_dict.sort(key=lambda x: x[2], reverse=True)

    # ── run shared inference engine ──────────────────────────────────────────
    predictions, emb_res, emb_prot, inf_fail_ids = run_prostt5_inference(
        seq_dict,
        model, vocab, predictor, device,
        max_residues=max_residues,
        max_seq_len=max_seq_len,
        max_batch=max_batch,
        output_probs=output_probs,
        save_per_residue_embeddings=save_per_residue_embeddings,
        save_per_protein_embeddings=save_per_protein_embeddings,
        desc="Predicting 3Di",
    )
    fail_ids.extend(inf_fail_ids)

    # restore original key order
    predictions = {k: predictions[k] for k in original_keys if k in predictions}

    # ── write outputs ────────────────────────────────────────────────────────
    if fail_ids:
        write_fail_ids(fail_ids, Path(out_path) / "fails.tsv")

    write_predictions(
        hypotheticals, predictions, output_3di, mask_threshold, has_duplicate_locus
    )

    if save_per_residue_embeddings:
        write_embeddings(emb_res, output_h5_per_residue)
    if save_per_protein_embeddings:
        write_embeddings(emb_prot, output_h5_per_protein)

    mean_probs_path = Path(out_path) / f"{prefix}_prostT5_3di_mean_probabilities.csv"
    all_probs_path = (
        Path(out_path) / f"{prefix}_prostT5_3di_all_probabilities.json"
        if output_probs else None
    )
    write_probs(predictions, mean_probs_path, all_probs_path, original_keys)

    return predictions
