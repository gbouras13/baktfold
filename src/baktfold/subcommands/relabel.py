#!/usr/bin/env python3
"""``baktfold json --relabel``: rerun label selection on the hits stored in a JSON.

The Foldseek hits (``pstc``) in a Baktfold JSON are kept as they are; only the
choice of label is redone, with the same code ``run``/``compare`` use
(annotation.combine_annotation). This brings outputs written by earlier
versions up to the current label rules without searching again.
"""

from pathlib import Path
from typing import Optional, Sequence

from loguru import logger

import baktfold.bakta.annotation as anno
import baktfold.bakta.pstc as pstc
from baktfold.utils.util import get_version


def relabel_features(
    features: Sequence[dict],
    fast: bool,
    has_duplicate_locus: bool,
    database: Optional[Path] = None,
    pdb_tophit: Optional[Path] = None,
    prostt5_probs: Optional[Path] = None,
) -> dict:
    """Rerun label selection in place on every CDS that has stored hits.

    Args:
      features: CDS/protein features (each distinct feature once).
      fast: the run skipped AFDB (from the JSON provenance).
      has_duplicate_locus: hits are keyed by feature ``id`` rather than ``locus``.
      database: Baktfold database directory; if given, descriptions are looked
        up again (chain-level PDB names need a database with ``pdb_chain``).
      pdb_tophit: ``<prefix>_pdb_tophit.tsv`` of the original run, for PDB
        chains that older JSONs don't record.
      prostt5_probs: ``<prefix>_prostT5_3di_mean_probabilities.csv`` of the
        original run, for JSONs that don't carry ``prostt5_confidence``.

    Returns:
      dict: provenance and counts, stored as ``baktfold_relabel`` in the JSON.
    """
    with_hits = [f for f in features if f.get('pstc')]
    before = {id(f): f.get('product') for f in with_hits}
    structures = any(
        isinstance(p, dict) and 'tmscore' in p
        for f in with_hits
        for p in (f['pstc'] if isinstance(f['pstc'], list) else [f['pstc']])
    )

    # 1. undo the earlier labelling (needs the hits as they were written)
    n_reset = sum(anno.reset_baktfold_annotation(f) for f in with_hits)

    # 2. optionally refresh descriptions / chains from the database
    refresh = {'pdb_chains_added': 0, 'pdb_no_chain_kept': 0}
    if database is not None:
        chains = pstc.read_pdb_tophit_chains(pdb_tophit) if pdb_tophit else None
        refresh = pstc.refresh_descriptions(with_hits, database, chains, has_duplicate_locus)
        if refresh['pdb_no_chain_kept']:
            logger.warning(
                f"{refresh['pdb_no_chain_kept']} PDB hits have no recorded chain and keep their "
                "PDB entry title. Pass --pdb-tophit <prefix>_pdb_tophit.tsv from the original "
                "run for chain-level names."
            )
    elif pdb_tophit is not None:
        logger.warning("--pdb-tophit has no effect without -d/--database; it is ignored.")

    # 3. ProstT5 confidence, needed for the confidence tier
    if prostt5_probs is not None:
        anno.attach_prostt5_confidence(with_hits, prostt5_probs, has_duplicate_locus)

    # 4. label selection, exactly as run/compare do it
    for f in with_hits:
        anno.combine_annotation(f, fast, structures=structures)

    # A ProstT5 run without stored ProstT5 confidence would be tiered with the
    # structure-input thresholds; leave the tier out rather than mislabel it.
    no_conf = 0
    if not structures:
        for f in with_hits:
            if 'annotation_confidence' in f and f.get('prostt5_confidence') is None:
                f.pop('annotation_confidence')
                no_conf += 1
        if no_conf:
            logger.warning(
                f"{no_conf} CDS have no ProstT5 confidence in the JSON, so no Annotation_Confidence "
                "is written for them. Pass --prostt5-probs <prefix>_prostT5_3di_mean_probabilities.csv "
                "from the original run to restore it."
            )

    changed = sum(1 for f in with_hits if f.get('product') != before[id(f)])
    functional = sum(1 for f in with_hits if 'hypothetical' not in f)
    logger.info(
        f"Relabelled {len(with_hits)} CDS with stored hits ({n_reset} previously labelled by Baktfold): "
        f"{changed} labels changed, {functional} now carry a functional label"
    )
    return {
        'baktfold': get_version(),
        'database': str(database) if database is not None else None,
        'pdb_tophit': str(pdb_tophit) if pdb_tophit is not None else None,
        'prostt5_probs': str(prostt5_probs) if prostt5_probs is not None else None,
        'cds_with_hits': len(with_hits),
        'labels_changed': changed,
        'functional': functional,
        'pdb_chains_added': refresh['pdb_chains_added'],
        'pdb_no_chain_kept': refresh['pdb_no_chain_kept'],
        'confidence_dropped': no_conf,
    }
