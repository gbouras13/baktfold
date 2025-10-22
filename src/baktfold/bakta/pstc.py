import logging
import pandas as pd
# import sqlite3

# from concurrent.futures import ThreadPoolExecutor
from typing import Sequence, Tuple

# import baktfold.bakta.config as cfg
import baktfold.bakta.constants as bc
from loguru import logger
from pathlib import Path

import csv


def parse(features: Sequence[dict], foldseek_df: pd.DataFrame, db_name: str = 'swissprot') -> None:
    """Update CDS in place with PSTC hits from foldseek_df if they pass filters."""

    # Convert foldseek_df to a lookup table keyed by query ID
    foldseek_hits = {row['query']: row for _, row in foldseek_df.iterrows()}

    updated_count = 0
    for cds in features:
        aa_identifier = cds.get('locus')
        if aa_identifier not in foldseek_hits:
            continue  # no hit, leave cds as-is

        row = foldseek_hits[aa_identifier]
        query_cov = float(row['qCov'])
        subject_cov = float(row['tCov'])
        identity = float(row['fident'])
        evalue = float(row['evalue'])
        bitscore = float(row['bitscore'])
        target_id = row['target']

        # keep only the accession

        if db_name == "swissprot" or db_name == "afdb":
            accession = target_id.split('-')[1]
        elif db_name == "pdb":
            accession = target_id.split('-')[0]
        else: # cath here as well as custom
            accession = target_id

        if (
            query_cov >= bc.MIN_PSTC_QCOVERAGE
            and subject_cov >= bc.MIN_PSTC_TCOVERAGE
            and identity >= bc.MIN_PSTC_IDENTITY
        ):
            new_pstc = {
                'source': db_name,
                'id': accession,
                'query_cov': query_cov,
                'subject_cov': subject_cov,
                'identity': identity,
                'score': bitscore,
                'evalue': evalue,
            }

            # Append if 'pstc' exists, else create
            if 'pstc' in cds:
                # ensure it's a list
                if isinstance(cds['pstc'], dict):
                    cds['pstc'] = [cds['pstc'], new_pstc]
                elif isinstance(cds['pstc'], list):
                    cds['pstc'].append(new_pstc)
                else:
                    # fallback: overwrite if type is unexpected
                    cds['pstc'] = [new_pstc]
            else:
                cds['pstc'] = new_pstc
        updated_count += 1

    # this gets the running tally for any DB
    # logger.info(f"PSTC for {db_name} updated in place for {sum('pstc' in cds for cds in features)} CDSs")
    logger.info(f"PSTC for {db_name} updated in place for {updated_count} CDSs")

    return features



def lookup(features: Sequence[dict], baktfold_db: Path, custom_annotations: Path):
    """Lookup PSTC information in swissprot"""
    no_pscc_lookups = 0

    # simple dictionary of accessions and protein_name
    swissprot_dict = {}
    with open(f"{baktfold_db}/swissprot.tsv", "r") as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            if len(row) >= 2:
                swissprot_dict[row[0]] = row[1]

    afdb_dict = {}
    with open(f"{baktfold_db}/AFDBClusters.tsv", "r") as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            if len(row) >= 2:
                afdb_dict[row[0]] = row[1]

    pdb_dict = {}
    with open(f"{baktfold_db}/pdb.tsv", "r") as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            if len(row) >= 2:
                pdb_dict[row[0]] = row[1]

    cath_dict = {}
    with open(f"{baktfold_db}/cath.tsv", "r") as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            if len(row) >= 3:
                pdb_dict[row[0]] = row[2] # 3 columns - the second is the CATH code

    # custom
    if custom_annotations:
        custom_dict = {}
        with open(f"{custom_annotations}", "r") as f:
            reader = csv.reader(f, delimiter="\t")
            for row in reader:
                if len(row) >= 2:
                    custom_dict[row[0]] = row[1]

    for feat in features:
        pstc = feat.get('pstc')
        if not pstc:
            continue

        # Normalize to list for consistent handling
        pstc_entries = pstc if isinstance(pstc, list) else [pstc]

        for entry in pstc_entries:
            accession = entry.get('id')
            source = entry.get('source')
            if source == 'swissprot' and accession in swissprot_dict:
                entry['description'] = swissprot_dict[accession]
            elif source == 'afdb' and accession in afdb_dict:
                entry['description'] = afdb_dict[accession]
            elif source == 'pdb' and accession in pdb_dict:
                entry['description'] = pdb_dict[accession]
            elif source == 'cath' and accession in cath_dict:
                entry['description'] = cath_dict[accession]
            elif source == 'custom_db':
                if accession in custom_dict:
                    entry['description'] = custom_dict[accession]
                else:
                    entry['description'] = accession # mark as accession if no annotation given for custom for now
            else:
                # Keep "hypothetical protein" for missing
                entry['description'] = "hypothetical protein"

        # Write back normalized list or single entry
        feat['pstc'] = pstc_entries if isinstance(pstc, list) else pstc_entries[0]

    return features





# def lookup(features: Sequence[dict], pseudo: bool = False):
#     """Lookup PSCC information"""
#     no_pscc_lookups = 0
#     try:
#         rec_futures = []
#         with sqlite3.connect(f"file:{cfg.db_path.joinpath('bakta.db')}?mode=ro&nolock=1&cache=shared", uri=True, check_same_thread=False) as conn:
#             conn.execute('PRAGMA omit_readlock;')
#             conn.row_factory = sqlite3.Row
#             with ThreadPoolExecutor(max_workers=max(10, cfg.threads)) as tpe:  # use min 10 threads for IO bound non-CPU lookups
#                 for feature in features:
#                     uniref50_id = None
#                     if(pseudo):  # if pseudogene use pseudogene info
#                         if('psc' in feature[bc.PSEUDOGENE]):
#                             uniref50_id = feature[bc.PSEUDOGENE]['psc'].get(DB_PSCC_COL_UNIREF50, None)
#                     else:
#                         if('psc' in feature):
#                             uniref50_id = feature['psc'].get(DB_PSCC_COL_UNIREF50, None)
#                         elif('pscc' in feature):
#                             uniref50_id = feature['pscc'].get(DB_PSCC_COL_UNIREF50, None)
#                     if(uniref50_id is not None):
#                         if(bc.DB_PREFIX_UNIREF_50 in uniref50_id):
#                             uniref50_id = uniref50_id[9:]  # remove 'UniRef50_' prefix
#                         future = tpe.submit(fetch_db_pscc_result, conn, uniref50_id)
#                         rec_futures.append((feature, future))

#         for (feature, future) in rec_futures:
#             rec = future.result()
#             if(rec is not None):
#                 pscc = parse_annotation(rec)
#                 if(pseudo):
#                     feature[bc.PSEUDOGENE]['pscc'] = pscc
#                 else:
#                     if('pscc' in feature):
#                         feature['pscc'] = {**feature['pscc'], **pscc}  # merge dicts, add PSCC annotation info to PSCC alignment info
#                     else:
#                         feature['pscc'] = pscc  # add PSCC annotation info
#                 no_pscc_lookups += 1
#                 log.debug(
#                     'lookup: seq=%s, start=%i, stop=%i, strand=%s, UniRef50=%s, product=%s',
#                     feature['sequence'], feature['start'], feature['stop'], feature['strand'], pscc.get(DB_PSCC_COL_UNIREF50, ''), pscc.get(DB_PSCC_COL_PRODUCT, '')
#                 )
#             else:
#                 log.debug('lookup: ID not found! uniref50_id=%s', uniref50_id)
#     except Exception as ex:
#         log.exception('Could not read PSCCs from db!', ex)
#         raise Exception('SQL error!', ex)
#     log.info('looked-up=%i', no_pscc_lookups)


# def fetch_db_pscc_result(conn: sqlite3.Connection, uniref50_id: str):
#     c = conn.cursor()
#     c.execute('select * from pscc where uniref50_id=?', (uniref50_id,))
#     rec = c.fetchone()
#     c.close()
#     return rec


# def parse_annotation(rec) -> dict:
#     uniref_full_id = bc.DB_PREFIX_UNIREF_50 + rec[DB_PSCC_COL_UNIREF50]
#     pscc = {
#         DB_PSCC_COL_UNIREF50: uniref_full_id,  # must not be NULL/None
#         'db_xrefs': [
#             'SO:0001217',
#             f'{bc.DB_XREF_UNIREF}:{uniref_full_id}'
#         ]
#     }
#     # add non-empty PSCC annotations and attach database prefixes to identifiers
#     if(rec[DB_PSCC_COL_PRODUCT]):
#         pscc[DB_PSCC_COL_PRODUCT] = rec[DB_PSCC_COL_PRODUCT]
#     return pscc