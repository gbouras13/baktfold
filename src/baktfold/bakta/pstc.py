import logging
import pandas as pd
import sqlite3

from concurrent.futures import ThreadPoolExecutor
from typing import Sequence, Tuple

# import baktfold.bakta.config as cfg
import baktfold.bakta.constants as bc
from loguru import logger
from pathlib import Path
from collections import defaultdict
import csv


def parse(features: Sequence[dict], foldseek_df: pd.DataFrame, db_name: str = 'swissprot') -> None:
    """Update CDS in place with PSTC hits from foldseek_df if they pass filters."""

    # Convert foldseek_df to a lookup table keyed by query ID
    foldseek_hits = {row['query']: row for _, row in foldseek_df.iterrows()}

    # each query maps to a list of rows now (to handle multiple CATH greedy tophits for multidomain proteins)
    foldseek_hits = defaultdict(list)
    for _, row in foldseek_df.iterrows():
        foldseek_hits[row['query']].append(row)

    updated_count = 0


    for cds in features:
        aa_identifier = cds.get('locus')
        if aa_identifier not in foldseek_hits:
            continue  # no hits, skip

        cds_updated = False  

        # Iterate over *all* hits for this query
        for row in foldseek_hits[aa_identifier]:
            query_cov = float(row['qCov'])
            subject_cov = float(row['tCov'])
            identity = float(row['fident'])
            evalue = float(row['evalue'])
            bitscore = float(row['bitscore'])
            target_id = row['target']

            # Extract accession depending on database
            if db_name in {"swissprot", "afdb"}:
                accession = target_id.split('-')[1]
            elif db_name == "pdb":
                accession = target_id.split('-')[0]
            else:  # cath and custom
                accession = target_id

            # Apply your filters
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

                # Append or initialize 'pstc'
                if 'pstc' in cds:
                    if isinstance(cds['pstc'], dict):
                        cds['pstc'] = [cds['pstc'], new_pstc]
                    elif isinstance(cds['pstc'], list):
                        cds['pstc'].append(new_pstc)
                    else:
                        cds['pstc'] = [new_pstc]
                else:
                    cds['pstc'] = [new_pstc]  # ← ensure list, since we may have many hits

                
                cds_updated = True  

        # Increment only once per CDS that had at least one valid hit (CATH might have multiple)
        if cds_updated:
            updated_count += 1

    logger.info(f"PSTC for {db_name} updated in place for {updated_count} CDSs")
    return features


def lookup_custom(features: Sequence[dict], baktfold_db: Path, custom_annotations: Path):
    """Lookup PSTC information from custom db """
    no_pstc_lookups = 0

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
            if source == 'custom_db':
                if accession in custom_dict:
                    entry['description'] = custom_dict[accession]
                else:
                    entry['description'] = accession # mark as accession if no annotation given for custom for now

        # Write back normalized list or single entry
        feat['pstc'] = pstc_entries if isinstance(pstc, list) else pstc_entries[0]

    return features


# def lookup(features: Sequence[dict], baktfold_db: Path, custom_annotations: Path):
#     """Lookup PSTC information"""
#     no_pscc_lookups = 0

#     # simple dictionary of accessions and protein_name
#     swissprot_dict = {}
#     with open(f"{baktfold_db}/swissprot.tsv", "r") as f:
#         reader = csv.reader(f, delimiter="\t")
#         for row in reader:
#             if len(row) >= 2:
#                 swissprot_dict[row[0]] = row[1]

#     afdb_dict = {}
#     with open(f"{baktfold_db}/AFDBClusters.tsv", "r") as f:
#         reader = csv.reader(f, delimiter="\t")
#         for row in reader:
#             if len(row) >= 2:
#                 afdb_dict[row[0]] = row[1]

#     pdb_dict = {}
#     with open(f"{baktfold_db}/pdb.tsv", "r") as f:
#         reader = csv.reader(f, delimiter="\t")
#         for row in reader:
#             if len(row) >= 2:
#                 pdb_dict[row[0]] = row[1]

#     cath_dict = {}
#     with open(f"{baktfold_db}/cath.tsv", "r") as f:
#         reader = csv.reader(f, delimiter="\t")
#         for row in reader:
#             if len(row) >= 3:
#                 pdb_dict[row[0]] = row[2] # 3 columns - the second is the CATH code

#     # custom
#     if custom_annotations:
#         custom_dict = {}
#         with open(f"{custom_annotations}", "r") as f:
#             reader = csv.reader(f, delimiter="\t")
#             for row in reader:
#                 if len(row) >= 2:
#                     custom_dict[row[0]] = row[1]

#     for feat in features:
#         pstc = feat.get('pstc')
#         if not pstc:
#             continue

#         # Normalize to list for consistent handling
#         pstc_entries = pstc if isinstance(pstc, list) else [pstc]

#         for entry in pstc_entries:
#             accession = entry.get('id')
#             source = entry.get('source')
#             if source == 'swissprot' and accession in swissprot_dict:
#                 entry['description'] = swissprot_dict[accession]
#             elif source == 'afdb' and accession in afdb_dict:
#                 entry['description'] = afdb_dict[accession]
#             elif source == 'pdb' and accession in pdb_dict:
#                 entry['description'] = pdb_dict[accession]
#             elif source == 'cath' and accession in cath_dict:
#                 entry['description'] = cath_dict[accession]
#             elif source == 'custom_db':
#                 if accession in custom_dict:
#                     entry['description'] = custom_dict[accession]
#                 else:
#                     entry['description'] = accession # mark as accession if no annotation given for custom for now
#             else:
#                 # Keep "hypothetical protein" for missing
#                 entry['description'] = "hypothetical protein"

#         # Write back normalized list or single entry
#         feat['pstc'] = pstc_entries if isinstance(pstc, list) else pstc_entries[0]

#     return features




def fetch_sql_description(conn, source, accession):
    table_map = {
        'swissprot': 'swissprot',
        'afdb': 'afdbclusters',
        'pdb': 'pdb',
        'cath': 'cath',
    }

    table = table_map.get(source)
    if table is None:
        return None

    # special case for cath, which can have multiple top hits (greedy) - multidomain proteins
    if table == 'cath':
        cursor = conn.execute("SELECT product FROM cath WHERE id = ?", (accession,))
    else:
        cursor = conn.execute(f"SELECT product FROM {table} WHERE id = ?", (accession,))
    
    row = cursor.fetchone()
    return row[0] if row else None


def fetch_sql_description_threadsafe(db_path, source, accession):
    """
    makes new connection every time so don't have 2 CATH accessions colliding (for multi domain proteins)
    """
    import sqlite3
    conn = sqlite3.connect(db_path, uri=True, check_same_thread=False)
    try:
        result = fetch_sql_description(conn, source, accession)
    finally:
        conn.close()
    return result

# add custom later
def lookup_sql(features: Sequence[dict], baktfold_db: Path, threads: int):
    """Lookup PSTC information"""
    
    no_pstc_lookups = 0
    # try:
    rec_futures = []
    logger.info("Looking up PSTC descriptions")
    # with sqlite3.connect(f"file:{baktfold_db.joinpath('baktfold.db')}?mode=ro&nolock=1&cache=shared", uri=True, check_same_thread=False) as conn:
    #     conn.execute('PRAGMA omit_readlock;')
    #     conn.row_factory = sqlite3.Row
    with ThreadPoolExecutor(max_workers=max(10, threads)) as tpe:  # use min 10 threads for IO bound non-CPU lookups
        for feat in features:
            pstc = feat.get('pstc')
            if not pstc:
                continue

            # Normalize to list for consistent handling
            pstc_entries = pstc if isinstance(pstc, list) else [pstc]
        
            rec_futures = []
            for entry in pstc_entries:
                accession = entry.get('id')

                source = entry.get('source')

                # submit database query as a future
                future = tpe.submit(fetch_sql_description_threadsafe, baktfold_db.joinpath('baktfold.db'), source, accession)
                rec_futures.append((entry, future))


            # Collect results
            for entry, future in rec_futures:
                desc = future.result()
                if desc:
                    entry['description'] = desc
                else:
                    if entry.get('source') == 'custom_db':
                        entry['description'] = accession  # keep accession if custom_db but missing
                    else:
                        entry['description'] = "hypothetical protein"

        # Write back normalized list or single entry
        feat['pstc'] = pstc_entries if isinstance(pstc, list) else pstc_entries[0]
    
    # except Exception as ex:
    #     logger.error('Could not read PSTCs from db!')
    #     raise Exception('SQL error!', ex)
    # log.info('looked-up=%i', no_pstc_lookups)

    return features

def fetch_db_pscc_result(conn: sqlite3.Connection, uniref50_id: str):
    c = conn.cursor()
    c.execute('select * from pscc where uniref50_id=?', (uniref50_id,))
    rec = c.fetchone()
    c.close()
    return rec


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