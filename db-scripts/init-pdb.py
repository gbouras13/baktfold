
import argparse
import hashlib
import logging
import sqlite3
import csv

from pathlib import Path

from alive_progress import alive_bar
from Bio import SeqIO
from xml.etree import ElementTree as et
from xopen import xopen

# do i want to filter to bacterial?

parser = argparse.ArgumentParser(
    description='Add PDB entries and create initial PSTC db.'
)
parser.add_argument('--tsv', action='store', help='Path to pdb.tsv')
parser.add_argument('--db', action='store', help='Path to Baktfold sqlite3 db file.')
args = parser.parse_args()

DISCARDED_PRODUCTS = [
    'hypothetical protein',
    'hypothetical conserved protein',
    'uncharacterized protein',
    'hypothetical membrane protein',
    'hypothetical cytosolic Protein'
]

pdb_tsv_path = Path(args.tsv).resolve()
db_path = Path(args.db)

logging.basicConfig(
    filename='baktfold.db.log',
    filemode='a',
    format='%(name)s - UniRef100 - %(levelname)s - %(message)s',
    level=logging.DEBUG
)
log_pstc = logging.getLogger('PSTC')


print('parse & store PSTC information...')
seq_hashes = set()  # unique protein sequences
seen_accessions = set()

with sqlite3.connect(str(db_path), isolation_level='EXCLUSIVE') as conn:
    conn.execute('PRAGMA page_size = 4096;')
    conn.execute('PRAGMA cache_size = 100000;')
    conn.execute('PRAGMA locking_mode = EXCLUSIVE;')
    conn.execute(f'PRAGMA mmap_size = {20 * 1024 * 1024 * 1024};')
    conn.execute('PRAGMA synchronous = OFF;')
    conn.execute('PRAGMA journal_mode = OFF')
    conn.execute('PRAGMA threads = 2;')
    conn.commit()
    conn.row_factory = sqlite3.Row

    db_updates = 0
    ups_seqs = 0
    pstc_seqs = 0

    with pdb_tsv_path.open(mode='r', encoding='utf-8') as fh_tsv, alive_bar(sum(1 for _ in fh_tsv), enrich_print=False) as bar:
        fh_tsv.seek(0)  # reset file pointer after counting lines
        reader = csv.reader(fh_tsv, delimiter='\t')
        for row in reader:
            if len(row) < 2:
                continue  # skip malformed lines
            accession, product = row[0], row[1]

            try:
                conn.execute(
                    "INSERT INTO pdb (id, product) VALUES (?, ?)",
                    (accession, product)
                )
                db_updates += 1
                pstc_seqs += 1
            except sqlite3.IntegrityError:
                print(f"{accession} duplicated")

            bar()

            if db_updates % 100000 == 0:
                conn.commit()

    # Final commit
    conn.commit()

    print(f"✅ Stored {pstc_seqs:,} PDB entries")

log_pstc.debug('summary: # PDB=%i', pstc_seqs)

