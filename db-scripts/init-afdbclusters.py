
import argparse
import hashlib
import logging
import sqlite3

from pathlib import Path

from alive_progress import alive_bar
from Bio import SeqIO
from xml.etree import ElementTree as et
from xopen import xopen

# do i want to filter to bacterial?

parser = argparse.ArgumentParser(
    description='Filter Uniprot\'s filtered trembl XML files to AFDB cluster sequences and create initial PSTC db.'
)
parser.add_argument('--taxonomy', action='store', help='Path to NCBI taxonomy node.dmp file.')
parser.add_argument('--xml', action='store', help='Path to filtered trembl xml file.')
parser.add_argument('--db', action='store', help='Path to Baktfold sqlite3 db file.')
args = parser.parse_args()

DISCARDED_PRODUCTS = [
    'hypothetical protein',
    'hypothetical conserved protein',
    'uncharacterized protein',
    'hypothetical membrane protein',
    'hypothetical cytosolic Protein'
]

taxonomy_path = Path(args.taxonomy).resolve()
xml_path = Path(args.xml).resolve()
db_path = Path(args.db)


logging.basicConfig(
    filename='baktfoldd.db.log',
    filemode='a',
    format='%(name)s - UniRef100 - %(levelname)s - %(message)s',
    level=logging.DEBUG
)
log_pstc = logging.getLogger('PSTC')

def is_taxon_child(child, LCA, taxonomy):
    parent = taxonomy.get(child, None)
    while(parent is not None and parent != '1'):
        if(parent == LCA):
            return True
        else:
            parent = taxonomy.get(parent, None)
    return False


print('parse & store NCBI taxonomy...')
taxonomy = {}
with taxonomy_path.open() as fh:
    for line in fh:
        cols = line.split('\t|\t', maxsplit=2)
        taxonomy[cols[0]] = cols[1]
print(f'\tstored tax ids: {len(taxonomy)}')


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

    with xopen(str(xml_path), mode='rb') as fh_xml, alive_bar(enrich_print=False) as bar:
        for event, elem in et.iterparse(fh_xml, events=('end',)):
            # namespace wildcard
            if elem.tag.endswith('entry'):

                ns0 = {"ns0": "https://uniprot.org/uniprot"}

                # Extract accession
                acc_elem = elem.find('./{*}accession')
                if acc_elem is None:
                    elem.clear()
                    continue
                accession = acc_elem.text.strip()

                if not accession in seen_accessions:

                    # Extract organism info
                    org_elem = elem.find('./{*}organism/{*}name[@type="scientific"]')
                    organism = org_elem.text if org_elem is not None else ""
                    tax_elem = elem.find('./{*}organism/{*}dbReference[@type="NCBI Taxonomy"]')
                    tax_id = tax_elem.get('id') if tax_elem is not None else "1"

                    # Extract protein name



                    # product_elem = elem.find('.//u:recommendedName/u:fullName', ns)

                    #et.dump(elem)

                    product_elem = None

                    recommended = elem.find(".//ns0:protein/ns0:recommendedName/ns0:fullName", ns0)
                    alternative = elem.find(".//ns0:protein/ns0:alternativeName/ns0:fullName", ns0)
                    submitted = elem.find(".//ns0:protein/ns0:submittedName/ns0:fullName", ns0)

                    if recommended is not None:
                        product_elem = recommended
                    elif alternative is not None:
                        product_elem = alternative
                    elif submitted is not None:
                        product_elem = submitted

                    product = product_elem.text if product_elem is not None else None
                    if product and product.lower() in DISCARDED_PRODUCTS:
                        product = "hypothetical protein" # None fails?

                    # Extract sequence
                    seq_elem = elem.find('./{*}sequence')
                    if seq_elem is None:
                        elem.clear()
                        continue
                    seq = seq_elem.text.replace('\n', '').strip().upper()
                    length = len(seq)

                    # Filter for bacteria or phage
                    if is_taxon_child(tax_id, '2', taxonomy) or 'phage' in organism.lower():
                        seq_hash = hashlib.md5(seq.encode())
                        seq_hash_hex = seq_hash.hexdigest()
                        seen_accessions.add(accession)
                        seq_hashes.add(seq_hash_hex)

                    # Store in SQLite

                    conn.execute(
                        "INSERT INTO afdbclusters (id, product) VALUES (?, ?)",
                        (accession, product)
                    )
                    db_updates += 1
                    pstc_seqs += 1

                    elem.clear()
                    bar()

                    if db_updates % 100000 == 0:
                        conn.commit()

                else:
                    print(f"{accession} duplicated")

        conn.commit()

    print(f"✅ Stored {pstc_seqs:,} AFDB Cluster entries")

log_pstc.debug('summary: # AFDBCluster=%i', pstc_seqs)

