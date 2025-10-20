
import argparse
import logging
import sqlite3
from pathlib import Path

parser = argparse.ArgumentParser(
    description='Create baktfold sql db.'
)
parser.add_argument('--db', action='store', help='Path to baktfold sqlite3 db file.')
args = parser.parse_args()


db_path = Path(args.db)


logging.basicConfig(
    filename='baktfold.db.log',
    filemode='w',
    format='%(name)s - INIT - %(levelname)s - %(message)s',
    level=logging.INFO
)

log_pstc = logging.getLogger('PSTC')


with sqlite3.connect(str(db_path), isolation_level='EXCLUSIVE') as conn:
    conn.execute('PRAGMA page_size = 4096;')
    conn.execute('PRAGMA cache_size = 100000;')
    conn.execute('PRAGMA locking_mode = EXCLUSIVE;')
    conn.execute(f'PRAGMA mmap_size = {20 * 1024 * 1024 * 1024};')
    conn.execute('PRAGMA synchronous = OFF;')
    conn.execute('PRAGMA journal_mode = OFF')
    conn.execute('PRAGMA threads = 2;')

    print('create SQL table AFDBClusters...')
    conn.execute('DROP TABLE IF EXISTS afdbclusters;')
    log_pstc.info('DROP TABLE IF EXISTS afdbclusters;')
    stmt = '''CREATE TABLE afdbclusters (
        id TEXT PRIMARY KEY,
        product INTEGER NOT NULL
        ) WITHOUT ROWID;'''
    stmt = ' '.join(stmt.replace('\n', '').split())
    conn.execute(stmt)
    log_pstc.info(stmt)
    conn.commit()
    print('\t...done')

    print('create SQL table Swissprot...')
    conn.execute('DROP TABLE IF EXISTS swissprot;')
    log_pstc.info('DROP TABLE IF EXISTS swissprot;')
    stmt = '''CREATE TABLE swissprot (
        id TEXT PRIMARY KEY,
        product INTEGER NOT NULL
        ) WITHOUT ROWID;'''
    stmt = ' '.join(stmt.replace('\n', '').split())
    conn.execute(stmt)
    log_pstc.info(stmt)
    conn.commit()
    print('\t...done')

    print('create SQL table PDB...')
    conn.execute('DROP TABLE IF EXISTS pdb;')
    log_pstc.info('DROP TABLE IF EXISTS pdb;')
    stmt = '''CREATE TABLE pdb (
        id TEXT PRIMARY KEY,
        product INTEGER NOT NULL
        ) WITHOUT ROWID;'''
    stmt = ' '.join(stmt.replace('\n', '').split())
    conn.execute(stmt)
    log_pstc.info(stmt)
    conn.commit()
    print('\t...done')


print(f'\nSQLite bakta db successfully created: {db_path}')
