# import atexit
# import json
# import logging
# import os
# import sys
# import binascii
# import gzip
# from typing import List, Tuple, TextIO
# from datetime import datetime
from pathlib import Path

from loguru import logger
from xopen import xopen

import baktfold.io.fasta as fasta
import baktfold.bakta.constants as bc
# import baktfold.bakta.config as cfg


def parse_protein_input(input_path, faa_path):
    """
    handles regular FASTA and gzipped 
    returns cds_dict
    """

    # handles regular FASTA and gzipped 

    try:
        if input_path == '':
            raise ValueError('File path argument must be non-empty')
        input_path = Path(input_path).resolve()
    except:
        logger.error(f'ERROR: annotation file {input_path} not valid!')


    try:
        logger.info('Parsing input protein sequences...')
        aas = fasta.import_sequences(input_path, False, False)
        logger.info(f'Imported sequences={len(aas)}')
    except:
        logger.error('ERROR: wrong file format or unallowed characters in amino acid sequences!')
    
    mock_start = 1
    for aa in aas:  # rename and mock feature attributes to reuse existing functions
        aa['type'] = bc.FEATURE_CDS
        aa['locus'] = aa['id']
        aa['sequence'] = '-'
        aa['start'] = mock_start
        aa['stop'] = mock_start + aa['length'] - 1
        aa['strand'] = bc.STRAND_UNKNOWN
        aa['frame'] = 1
        mock_start += 100

    # write hypothetical proteins to file
    with faa_path.open('wt') as fh:
        for aa in aas:
            fh.write(f">{aa['locus']}\n{aa['aa']}\n")

    logger.info('Parsing complete')

    return aas


   