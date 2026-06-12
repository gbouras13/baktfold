import sys
from pathlib import Path

from loguru import logger
from xopen import xopen

import baktfold.io.fasta as fasta
from baktfold.io.json_in import parse_json_input, log_for_other_genbank_tools
import baktfold.bakta.constants as bc
from baktfold.utils.validation import is_fasta
# import baktfold.bakta.config as cfg


def parse_protein_input(input_path, faa_path):
    """
    handles regular FASTA and gzipped 
    returns cds_dict
    """

    # handles regular FASTA and gzipped
    fasta_flag = False
    try:
        if input_path == '':
            raise ValueError('File path argument must be non-empty')
        input_path = Path(input_path).resolve()
        fasta_flag = is_fasta(input_path)
        if fasta_flag:
            logger.info('FASTA input format detected.')
        else:
            logger.info('Bakta JSON input format detected. Only hypothetical proteins from the Bakta JSON input file will be annotated.')
    except Exception as e:
        logger.error(f'ERROR: annotation file {input_path} not valid! {e}')
        sys.exit(1)

    aas = []
    hypotheticals = []
    bakta_version = {}
    try:
        if fasta_flag:
            logger.info('Attempting to parse input protein sequences as .faa format ...')
            aas = fasta.import_sequences(input_path, False, False)
        else:
            logger.info('Attempting to parse input protein sequences as Bakta JSON format ...')
            aas, hypotheticals, bakta_version = parse_json_input(input_path, False, False, protein_json_flag=True)
        logger.info(f'Imported sequences={len(aas)}')
    except Exception as e:
        logger.error(f'ERROR: wrong file format or unallowed characters in amino acid sequences! {e}')
        sys.exit(1)
    
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

    
    if fasta_flag:
        with faa_path.open('wt') as fh:
            for aa in aas:
                fh.write(f">{aa['locus']}\n{aa['aa']}\n")
    else: # write hypothetical proteins to file if JSON input
        with faa_path.open('wt') as fh:
            for aa in hypotheticals:
                fh.write(f">{aa['locus']}\n{aa['aa']}\n")

    logger.info('Parsing complete')

    return aas, bakta_version


   