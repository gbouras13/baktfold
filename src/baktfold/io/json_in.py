# import atexit
import json
# import logging
# import os
# import sys
# from datetime import datetime
from pathlib import Path

from loguru import logger
from xopen import xopen

import baktfold.bakta.constants as bc
import baktfold.bakta.config as cfg
# import baktfold.utils as bu
# import baktfold.io.fasta as fasta
# import baktfold.io.tsv as tsv
# import baktfold.io.gff as gff
# import baktfold.io.insdc as insdc
# import baktfold.plot as plot


def parse_json_input(input_path, faa_path, all_proteins, protein_json_flag):
    """
    Parses genome annotations from input JSON file.

    Args:
      input_path (str): Path to input JSON file.
      faa_path (str): Path to output file for hypothetical proteins.
      all_proteins (bool): Whether to keep all proteins or only hypothetical ones.
      protein_json_flag (bool): Whether input is protein JSON

    Returns:
      tuple: A tuple containing the data, features, and whether there are duplicate locus tags.

    Examples:
      >>> parse_json_input('input.json', 'hypotheticals.faa', False, False)
      (data, features, False, False)
    """

    

    ############################################################################
    # Checks and configurations
    # - check parameters and setup global configuration
    # - test database
    # - test binary dependencies
    ############################################################################

    try:
        if input_path == '':
            raise ValueError('File path argument must be non-empty')
        annotation_path = Path(input_path).resolve()
        cfg.check_readability('annotation', annotation_path)
        cfg.check_content_size('annotation', annotation_path)
    except:
        logger.error(f'ERROR: annotation file {annotation_path} not valid!')
    
    #print(f'baktfold v{cfg.version}')

    logger.info(f'Parsing annotations from input: {annotation_path}')
    with xopen(str(annotation_path), threads=0) as fh:
        data = json.load(fh)


    features = data['features']

    # features_by_sequence = {seq['id']: [] for seq in data['sequences']}
    # for feature in data['features']:
    #     seq_id = feature['sequence'] if 'sequence' in feature else feature['contig']  # <1.10.0 compatibility
    #     sequence_features = features_by_sequence.get(seq_id)
    #     sequence_features.append(feature)

    # keep all proteins
    if all_proteins:
        hypotheticals = [feat for feat in features if feat['type'] == bc.FEATURE_CDS ]
    else:
        hypotheticals = [feat for feat in features if feat['type'] == bc.FEATURE_CDS and 'hypothetical' in feat]


    if protein_json_flag: # this will also be only hypotheticals if protein mode (or else why not just run with the FASTA)
        return features, hypotheticals


    # check if dupe locus tags (euks can have multiple CDS same locus tag e.g. Cladocopium goreaui CAMXCT020000001.1)
    seen_loci = set()
    has_duplicate_locus = False

    for feat in hypotheticals:
        locus = feat['locus']
        if locus in seen_loci:
            has_duplicate_locus = True
            logger.warning("Multiple CDS per locus tag were detected in your input JSON.")
            logger.warning("CDS id (which is unique) rather than locus tag will be used for ProstT5+Foldseek searches.")
            break
        seen_loci.add(locus)

    

    # this is done after getting all the sequences into the dict for baktfold proteins

    if has_duplicate_locus:
        # write hypothetical proteins to file with id (not locus) as guaranteed exists and unique
        with faa_path.open('wt') as fh:
            for feat in hypotheticals:
                fh.write(f">{feat['id']}\n{feat['aa']}\n")

    else:
        # write hypothetical proteins to file - almost always
        with faa_path.open('wt') as fh:
            for feat in hypotheticals:
                fh.write(f">{feat['locus']}\n{feat['aa']}\n")

    # none of this is relevant for proteins
    try:
        genome_block = data.get("genome")

        if genome_block is None:
            logger.error("No 'genome' block found in input JSON. Please check.")
            translation_table = None
        else:
            if "translation_table" not in genome_block:
                logger.error("No translation table found in input JSON. Please check your input.")
            else:
                raw_value = genome_block["translation_table"]

                try:
                    translation_table = int(raw_value)
                    logger.info(
                        f"Translation table {translation_table} detected from input JSON"
                    )

                except (ValueError, TypeError):
                    translation_table = str(raw_value)
                    logger.warning(
                        f"Translation table '{raw_value}' is not an integer. "
                        f"Parsing it as a string."
                    )

    except Exception as e:
        logger.exception(
            f"Unexpected error while parsing translation table: {e}"
        )
        translation_table = None

    # input detection

    version = data.get("version", {})

    prokka = False
    other_genbank = False

    if "prokka" in version:
        prokka = True
        logger.info("Prokka input detected")
    if  "prokka"  not in version and "bakta" not in version:
        other_genbank = True

    logger.info('Parsing complete')

    return data, features, has_duplicate_locus, translation_table, prokka, other_genbank


def log_for_other_genbank_tools(cds_program,trna_program, rrna_program, tmrna_program, ncrna_program):

    logger.warning("Neither bakta nor prokka input detected")
    logger.info("If you would like to specify consituent inference tools for CDS, tRNA, rRNA, tmRNA and ncRNA")
    logger.info("Reminder: please use --cds-tool --trna-program --rrna-program  --tmrna-program --ncrna-program to modify them if you haven't already")
    logger.info(f"For this genome, they are --cds_program {cds_program} --trna-program {trna_program} --rrna_program {rrna_program} --tmrna-program {tmrna_program} --ncrna-program  {ncrna_program}")


    