import atexit
import json
import logging
import os
import sys
from loguru import logger

from datetime import datetime
from pathlib import Path

from xopen import xopen

import baktfold.bakta.constants as bc
import baktfold.bakta.config as cfg

# import baktfold.utils as bu
# import baktfold.io.fasta as fasta
# import baktfold.io.tsv as tsv
# import baktfold.io.gff as gff
# import baktfold.io.insdc as insdc
# import baktfold.plot as plot

def parse_json_input(input_path, faa_path, all_proteins):
    """
    Parses genome annotations from input JSON file.

    Args:
      input_path (str): Path to input JSON file.
      faa_path (str): Path to output file for hypothetical proteins.
      all_proteins (bool): Whether to keep all proteins or only hypothetical ones.

    Returns:
      tuple: A tuple containing the data, features, and whether there are duplicate locus tags.

    Examples:
      >>> parse_json_input('input.json', 'hypotheticals.faa', False)
      (data, features, False)
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

    logger.info(f'Parsing genome annotations from input: {annotation_path}')
    with xopen(str(annotation_path), threads=0) as fh:
        data = json.load(fh)
    features = data['features']
    features_by_sequence = {seq['id']: [] for seq in data['sequences']}
    for feature in data['features']:
        seq_id = feature['sequence'] if 'sequence' in feature else feature['contig']  # <1.10.0 compatibility
        sequence_features = features_by_sequence.get(seq_id)
        sequence_features.append(feature)

    # keep all proteins
    if all_proteins:
        hypotheticals = [feat for feat in features if feat['type'] == bc.FEATURE_CDS ]
    else:

        hypotheticals = [feat for feat in features if feat['type'] == bc.FEATURE_CDS and 'hypothetical' in feat]


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


    logger.info('Parsing complete')

    return data, features, has_duplicate_locus


    
#     # set global config objects based on information from imported JSON document
#     cfg.version = data['version']['baktfold']
#     cfg.db_info = {
#         'type': data['version']['db']['type'],
#         'major': data['version']['db']['version'].split('.')[0],
#         'minor': data['version']['db']['version'].split('.')[1]
#     }
#     cfg.translation_table = data['genome']['translation_table']
#     cfg.run_start = datetime.strptime(data['run']['start'], '%Y-%m-%d %H:%M:%S')
#     cfg.run_end = datetime.strptime(data['run']['end'], '%Y-%m-%d %H:%M:%S')

#     ############################################################################
#     # Write output files
#     # - write optional output files in GFF3/GenBank/EMBL formats
#     # - measure runtime
#     # - write comprehensive annotation results as JSON
#     # - remove temp directory
#     ############################################################################
#     print(f'\nExport annotation results to: {cfg.output_path}')
#     print('\thuman readable TSV...')
#     tsv_path = cfg.output_path.joinpath(f'{cfg.prefix}.tsv')
#     tsv.write_features(data['sequences'], features_by_sequence, tsv_path)

#     print('\tGFF3...')
#     gff3_path = cfg.output_path.joinpath(f'{cfg.prefix}.gff3')
#     gff.write_features(data, features_by_sequence, gff3_path)

#     print('\tINSDC GenBank & EMBL...')
#     genbank_path = cfg.output_path.joinpath(f'{cfg.prefix}.gbff')
#     embl_path = cfg.output_path.joinpath(f'{cfg.prefix}.embl')
#     insdc.write_features(data, features, genbank_path, embl_path)

#     print('\tgenome sequences...')
#     fna_path = cfg.output_path.joinpath(f'{cfg.prefix}.fna')
#     fasta.export_sequences(data['sequences'], fna_path, description=True, wrap=True)

#     print('\tfeature nucleotide sequences...')
#     ffn_path = cfg.output_path.joinpath(f'{cfg.prefix}.ffn')
#     fasta.write_ffn(features, ffn_path)

#     print('\ttranslated CDS sequences...')
#     faa_path = cfg.output_path.joinpath(f'{cfg.prefix}.faa')
#     fasta.write_faa(features, faa_path)

#     print('\tfeature inferences...')
#     tsv_path = cfg.output_path.joinpath(f'{cfg.prefix}.inference.tsv')
#     tsv.write_feature_inferences(data['sequences'], features_by_sequence, tsv_path)

#     print('\tcircular genome plot...')
#     plot.write(data, features, cfg.output_path)

#     hypotheticals = [feat for feat in features if feat['type'] == bc.FEATURE_CDS and 'hypothetical' in feat]
#     print('\thypothetical TSV...')
#     tsv_path = cfg.output_path.joinpath(f'{cfg.prefix}.hypotheticals.tsv')
#     tsv.write_hypotheticals(hypotheticals, tsv_path)
    
#     print('\ttranslated hypothetical CDS sequences...')
#     faa_path = cfg.output_path.joinpath(f'{cfg.prefix}.hypotheticals.faa')
#     fasta.write_faa(hypotheticals, faa_path)

#     print('\tGenome and annotation summary...')
#     summary_path = cfg.output_path.joinpath(f'{cfg.prefix}.txt')
#     with summary_path.open('w') as fh_out:
#         fh_out.write('Sequence(s):\n')
#         fh_out.write(f"Length: {data['stats']['size']:}\n")
#         fh_out.write(f"Count: {len(data['sequences'])}\n")
#         fh_out.write(f"GC: {100 * data['stats']['gc']:.1f}\n")
#         fh_out.write(f"N50: {data['stats']['n50']:}\n")
#         if('n90' in data['stats']):
#             fh_out.write(f"N90: {data['stats']['n90']:}\n")
#         fh_out.write(f"N ratio: {100 * data['stats']['n_ratio']:.1f}\n")
#         fh_out.write(f"coding density: {100 * data['stats']['coding_ratio']:.1f}\n")
#         fh_out.write('\nAnnotation:\n')
#         fh_out.write(f"tRNAs: {len([feat for feat in features if feat['type'] == bc.FEATURE_T_RNA])}\n")
#         fh_out.write(f"tmRNAs: {len([feat for feat in features if feat['type'] == bc.FEATURE_TM_RNA])}\n")
#         fh_out.write(f"rRNAs: {len([feat for feat in features if feat['type'] == bc.FEATURE_R_RNA])}\n")
#         fh_out.write(f"ncRNAs: {len([feat for feat in features if feat['type'] == bc.FEATURE_NC_RNA])}\n")
#         fh_out.write(f"ncRNA regions: {len([feat for feat in features if feat['type'] == bc.FEATURE_NC_RNA_REGION])}\n")
#         fh_out.write(f"CRISPR arrays: {len([feat for feat in features if feat['type'] == bc.FEATURE_CRISPR])}\n")
#         cdss = [feat for feat in features if feat['type'] == bc.FEATURE_CDS]
#         fh_out.write(f"CDSs: {len(cdss)}\n")
#         fh_out.write(f"pseudogenes: {len([cds for cds in cdss if 'pseudogene' in cds])}\n")
#         fh_out.write(f"hypotheticals: {len([cds for cds in cdss if 'hypothetical' in cds])}\n")
#         fh_out.write(f"sORFs: {len([feat for feat in features if feat['type'] == bc.FEATURE_SORF])}\n")
#         fh_out.write(f"gaps: {len([feat for feat in features if feat['type'] == bc.FEATURE_GAP])}\n")
#         fh_out.write(f"oriCs: {len([feat for feat in features if feat['type'] == bc.FEATURE_ORIC])}\n")
#         fh_out.write(f"oriVs: {len([feat for feat in features if feat['type'] == bc.FEATURE_ORIV])}\n")
#         fh_out.write(f"oriTs: {len([feat for feat in features if feat['type'] == bc.FEATURE_ORIT])}\n")
#         fh_out.write('\nbaktfold:\n')
#         fh_out.write(f'Software: v{cfg.version}\n')
#         fh_out.write(f"Database: v{cfg.db_info['major']}.{cfg.db_info['minor']}, {cfg.db_info['type']}\n")
#         fh_out.write('DOI: 10.1099/mgen.0.000685\n')
#         fh_out.write('URL: github.com/oschwengers/baktfold\n')


# if __name__ == '__main__':
#     main()
