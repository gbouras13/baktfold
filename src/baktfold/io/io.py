from loguru import logger
from pathlib import Path
import baktfold.bakta.config as cfg
import baktfold.io.gff as gff
import baktfold.io.tsv as tsv
import baktfold.io.insdc as insdc
import baktfold.io.fasta as fasta
import baktfold.io.json as json
import baktfold.bakta.constants as bc
from typing import Sequence
import pandas as pd

"""
wrapper script over all io output submodules
"""

def write_foldseek_tophit(tophit_df: pd.DataFrame, pdb_tophit_path: Path):
    logger.info(f"Saving foldseek tophits against PDB to {pdb_tophit_path}")
    tophit_df.to_csv(pdb_tophit_path, sep="\t", index=False)


def write_bakta_outputs(data: dict, features: Sequence[dict], features_by_sequence: Sequence[dict] , output: Path, prefix: str):

    #logger.info(f'selected features={len(features)}')

    logger.info('writing human readable TSV...')
    tsv_path: Path = Path(output) / f"{prefix}.tsv"
    tsv.write_features(data['sequences'], features_by_sequence, tsv_path)

    logger.info('writing GFF3...')
    gff3_path: Path = Path(output) / f"{prefix}.gff3"
    gff.write_features(data, features_by_sequence, gff3_path)

    logger.info('writing INSDC GenBank & EMBL...')
    genbank_path: Path = Path(output) / f"{prefix}.gbff"
    embl_path: Path = Path(output) / f"{prefix}.embl"
    insdc.write_features(data, features, genbank_path, embl_path)

    # add 3Di sequence here I think or debate how to handle this

    logger.info('writing genome sequences...')
    fna_path: Path = Path(output) / f"{prefix}.fna"
    fasta.export_sequences(data['sequences'], fna_path, description=True, wrap=True)

    logger.info('writing feature nucleotide sequences...')
    ffn_path: Path = Path(output) / f"{prefix}.ffn"
    fasta.write_ffn(features, ffn_path)

    logger.info('writing translated CDS sequences...')
    faa_path: Path = Path(output) / f"{prefix}.faa"
    fasta.write_faa(features, faa_path)

    # inference here is the different databases?
    annotations_path: Path = Path(output) / f"{prefix}.inference.tsv"
    header_columns = ['ID', 'Length', 'Product', 'Swissprot', 'AFDBClusters', 'PDB']
    logger.info(f'Exporting annotations (TSV) to: {annotations_path}')

    selected_features = []

    for seq_id, features in features_by_sequence.items():
        for feat in features:
            # get() ensures we don't crash if the key doesn't exist
            if 'hypothetical' in feat or 'baktfold' in feat:
                selected_features.append(feat)


    tsv.write_protein_features(selected_features, header_columns, annotations_path)
    
    

    cfg.skip_cds = False
    if(cfg.skip_cds is False):

        # no need to write the hypotheticals I think

        # hypotheticals = [feat for feat in features if feat['type'] == bc.FEATURE_CDS and 'hypothetical' in feat]


        # print('writing hypothetical TSV...')
        # tsv_path: Path = Path(output) / f"{prefix}.hypotheticals.tsv"
        # tsv.write_hypotheticals(hypotheticals, tsv_path)

        # print('writing translated hypothetical CDS sequences...')
        # print('writing translated CDS sequences...')
        # faa_path: Path = Path(output) / f"{prefix}.hypotheticals.faa"
        # fasta.write_faa(hypotheticals, faa_path)

        # calc & store runtime

        # run_duration = (cfg.run_end - cfg.run_start).total_seconds()
        # data['run'] = {
        #     'start': cfg.run_start.strftime('%Y-%m-%d %H:%M:%S'),
        #     'end': cfg.run_end.strftime('%Y-%m-%d %H:%M:%S'),
        #     'duration': f'{(run_duration / 60):.2f} min'
        # }

        logger.info('write machine readable JSON...')
        json_path: Path = Path(output) / f"{prefix}.json"
        json.write_json(data, features, json_path)
