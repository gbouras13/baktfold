from loguru import logger

from pathlib import Path
from types import LambdaType
from typing import Dict, Sequence
from datetime import datetime
import baktfold.bakta
import baktfold.bakta.config as cfg
import baktfold.bakta.constants as bc
# import bakta.ips as bips
# import bakta.ups as bups
# import bakta.psc as bpsc
# import bakta.pscc as bpscc




# log = logging.getLogger('TSV')


def write_features(sequences: Sequence[dict], features_by_sequence: Dict[str, dict], tsv_path: Path):
    """Export features in TSV format."""
    logger.info(f'write feature tsv: path={tsv_path}')

    with tsv_path.open('wt') as fh:
        fh.write('# Annotated with Bakta\n')
        fh.write(f'# Software: v{cfg.version}\n')
        fh.write(f"# Database: v{cfg.version}\n") # fix later
        #fh.write(f"# Database: v{cfg.db_info['major']}.{cfg.db_info['minor']}, {cfg.db_info['type']}\n")
        fh.write(f'# DOI: {bc.BAKTFOLD_DOI}\n')
        fh.write(f'# URL: {bc.BAKTFOLD_URL}\n')
        fh.write('#Sequence Id\tType\tStart\tStop\tStrand\tLocus Tag\tGene\tProduct\tDbXrefs\n')

        for seq in sequences:
            for feat in features_by_sequence[seq['id']]:
                seq_id = feat['sequence'] if 'sequence' in feat else feat['contig']  # <1.10.0 compatibility
                feat_type = feat['type']
                if(feat_type == bc.FEATURE_GAP):
                    feat_type = bc.INSDC_FEATURE_ASSEMBLY_GAP if feat['length'] >= 100 else bc.INSDC_FEATURE_GAP

                gene = feat['gene'] if feat.get('gene', None) else ''
                product = feat.get('product', '')
                if(bc.PSEUDOGENE in feat):
                    product = f"(pseudo) {product}"
                elif(feat.get('truncated', '') == bc.FEATURE_END_5_PRIME):
                    product = f"(5' truncated) {product}"
                elif(feat.get('truncated', '') == bc.FEATURE_END_3_PRIME):
                    product = f"(3' truncated) {product}"
                elif(feat.get('truncated', '') == bc.FEATURE_END_BOTH):
                    product = f"(partial) {product}"
                fh.write('\t'.join(
                    [
                        seq_id,
                        feat_type,
                        str(feat['start']),
                        str(feat['stop']),
                        feat['strand'],
                        feat.get('locus', ''),
                        gene,
                        product,
                        ', '.join(sorted(feat.get('db_xrefs', [])))
                    ])
                )
                fh.write('\n')
                if(feat_type == bc.FEATURE_CRISPR):
                    i = 0
                    while i < len(feat['spacers']):
                        repeat = feat['repeats'][i]
                        fh.write('\t'.join([seq_id, bc.FEATURE_CRISPR_REPEAT, str(repeat['start']), str(repeat['stop']), repeat['strand'], '', '', f"CRISPR repeat", '']))
                        fh.write('\n')
                        spacer = feat['spacers'][i]
                        fh.write('\t'.join([seq_id, bc.FEATURE_CRISPR_SPACER, str(spacer['start']), str(spacer['stop']), spacer['strand'], '', '', f"CRISPR spacer, sequence {spacer['sequence']}", '']))
                        fh.write('\n')
                        i += 1
                    if(len(feat['repeats']) - 1 == i):
                        repeat = feat['repeats'][i]
                        fh.write('\t'.join([seq_id, bc.FEATURE_CRISPR_REPEAT, str(repeat['start']), str(repeat['stop']), repeat['strand'], '', '', f"CRISPR repeat", '']))
                        fh.write('\n')
    return

DB_IPS_COL_UNIREF100 = 'uniref100_id'
DB_UPS_COL_UNIPARC = 'uniparc_id'
DB_PSC_COL_UNIREF90 = 'uniref90_id'
DB_PSCC_COL_UNIREF50 = 'uniref50_id'

def write_feature_inferences(sequences: Sequence[dict], features_by_sequence: Dict[str, dict], tsv_path: Path):
    """Export feature inference statistics in TSV format."""
    logger.info('write tsv: path=%s', tsv_path)

    with tsv_path.open('wt') as fh:
        fh.write('# Annotated with Bakta\n')
        fh.write(f'# Software: v{cfg.version}\n')
        fh.write(f"# Database: v{cfg.version}\n") # fix later
        #fh.write(f"# Database: v{cfg.db_info['major']}.{cfg.db_info['minor']}, {cfg.db_info['type']}\n")
        fh.write(f'# DOI: {bc.BAKTFOLD_DOI}\n')
        fh.write(f'# URL: {bc.BAKTFOLD_URL}\n')
        fh.write('#Sequence Id\tType\tStart\tStop\tStrand\tLocus Tag\tScore\tEvalue\tQuery Cov\tSubject Cov\tId\tAccession\n')

        for seq in sequences:
            for feat in features_by_sequence[seq['id']]:
                if(feat['type'] in [bc.FEATURE_CDS, bc.FEATURE_SORF]):
                    score, evalue, query_cov, subject_cov, identity, accession = None, None, None, None, None, '-'
                    if('ups' in feat or 'ips' in feat):
                        query_cov = 1
                        subject_cov = 1
                        identity = 1
                        evalue = 0
                        accession = f"{bc.DB_XREF_UNIREF}:{feat['ips'][DB_IPS_COL_UNIREF100]}" if 'ips' in feat else f"{bc.DB_XREF_UNIPARC}:{feat['ups'][DB_UPS_COL_UNIPARC]}"
                    elif('psc' in feat or 'pscc' in feat):
                        psc_type = 'psc' if 'psc' in feat else 'pscc'
                        query_cov = feat[psc_type]['query_cov']
                        subject_cov = feat[psc_type].get('subject_cov', -1)
                        identity = feat[psc_type]['identity']
                        score = feat[psc_type].get('score', -1)
                        evalue = feat[psc_type].get('evalue', -1)
                        accession = f"{bc.DB_XREF_UNIREF}:{feat['psc'][DB_PSC_COL_UNIREF90]}" if 'psc' in feat else f"{bc.DB_XREF_UNIREF}:{feat['pscc'][DB_PSCC_COL_UNIREF50]}"
                    fh.write('\t'.join(
                        [
                            feat['sequence'] if 'sequence' in feat else feat['contig'],  # <1.10.0 compatibility
                            feat['type'],
                            str(feat['start']),
                            str(feat['stop']),
                            feat['strand'],
                            feat['locus'],
                            f"{score:0.1f}" if score != None else '-',
                            ('0.0' if evalue == 0 else f"{evalue:1.1e}") if evalue != None else '-',
                            ('1.0' if query_cov == 1 else f"{query_cov:0.3f}") if query_cov != None else '-',
                            ('1.0' if subject_cov == 1 else f"{subject_cov:0.3f}") if subject_cov != None else '-',
                            ('1.0' if identity == 1 else f"{identity:0.3f}") if identity != None else '-',
                            accession
                        ])
                    )
                    fh.write('\n')
                elif(feat['type'] in [bc.FEATURE_T_RNA, bc.FEATURE_R_RNA, bc.FEATURE_NC_RNA, bc.FEATURE_NC_RNA_REGION]):
                    accession = '-' if feat['type'] == bc.FEATURE_T_RNA else [xref for xref in feat['db_xrefs'] if bc.DB_XREF_RFAM in xref][0]
                    fh.write('\t'.join(
                        [
                            feat['sequence'] if 'sequence' in feat else feat['contig'],  # <1.10.0 compatibility
                            feat['type'],
                            str(feat['start']),
                            str(feat['stop']),
                            feat['strand'],
                            feat['locus'] if 'locus' in feat else '-',
                            f"{feat['score']:0.1f}",
                            ('0.0' if feat['evalue'] == 0 else f"{feat['evalue']:1.1e}") if 'evalue' in feat else '-',
                            ('1.0' if feat['query_cov'] == 1 else f"{feat['query_cov']:0.3f}") if 'query_cov' in feat else '-',
                            ('1.0' if feat['subject_cov'] == 1 else f"{feat['subject_cov']:0.3f}") if 'subject_cov' in feat else '-',
                            ('1.0' if feat['identity'] == 1 else f"{feat['identity']:0.3f}") if 'identity' in feat else '-',
                            accession
                        ])
                    )
                    fh.write('\n')
    return

def map_aa_columns(feat: dict, custom_db: bool) -> Sequence[str]:
    # no gene here for now
    # gene = feat.get('gene', None)
    # if(gene is None):
    #     gene = ''

    # header_columns = ['Locus', 'Length', 'Product', 'Swissprot', 'AFDBClusters', 'PDB']

    # for the GenBank
    if 'length' not in feat:
        feat['length'] = int(len(feat['nt'])/3)

    if custom_db:

        return [
            feat['locus'],
            str(feat['length']),
            #gene,
            feat['product'],
            ','.join([dbxref.replace('afdb_v6:', '') for dbxref in feat['db_xrefs'] if 'swissprot' in dbxref]),
            ','.join([dbxref.replace('afdb_v6:', '') for dbxref in feat['db_xrefs'] if 'afdbclusters_' in dbxref]),
            ','.join([dbxref.replace('pdb:', '') for dbxref in feat['db_xrefs'] if 'pdb:' in dbxref]),
            ','.join([dbxref.replace('cath:', '') for dbxref in feat['db_xrefs'] if 'cath:' in dbxref]),
            ','.join([dbxref.replace('custom:custom_', '') for dbxref in feat['db_xrefs'] if 'custom:' in dbxref]),
        ]
    else:
        return [
            feat['locus'],
            str(feat['length']),
            #gene,
            feat['product'],
            ','.join([dbxref.replace('afdb_v6:', '') for dbxref in feat['db_xrefs'] if 'swissprot' in dbxref]),
            ','.join([dbxref.replace('afdb_v6:', '') for dbxref in feat['db_xrefs'] if 'afdbclusters_' in dbxref]),
            ','.join([dbxref.replace('pdb:', '') for dbxref in feat['db_xrefs'] if 'pdb:' in dbxref]),
            ','.join([dbxref.replace('cath:', '') for dbxref in feat['db_xrefs'] if 'cath:' in dbxref]),
        ]

def write_protein_features(features: Sequence[dict], header_columns: Sequence[str], tsv_path: Path, custom_db: bool):
    """Export protein features in TSV format."""
    logger.info(f'write protein feature tsv: path={tsv_path}')

    with tsv_path.open('wt') as fh:
        fh.write(f'#Annotated with Baktfold (v{cfg.version}): https://github.com/gbouras13/baktfold\n')
        #fh.write(f"#Database (v{cfg.db_info['major']}.{cfg.db_info['minor']}): https://doi.org/10.5281/zenodo.4247252\n")
        fh.write('\t'.join(header_columns))
        fh.write('\n')
        for feat in features:
            columns = map_aa_columns(feat, custom_db)
            fh.write('\t'.join(columns))
            fh.write('\n')
    return


def write_hypotheticals(hypotheticals: Sequence[dict], tsv_path: Path):
    """Export hypothetical information in TSV format."""
    logger.info('write hypothetical tsv: path=%s', tsv_path)

    with tsv_path.open('wt') as fh:
        fh.write(f'#Annotated with Bakta v{cfg.version}, https://github.com/oschwengers/bakta\n')
        #fh.write(f"#Database v{cfg.db_info['major']}.{cfg.db_info['minor']}, https://doi.org/10.5281/zenodo.4247252\n")
        fh.write('#Sequence Id\tStart\tStop\tStrand\tLocus Tag\tMol Weight [kDa]\tIso El. Point\tPfam hits\tDbxrefs\n')
        for hypo in hypotheticals:
            pfams = [f"{pfam['id']}|{pfam['name']}" for pfam in hypo.get('pfams', [])]
            seq_stats = hypo['seq_stats']
            mol_weight = f"{(seq_stats['molecular_weight']/1000):.1f}" if seq_stats['molecular_weight'] else 'NA'
            iso_point = f"{seq_stats['isoelectric_point']:.1f}" if seq_stats['isoelectric_point'] else 'NA'
            seq_id = hypo['sequence'] if 'sequence' in hypo else hypo['contig']  # <1.10.0 compatibility
            fh.write(f"{seq_id}\t{hypo['start']}\t{hypo['stop']}\t{hypo['strand']}\t{hypo.get('locus', '')}\t{mol_weight}\t{iso_point}\t{', '.join(sorted(pfams))}\t{', '.join(sorted(hypo.get('db_xrefs', [])))}\n")
    return
