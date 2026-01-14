import re
# from datetime import date
from pathlib import Path
from typing import Sequence, Tuple

from loguru import logger
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation, AfterPosition, BeforePosition

import baktfold.bakta.config as cfg
import baktfold.bakta.constants as bc
import baktfold.bakta.annotation as ba
import baktfold.bakta.so as so

# log = logging.getLogger('INSDC')

def move_product_to_note_if_exists(qualifiers):
    """
    If a 'product' qualifier exists, append it to 'note' and remove 'product'.

    Designed for the eukaryotic entries

    Parameters
    ----------
    qualifiers : dict
        Feature qualifiers dictionary (values are usually lists).

    Returns
    -------
    None
        Modifies qualifiers in place.
    """
    product = qualifiers.get("product")
    if not product:
        return

    # Ensure note exists and is a list
    if "note" not in qualifiers:
        qualifiers["note"] = []

    if isinstance(product, list):
        qualifiers["note"].extend(product)
    else:
        qualifiers["note"].append(product)

    qualifiers.pop("product", None)

def build_biopython_sequence_list(data: dict, features: Sequence[dict], prokka, euk):

    # need to pass the transl table from bakta or prokka not the config
    # best to just detect it from the input json
    # if prokka:
    #     translation_table

    # need to give euk to not give trnas as trnascan - do

    sequence_list = []
    for seq in data['sequences']:
        sequence_features = []
        if len(features) > 0:
            sequence_features = [feat for feat in features if feat['sequence'] == seq['id']] if 'sequence' in features[0] else [feat for feat in features if feat['contig'] == seq['id']]  # <1.10.0 compatibility
        comment = (
            'Annotated with Baktfold',
            f"Software: v{cfg.version}\n",
            f"Database: v{cfg.version}\n", # fix later
            #f"Database: v{cfg.db_info['major']}.{cfg.db_info['minor']}, {cfg.db_info['type']}\n",
            f'DOI: {bc.BAKTFOLD_DOI}\n',
            f'URL: {bc.BAKTFOLD_URL}\n',
            '\n',
            '##Genome Annotation Summary:##\n',
            #f"{'Annotation Date':<30} :: {cfg.run_end.strftime('%m/%d/%Y, %H:%M:%S')}\n",
            f"{'CDSs':<30} :: {len([feat for feat in sequence_features if feat['type'] == bc.FEATURE_CDS or feat['type'] == bc.FEATURE_SORF]):5,}\n",
            f"{'tRNAs':<30} :: {len([feat for feat in sequence_features if feat['type'] == bc.FEATURE_T_RNA]):5,}\n",
            f"{'tmRNAs':<30} :: {len([feat for feat in sequence_features if feat['type'] == bc.FEATURE_TM_RNA]):5,}\n",
            f"{'rRNAs':<30} :: {len([feat for feat in sequence_features if feat['type'] == bc.FEATURE_R_RNA]):5,}\n",
            f"{'ncRNAs':<30} :: {len([feat for feat in sequence_features if feat['type'] == bc.FEATURE_NC_RNA]):5,}\n",
            f"{'regulatory ncRNAs':<30} :: {len([feat for feat in sequence_features if feat['type'] == bc.FEATURE_NC_RNA_REGION]):5,}\n",
            f"{'CRISPR Arrays':<30} :: {len([feat for feat in sequence_features if feat['type'] == bc.FEATURE_CRISPR]):5,}",
            f"{'oriCs/oriVs':<30} :: {len([feat for feat in sequence_features if feat['type'] == bc.FEATURE_ORIC or feat['type'] == bc.FEATURE_ORIV]):5,}",
            f"{'oriTs':<30} :: {len([feat for feat in sequence_features if feat['type'] == bc.FEATURE_ORIT]):5,}",
            f"{'gaps':<30} :: {len([feat for feat in sequence_features if feat['type'] == bc.FEATURE_GAP]):5,}",
            f"{'pseudogenes':<30} :: {len([feat for feat in sequence_features if feat['type'] == bc.FEATURE_CDS and bc.PSEUDOGENE in feat]):5,}\n"
        )
        sequence_annotations = {
            'molecule_type': 'DNA',
            'source': data['genome'].get('taxon', ''),
            #'date': cfg.run_end.strftime('%d-%b-%Y').upper(),
            'topology': seq['topology'],
            'data_file_division': 'HGT' if seq['type'] == bc.REPLICON_CONTIG else 'BCT',
            # 'accession': '*',  # hold back until EMBL output bug is fixed in BioPython (https://github.com/biopython/biopython/pull/3572)
            'comment': comment
            # TODO: taxonomy
        }
        if euk:

            source_qualifiers = {}

            # Optional fields – add only if present
            optional_fields = [
                'mol_type',
                'organism', #  should be present
                'strain', #  should be present
                'isolate', #  should be present
                'db_xref', #  should be present
                'note', #  should be present
                'collection_date',
                'country',
                'host',
                'isolation_source'
            ]

            for field in optional_fields:
                value = seq.get(field)
                if value:  # only add if exists and is truth
                    source_qualifiers[field] = value
                else:
                    if field == 'mol_type':
                        # Always include mol_type
                        source_qualifiers['mol_type'] = 'genomic DNA'
            
        else:
            
            source_qualifiers = {
                'mol_type': 'genomic DNA'
                # 'molecule_type': 'DNA' #  might be necessary in BioPython > 1.78 along with removal of Seq(..., generic_dna)
            }

        description = ''
        if(data['genome'].get('taxon', None)):
            sequence_annotations['organism'] = data['genome']['taxon']
            source_qualifiers['organism'] = data['genome']['taxon']
            description = data['genome']['taxon']
        if(data['genome']['strain']):
            source_qualifiers['strain'] = data['genome']['strain']

        if(seq['type'] == bc.REPLICON_PLASMID):
            source_qualifiers['plasmid'] = seq['name'] if seq.get('name', None) else 'unnamed'
            description = f"{description} plasmid {seq.get('name', 'unnamed')}"
            description += ', complete sequence' if seq['complete'] else ', whole genome shotgun sequence'
        elif(seq['type'] == bc.REPLICON_CHROMOSOME):
            if seq.get('name', None):
                source_qualifiers['chromosome'] = seq['name']
            description = f'{description} chromosome, complete genome' if seq['complete'] else f"{description} chromosome {seq['id']}, whole genome shotgun sequence"
        else:
            description += f" {seq['id']}, whole genome shotgun sequence"

        if(len(description) > 0 and description[0] == ' '):  # discard potential leading whitespace
            description = description[1:]

        seq_bio = Seq(seq['nt']) if 'nt' in seq else Seq(seq['sequence'])  # <1.10.0 compatibility
        sequence_record = SeqIO.SeqRecord(id=seq['id'], name=seq['id'], description=description, annotations=sequence_annotations, seq=seq_bio)

        source = SeqFeature(FeatureLocation(0, seq['length'], strand=+1), type='source', qualifiers=source_qualifiers)
        seq_feature_list = [source]

        for feature in sequence_features:
            insdc_feature_type = None
            qualifiers = {
                'note': []
            }
            if('db_xrefs' in feature):
                if(cfg.compliant):
                    qualifiers['db_xref'], qualifiers['note'] = revise_dbxref_insdc(feature['db_xrefs'])
                else:
                    qualifiers['db_xref'] = feature['db_xrefs']
            if('product' in feature):
                qualifiers['product'] = feature['product']
            if('locus' in feature):
                qualifiers['locus_tag'] = feature['locus']
            if('standard_name' in feature): # euk mrnas
                qualifiers['standard_name'] = feature['standard_name'] 
            if euk:
                if('pseudo' in feature): # euk genes, mrna, CDS pseudo
                    qualifiers['pseudo'] = ['']        

            accompanying_features=[]
            if(feature['type'] == bc.FEATURE_GAP):
                insdc_feature_type = bc.INSDC_FEATURE_GAP
                qualifiers['estimated_length'] = feature['length']
            elif(feature['type'] == bc.FEATURE_ORIC or feature['type'] == bc.FEATURE_ORIV):
                # TODO: Add fuzzy positions for oriC/oriV
                insdc_feature_type = bc.INSDC_FEATURE_ORIGIN_REPLICATION
                qualifiers['inference'] = 'similar to DNA sequence'
                qualifiers['note'].append(feature['product'])
                if('product' in qualifiers):
                    qualifiers['note'] = feature['product']
                    del qualifiers['product']
            elif(feature['type'] == bc.FEATURE_ORIT):
                # TODO: Add fuzzy positions for oriT
                insdc_feature_type = bc.INSDC_FEATURE_ORIGIN_TRANSFER
                qualifiers['inference'] = 'similar to DNA sequence'
                qualifiers['note'].append(feature['product'])
                if('product' in qualifiers):
                    del qualifiers['product']

            elif(feature['type'] == bc.FEATURE_GENE):
                insdc_feature_type = bc.FEATURE_GENE
            elif(feature['type'] == bc.FEATURE_MRNA):
                insdc_feature_type = bc.FEATURE_MRNA
                

            elif(feature['type'] == bc.FEATURE_CDS) or (feature['type'] == bc.FEATURE_SORF):
                if(bc.PSEUDOGENE in feature):
                    qualifiers[bc.INSDC_FEATURE_PSEUDOGENE] = bc.INSDC_FEATURE_PSEUDOGENE_TYPE_UNPROCESSED if feature[bc.PSEUDOGENE]['paralog'] else bc.INSDC_FEATURE_PSEUDOGENE_TYPE_UNITARY
                    qualifiers['note'].append(feature[bc.PSEUDOGENE]['description'])
                else:
                    # add protein id for euk
                    if euk:
                        if feature.get('protein_id'):
                            qualifiers['protein_id'] = feature['protein_id']
                        else:
                            qualifiers['protein_id'] = f"gnl|Baktfold|{feature['locus']}"  
                    else:
                        qualifiers['protein_id'] = f"gnl|Baktfold|{feature['locus']}"
                    qualifiers['translation'] = feature['aa']
                qualifiers['codon_start'] = 1
                if not euk:
                    qualifiers['transl_table'] = cfg.translation_table
                insdc_feature_type = bc.INSDC_FEATURE_CDS
                inference = []
                if(feature['type'] == bc.FEATURE_CDS):
                    # add inference
                    if prokka:
                        inference.append('ab initio prediction:Prodigal:2.6')
                    if(feature.get('source', None) == bc.CDS_SOURCE_USER):
                        inference.append('EXISTENCE:non-experimental evidence, no additional details recorded')
                    else:
                        if not euk: # if euk just dont have anything
                            inference.append('ab initio prediction:Prodigal:2.6')
                else:
                    inference.append(f"ab initio prediction:Bakta:{'.'.join(cfg.version.split('.')[0:2])}")
                if('ncbi_nrp_id' in feature.get('ups', {})):
                    nrp_id = feature['ups']['ncbi_nrp_id']
                    inference.append(f'similar to AA sequence:{bc.DB_XREF_REFSEQ_NRP}:{nrp_id}')
                elif('uniref100_id' in feature.get('ips', {})):
                    ips_subject_id = feature['ips']['uniref100_id']
                    inference.append(f'similar to AA sequence:{bc.DB_XREF_UNIREF}:{ips_subject_id}')
                elif('uniparc_id' in feature.get('ups', {})):
                    uniparc_id = feature['ups']['uniparc_id']
                    inference.append(f'similar to AA sequence:{bc.DB_XREF_UNIPARC}:{uniparc_id}')
                elif('uniref90_id' in feature.get('psc', {}) and feature.get('psc', {}).get('valid', False)):
                    psc_subject_id = feature['psc']['uniref90_id']
                    inference.append(f'similar to AA sequence:{bc.DB_XREF_UNIREF}:{psc_subject_id}')
                elif('uniref50_id' in feature.get('pscc', {})):
                    pscc_subject_id = feature['pscc']['uniref50_id']
                    inference.append(f'similar to AA sequence:{bc.DB_XREF_UNIREF}:{pscc_subject_id}')
                qualifiers['inference'] = inference
                if(cfg.compliant):
                    qualifiers['note'], ec_number = extract_ec_from_notes_insdc(qualifiers, 'note')
                    if(ec_number is not None):
                            qualifiers['EC_number'] = ec_number
                if('exception' in feature):
                    ex = feature['exception']
                    pos = f"{ex['start']}..{ex['stop']}"
                    if(feature['strand'] == bc.STRAND_REVERSE):
                        pos = f"complement({pos})"
                    qualifiers['transl_except']=f"(pos:{pos},aa:{ex['aa']})"
                    qualifiers['note'].append(f"codon on position {ex['codon_position']} is a {ex['type']} codon")
                if(bc.FEATURE_SIGNAL_PEPTIDE in feature) and (not euk):
                    sigpep_qualifiers = {}
                    sigpep_qualifiers['locus_tag'] = feature['locus']
                    sigpep_qualifiers['inference'] = 'ab initio prediction:DeepSig:1.2'
                    sigpep = feature[bc.FEATURE_SIGNAL_PEPTIDE]
                    sigpep_strand = 1 if feature['strand'] == bc.STRAND_FORWARD else -1 if feature['strand'] == bc.STRAND_REVERSE else 0
                    sigpep_location = FeatureLocation(sigpep['start'] - 1, sigpep['stop'], strand=sigpep_strand)
                    sigpep_feature = SeqFeature(sigpep_location, type=bc.INSDC_FEATURE_SIGNAL_PEPTIDE, qualifiers=sigpep_qualifiers)
                    accompanying_features.append(sigpep_feature)
            elif(feature['type'] == bc.FEATURE_T_RNA):
                if('amino_acid' in feature and 'anti_codon' in feature):
                    if('anti_codon_pos' in feature):
                        anti_codon_pos = feature['anti_codon_pos']
                        qualifiers['anticodon'] = f"(pos:{anti_codon_pos[0]}..{anti_codon_pos[1]},aa:{feature['amino_acid']},seq:{feature['anti_codon']})"
                    else:
                        qualifiers['note'].append(f"tRNA-{feature['amino_acid']} ({feature['anti_codon']})")
                if prokka:
                    qualifiers['inference'] = 'profile:tRNAscan:2.0'
                elif euk:
                    qualifiers['inference'] = 'profile:other:unknown'
                else:
                    qualifiers['inference'] = 'profile:aragorn:1.2'
                insdc_feature_type = bc.INSDC_FEATURE_T_RNA
                if(bc.PSEUDOGENE in feature):
                    qualifiers[bc.INSDC_FEATURE_PSEUDOGENE] = bc.INSDC_FEATURE_PSEUDOGENE_TYPE_UNKNOWN
            elif(feature['type'] == bc.FEATURE_TM_RNA):
                qualifiers['inference'] = 'profile:aragorn:1.2'
                if('tag' in feature):
                    qualifiers['tag_peptide'] = f"{feature['tag']['start']}..{feature['tag']['stop']}"
                    if feature['strand'] == bc.STRAND_REVERSE:
                        qualifiers['tag_peptide'] = f"complement({qualifiers['tag_peptide']})"
                insdc_feature_type = bc.INSDC_FEATURE_TM_RNA
            elif(feature['type'] == bc.FEATURE_R_RNA):
                if prokka:
                    qualifiers['inference'] = 'profile:barrnap:0.9'
                for rfam_id in [dbxref.split(':')[1] for dbxref in feature['db_xrefs'] if dbxref.split(':')[0] == bc.DB_XREF_RFAM]:
                    qualifiers['inference'] = f'profile:Rfam:{rfam_id}'
                insdc_feature_type = bc.INSDC_FEATURE_R_RNA
            elif(feature['type'] == bc.FEATURE_NC_RNA):
                if prokka:
                    qualifiers['inference'] = 'profile:aragorn:1.2'
                for rfam_id in [dbxref.split(':')[1] for dbxref in feature['db_xrefs'] if dbxref.split(':')[0] == bc.DB_XREF_RFAM]:
                    qualifiers['inference'] = f'profile:Rfam:{rfam_id}'
                qualifiers[bc.INSDC_FEATURE_NC_RNA_CLASS] = select_ncrna_class(feature)
                insdc_feature_type = bc.INSDC_FEATURE_NC_RNA
            elif(feature['type'] == bc.FEATURE_NC_RNA_REGION):
                for rfam_id in [dbxref.split(':')[1] for dbxref in feature['db_xrefs'] if dbxref.split(':')[0] == bc.DB_XREF_RFAM]:
                    qualifiers['inference'] = f'profile:Rfam:{rfam_id}'
                qualifiers[bc.INSDC_FEATURE_REGULATORY_CLASS] = select_regulatory_class(feature)
                insdc_feature_type = bc.INSDC_FEATURE_REGULATORY
                qualifiers['note'].append(feature['product'])
                qualifiers.pop('product', None)
            elif(feature['type'] == bc.FEATURE_CRISPR):
                qualifiers[bc.INSDC_FEATURE_REPEAT_FAMILY] = 'CRISPR'
                qualifiers[bc.INSDC_FEATURE_REPEAT_TYPE] = 'direct'
                # prokka wont have this 
                if('repeat_consensus' in feature):
                    qualifiers[bc.INSDC_FEATURE_REPEAT_UNIT_SEQ] = feature['repeat_consensus']
                    qualifiers['inference'] = 'COORDINATES:alignment:pilercr:1.02'
                # prokka uses minced - maybe for later here
                insdc_feature_type = bc.INSDC_FEATURE_REPEAT_REGION
                qualifiers['note'].append(feature['product'])
                qualifiers.pop('product', None)
            elif(feature['type'] == bc.FEATURE_REPEAT) and (euk): # for euks
                insdc_feature_type = bc.INSDC_FEATURE_REPEAT_REGION
                qualifiers['rpt_family'] = feature['family']
                qualifiers['note'].append(feature['product'])
                qualifiers.pop('product', None)
            elif(feature['type'] == bc.FEATURE_5UTR) and (euk): # for euks
                insdc_feature_type = bc.INSDC_FEATURE_5UTR
                move_product_to_note_if_exists(qualifiers)
            elif(feature['type'] == bc.FEATURE_3UTR) and (euk): # for euks
                insdc_feature_type = bc.INSDC_FEATURE_3UTR
                move_product_to_note_if_exists(qualifiers)
            elif(feature['type'] == bc.FEATURE_MISC_RNA) and (euk): # for euks
                insdc_feature_type = bc.INSDC_FEATURE_MISC_RNA
                move_product_to_note_if_exists(qualifiers)
            elif(feature['type'] == bc.FEATURE_EXON) and (euk): # for euks
                insdc_feature_type = bc.INSDC_FEATURE_EXON
                move_product_to_note_if_exists(qualifiers)
            elif(feature['type'] == bc.FEATURE_MAT_PEPTIDE) and (euk): # for euks
                insdc_feature_type = bc.INSDC_FEATURE_MAT_PEPTIDE
                move_product_to_note_if_exists(qualifiers)
            elif(feature['type'] == bc.FEATURE_MOBILE_ELEMENT) and (euk): # for euks
                insdc_feature_type = bc.INSDC_FEATURE_MOBILE_ELEMENT
                move_product_to_note_if_exists(qualifiers)
            elif(feature['type'] == bc.INSDC_FEATURE_MISC_FEATURE) and (euk): # for euks
                insdc_feature_type = bc.INSDC_FEATURE_MISC_FEATURE
                move_product_to_note_if_exists(qualifiers)
            elif(feature['type'] == bc.FEATURE_PRECURSOR_RNA) and (euk): # for euks
                insdc_feature_type = bc.INSDC_FEATURE_PRECURSOR_RNA
                move_product_to_note_if_exists(qualifiers)          
            elif(feature['type'] == bc.FEATURE_PROPEPTIDE) and (euk): # for euks
                insdc_feature_type = bc.INSDC_FEATURE_PROPEPTIDE
                move_product_to_note_if_exists(qualifiers)
            elif(feature['type'] == bc.FEATURE_PROTEIN_BIND) and (euk): # for euks
                insdc_feature_type = bc.INSDC_FEATURE_PROTEIN_BIND
                move_product_to_note_if_exists(qualifiers)
            elif(feature['type'] == bc.INSDC_FEATURE_REGULATORY) and (euk): # for euks
                insdc_feature_type = bc.INSDC_FEATURE_REGULATORY
                qualifiers['note'].append(feature['product']) 
                qualifiers.pop('product', None)  
            elif(feature['type'] == bc.FEATURE_SIGNAL_PEPTIDE) and (euk): # for euks
                insdc_feature_type = bc.INSDC_FEATURE_SIGNAL_PEPTIDE
                move_product_to_note_if_exists(qualifiers)
            elif(feature['type'] == bc.FEATURE_TRANSIT_PEPTIDE): # for euks
                insdc_feature_type = bc.INSDC_FEATURE_TRANSIT_PEPTIDE
                move_product_to_note_if_exists(qualifiers)           
            

            strand = None
            if(feature['strand'] == bc.STRAND_FORWARD):
                strand = 1
            elif(feature['strand'] == bc.STRAND_REVERSE):
                strand = -1
            elif(feature['strand'] == bc.STRAND_UNKNOWN):
                strand = 0

            # euk
            if (
                feature.get("starts") 
                and feature.get("stops") 
                and isinstance(feature["starts"], list) 
                and isinstance(feature["stops"], list)
                and len(feature["starts"]) == len(feature["stops"])
            ):
                # Multi-exon CDS/mRNA – use the FIRST and LAST coords as overall range
                start = feature["starts"][0] - 1       # 0-based for Biopython
                stop = feature["stops"][-1]            # inclusive

                # also build in here the parts

                coords = list(zip(feature.get("starts"), feature.get("stops")))  # multi-part

                # Build Biopython FeatureLocation or CompoundLocation
                strand = 1 if feature['strand'] == bc.STRAND_FORWARD else -1
                parts = []

                for s, e in coords:
                    parts.append(
                        FeatureLocation(s - 1, e, strand=strand)
                    )

                feature_location = (
                    CompoundLocation(parts)
                    if len(parts) > 1
                    else parts[0]
                )


            else:

                start = feature['start'] - 1
                stop = feature['stop']
                
                # will be for all non euk ones here
                if('edge' in feature):
                    fl_1 = FeatureLocation(start, seq['length'], strand=strand)
                    fl_2 = FeatureLocation(0, stop, strand=strand)
                    if(feature['strand'] == bc.STRAND_REVERSE):
                        feature_location = CompoundLocation([fl_2, fl_1])
                    else:
                        feature_location = CompoundLocation([fl_1, fl_2])
                else:
                    if('truncated' in feature):
                        if(bc.PSEUDOGENE not in feature):  # only add /pseudo qualifier if /pseudogene is not already set
                            qualifiers[bc.INSDC_FEATURE_PSEUDO] = None
                        if(feature['truncated'] == bc.FEATURE_END_5_PRIME):
                            qualifiers['note'].append("(5' truncated)")
                            if(feature['strand'] == bc.STRAND_FORWARD):
                                start = BeforePosition(start)
                            else:
                                stop = AfterPosition(stop)
                        elif(feature['truncated'] == bc.FEATURE_END_3_PRIME):
                            qualifiers['note'].append("(3' truncated)")
                            if(feature['strand'] == bc.STRAND_FORWARD):
                                stop = AfterPosition(stop)
                            else:
                                start = BeforePosition(start)
                        elif(feature['truncated'] == bc.FEATURE_END_BOTH):
                            qualifiers['note'].append('(partial)')
                            start = BeforePosition(start)
                            stop = AfterPosition(stop)
                    feature_location = FeatureLocation(start, stop, strand=strand)
            if(feature.get('locus', None)):
                gene_qualifier = {
                    'locus_tag': feature['locus']
                }
                if(feature.get('gene', None)):
                    if(cfg.compliant):
                        if(ba.RE_GENE_SYMBOL.fullmatch(feature['gene'])):  # discard non-standard gene symbols
                            gene_qualifier['gene'] = feature['gene']
                        else:
                            qualifiers.pop('gene', None)
                    else:
                        qualifiers['gene'] = feature['gene']
                        gene_qualifier['gene'] = feature['gene']
                if(bc.PSEUDOGENE in feature):
                    if(feature['type'] == bc.FEATURE_CDS):
                        gene_qualifier[bc.INSDC_FEATURE_PSEUDOGENE] = bc.INSDC_FEATURE_PSEUDOGENE_TYPE_UNPROCESSED if feature[bc.PSEUDOGENE]['paralog'] else bc.INSDC_FEATURE_PSEUDOGENE_TYPE_UNITARY
                    else:
                        gene_qualifier[bc.INSDC_FEATURE_PSEUDOGENE] = bc.INSDC_FEATURE_PSEUDOGENE_TYPE_UNKNOWN
                elif('truncated' in feature):
                    gene_qualifier[bc.INSDC_FEATURE_PSEUDO] = None
                gen_seqfeat = SeqFeature(feature_location, type='gene', qualifiers=gene_qualifier)
                # dont append gene feature for CDS MRNA GENE euks
                if not euk:
                    seq_feature_list.append(gen_seqfeat)
            feat_seqfeat = SeqFeature(feature_location, type=insdc_feature_type, qualifiers=qualifiers)
            seq_feature_list.append(feat_seqfeat)
            for acc_feature in accompanying_features:  # add accompanying features, e.g. signal peptides
                seq_feature_list.append(acc_feature)
        sequence_record.features = seq_feature_list
        sequence_list.append(sequence_record)
    return sequence_list


def write_features(data: dict, features: Sequence[dict], genbank_output_path: Path, embl_output_path: Path, prokka: bool = False, euk: bool = False):
    logger.info(f'prepare: genbank={genbank_output_path}, embl={embl_output_path}')

    # fix later
    sequence_list = build_biopython_sequence_list(data, features, prokka, euk)
    
    with genbank_output_path.open('wt', encoding='utf-8') as fh:
        logger.info(f'write GenBank: path={genbank_output_path}')
        SeqIO.write(sequence_list, fh, format='genbank')

    with embl_output_path.open('wt', encoding='utf-8') as fh:
        logger.info(f'write EMBL: path={embl_output_path}')
        SeqIO.write(sequence_list, fh, format='embl')


def select_ncrna_class(feature: dict) -> str:
    feature_class = feature['class']
    if(feature_class is None):
        return bc.INSDC_FEATURE_NC_RNA_CLASS_OTHER
    else:
        if(isinstance(feature_class, list)):  #  workaround for JSON-imported features
            feature_class = so.SO(feature_class[0], feature_class[1])
            feature['class'] = feature_class

        if(feature_class.id == so.SO_NCRNA_GENE_ANTISENSE.id):
            return bc.INSDC_FEATURE_NC_RNA_CLASS_ANTISENSE
        elif(feature_class.id == so.SO_NCRNA_GENE_RIBOZYME.id):
            return bc.INSDC_FEATURE_NC_RNA_CLASS_RIBOZYME
        elif(feature_class.id == so.SO_NCRNA_GENE_RNASEP.id):
            return bc.INSDC_FEATURE_NC_RNA_CLASS_RNASEP
        else:
            return bc.INSDC_FEATURE_NC_RNA_CLASS_OTHER


def select_regulatory_class(feature: dict) -> str:
    feature_class = feature['class']
    if(feature_class is None):
        return bc.INSDC_FEATURE_REGULATORY_CLASS_OTHER
    else:
        if(isinstance(feature_class, list)):  #  workaround for JSON-imported features
            feature_class = so.SO(feature_class[0], feature_class[1])
            feature['class'] = feature_class

        if(feature_class.id == so.SO_CIS_REG_ATTENUATOR.id):
            return bc.INSDC_FEATURE_REGULATORY_CLASS_ATTENUATOR
        elif(feature_class.id == so.SO_CIS_REG_RIBOSWITCH.id):
            return bc.INSDC_FEATURE_REGULATORY_CLASS_RIBOSWITCH
        elif(feature_class.id == so.SO_CIS_REG_THERMOMETER.id):
            return bc.INSDC_FEATURE_REGULATORY_CLASS_RESPONSE_ELEMENT
        elif(feature_class.id == so.SO_CIS_REG_RECODING_STIMULATION_REGION.id or feature_class.id == so.SO_CIS_REG_FRAMESHIFT.id):
            return bc.INSDC_FEATURE_REGULATORY_CLASS_RECODING_STIMULATORY_REGION
        elif(feature_class.id == so.SO_CIS_REG_RIBOSOME_BINDING_SITE.id):
            return bc.INSDC_FEATURE_REGULATORY_CLASS_RIBOSOME_BINDING_SITE
        else:
            return bc.INSDC_FEATURE_REGULATORY_CLASS_OTHER


def revise_product_insdc(product: str):
    """Revise product name for INSDC compliant submissions"""

    old_product = product
    if(re.search(r'(uncharacteri[sz]ed)', product, flags=re.IGNORECASE)):  # replace putative synonyms)
        product = re.sub(r'(uncharacteri[sz]ed)', 'putative', product, flags=re.IGNORECASE)
        logger.info('fix product: replace putative synonyms. new=%s, old=%s', product, old_product)

    old_product = product
    if(product.count('(') != product.count(')')):  # remove unbalanced parentheses
        product = product.replace('(', '').replace(')', '')  # ToDo: find and replace only legend parentheses
        logger.info('fix product: remove unbalanced parantheses. new=%s, old=%s', product, old_product)

    old_product = product
    if(product.count('[') != product.count(']')):  # remove unbalanced brackets
        product = product.replace('[', '').replace(']', '')  # ToDo: find and replace only legend bracket
        logger.info('fix product: remove unbalanced brackets. new=%s, old=%s', product, old_product)

    return product


def revise_dbxref_insdc(dbxrefs: Sequence[str]) -> Tuple[Sequence[str], Sequence[str]]:
    """Remove INSDC non-compliant DbXrefs."""
    insdc_valid_dbxrefs = [bc.DB_XREF_UNIPROTKB, bc.DB_XREF_GO, bc.DB_XREF_PFAM, bc.DB_XREF_RFAM]
    valid_dbxrefs = []
    invalid_dbxrefs = []
    for dbxref in dbxrefs:
        if(dbxref.split(':')[0] in insdc_valid_dbxrefs):
            valid_dbxrefs.append(dbxref)
        else:
            invalid_dbxrefs.append(dbxref)
    return valid_dbxrefs, invalid_dbxrefs


def extract_ec_from_notes_insdc(qualifiers: dict, note_key: str):
    note_tmp = []
    ec_number = None
    for note in qualifiers[note_key]:  # move EC numbers from note to EC_number qualifier
        if(note.split(':')[0] == bc.DB_XREF_EC):
            ec_number = note.replace('EC:', '')
        else:
            note_tmp.append(note)
    return note_tmp, ec_number
