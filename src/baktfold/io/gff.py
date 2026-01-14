# import logging

from pathlib import Path
from typing import Dict, Sequence, Union

from loguru import logger

import baktfold.bakta
import baktfold.bakta.config as cfg
import baktfold.bakta.constants as bc
import baktfold.io.fasta as fasta
import baktfold.io.insdc as insdc
import baktfold.bakta.annotation as ba
import baktfold.bakta.so as so


def write_gene_feature(fh, seq_id, feat):
    """Write a 'gene' feature including fuzzy boundaries."""
    start = int(feat['start'])
    stop  = int(feat['stop'])
    strand = feat['strand']

    locus = feat['locus']

    attrs = {
        "ID": f"{locus}"
    }

    if feat.get('gene') is not None:
        attrs["Name"] = feat.get('gene')

    attr_str = ";".join(f"{k}={v}" for k, v in attrs.items())

    fh.write(f"{seq_id}\tbaktfold\tgene\t{start}\t{stop}\t.\t{strand}\t.\t{attr_str}\n")


def write_mrna_feature(fh, seq_id, feat):
    """Write mRNA + implied exons based on join() structure."""

    start = int(feat['start'])
    stop  = int(feat['stop'])
    strand = feat['strand']

    locus = feat.get("locus")

    mrna_id = f"{locus}-T1"

    # Top-level mRNA line
    attrs = {
        "ID": mrna_id,
        "Parent": f"{locus}",
    }

    product = feat.get("product", [])

    if product:

        key = "product"         
        if isinstance(product, list):
            if len(product) == 1:
                attrs[key] = str(product[0])
            else:
                attrs[key] = ",".join(str(v) for v in product)
        else:
            attrs[key] = str(product)


    # Ensure db_xrefs exists and is a list
    db_xrefs = feat.get("db_xrefs", [])

    # Access note safely
    note = feat.get("note", None)


    if db_xrefs:

        key = "Dbxref"         
        if isinstance(db_xrefs, list):
            if len(db_xrefs) == 1:
                attrs[key] = str(db_xrefs[0])
            else:
                attrs[key] = ",".join(str(v) for v in db_xrefs)
        else:
            # if somehow not a list, just convert to string
            attrs[key] = str(db_xrefs)

    if note:

        key = "note"         # <-- you must define this
        if isinstance(db_xrefs, list):
            if len(db_xrefs) == 1:
                attrs[key] = str(db_xrefs[0])
            else:
                attrs[key] = ",".join(str(v) for v in db_xrefs)
        else:
            # if somehow not a list, just convert to string
            attrs[key] = str(db_xrefs)

    
    attr_str = ";".join(f"{k}={v}" for k, v in attrs.items())

    fh.write(f"{seq_id}\tbaktfold\tmRNA\t{start}\t{stop}\t.\t{strand}\t.\t{attr_str}\n")

    starts = feat.get("starts")
    stops  = feat.get("stops")
    strand = feat.get("strand")
    seq_id = feat.get("sequence")

    if (
        isinstance(starts, list)
        and isinstance(stops, list)
        and len(starts) == len(stops)
        and len(starts) > 0
    ):
        # For minus strand, exons must be written in reverse order (5'→3')
        if strand == "-":
            exon_parts = list(zip(starts, stops))
        else:
            exon_parts = list(zip(starts, stops))

        # Exons must be numbered in biological order (5' to 3')
        if strand == "-":
            exon_parts = exon_parts[::-1]   # reverse order

        # Write each exon to GFF
        for idx, (ex_start, ex_stop) in enumerate(exon_parts, start=1):
            exon_id = f"{mrna_id}.exon{idx}"
            exon_attrs = f"ID={exon_id};Parent={mrna_id}"
            fh.write(
                f"{seq_id}\tbaktfold\texon\t{ex_start}\t{ex_stop}\t.\t{strand}\t.\t{exon_attrs}\n"
            )
    else:
        # Single exon (no starts/stops provided)
        exon_start = feat["start"]
        exon_stop = feat["stop"]
        exon_id = f"{mrna_id}.exon1"
        exon_attrs = f"ID={exon_id};Parent={mrna_id}"

        fh.write(
            f"{seq_id}\tbaktfold\texon\t{exon_start}\t{exon_stop}"
            f"\t.\t{feat['strand']}\t.\t{exon_attrs}\n"
        )



def write_euk_cds_feature(fh, seq_id, feat):
    """
    Write a eukaryotic CDS feature to GFF3 with multiple CDS parts.

    Parameters
    ----------
    fh : file-handle
    seq_id : str
    feat : dict-like feature with keys:
            "start", "stop", "strand", "locus", "starts", "stops"
    """

    strand = feat.get("strand", "+")
    locus = feat.get("locus", "unknown")

    transcript_id = f"{locus}-T1"
    cds_id = f"{transcript_id}.cds"

    starts = feat.get("starts")
    stops = feat.get("stops")

    # -------------------------------
    # 1. Determine CDS sub-coordinates
    # -------------------------------
    if (
        isinstance(starts, list)
        and isinstance(stops, list)
        and len(starts) == len(stops)
        and len(starts) > 0
    ):
        cds_coords = list(zip(starts, stops))
    else:
        cds_coords = [(feat["start"], feat["stop"])]

    # -------------------------------
    # 2. Reverse order for negative strand
    # -------------------------------
    if strand == "-":
        cds_coords.reverse()

    # -------------------------------
    # 3. Emit CDS lines with correct phase
    # -------------------------------
    offset = 0

    for i, (cds_start, cds_stop) in enumerate(cds_coords, start=1):

        length = cds_stop - cds_start + 1
        phase = offset % 3
        offset += length

        attr = f"ID={cds_id}-{i};Parent={transcript_id}"

        fh.write(
            f"{seq_id}\tbaktfold\tCDS\t{cds_start}\t{cds_stop}"
            f"\t.\t{strand}\t{phase}\t{attr}\n"
        )

def write_euk_repeat_region_feature(fh, seq_id, feat):
    """
    Writes a repeat region feature to a file.

    Args:
      fh (file): The file handle to write to.
      seq_id (str): The sequence ID.
      feat (dict): A dictionary containing the feature information.

    Returns:
      None

    Examples:
      >>> write_euk_repeat_region_feature(fh, 'DS572673.1', {
          "type": "repeat_region",
          "sequence": "DS571531.1",
          "start": 1470,
          "stop": 1716,
          "strand": "?",
          "family": "LINE2",
          "rpt_type": null,
          "repeat_unit": null,
          "product": null,
          "nt": "AATAAAATCATATCAGAAATAAAAAGAATGAAAATAAACAAATTAAAGAAAATAATTATAAAATTAATAAACGATATTTAAATGAAAGAAAATAGAGAATATGTAATAAGTACAAATGGTTCATTCATTAATAAGAAATTAACAATAATAAAATAGAGAATATTGATTATAAAAAGAAATATATTTCTCAAAACAGTAGAGATACAAAAAGAATAGATATGAAATAAATATTAATTCTAAAATACTC",
          "id": "EHICP_3230",
          "db_xrefs": [
              "SO:0000657"
          ]
      })
    """

    start = int(feat['start'])
    stop  = int(feat['stop'])
    strand = feat['strand']

    id = feat['sequence']

    attrs = {
        "ID": f"{id}:{start}..{stop}",
        "gbkey": "repeat_region"
    }

    if feat.get('family') is not None:
        attrs["rpt_family"] = feat.get('family')

    attr_str = ";".join(f"{k}={v}" for k, v in attrs.items())

    fh.write(f"{seq_id}\tbaktfold\trepeat_region\t{start}\t{stop}\t.\t{strand}\t.\t{attr_str}\n")

    # DS572673.1	Genbank	repeat_region	1	364	.	-	.	ID=id-DS572673.1:1..364;gbkey=repeat_region;rpt_family=LINE1

        #     "type": "repeat_region",
        #     "sequence": "DS571531.1",
        #     "start": 1470,
        #     "stop": 1716,
        #     "strand": "?",
        #     "family": "LINE2",
        #     "rpt_type": null,
        #     "repeat_unit": null,
        #     "product": null,
        #     "nt": "AATAAAATCATATCAGAAATAAAAAGAATGAAAATAAACAAATTAAAGAAAATAATTATAAAATTAATAAACGATATTTAAATGAAAGAAAATAGAGAATATGTAATAAGTACAAATGGTTCATTCATTAATAAGAAATTAACAATAATAAAATAGAGAATATTGATTATAAAAAGAAATATATTTCTCAAAACAGTAGAGATACAAAAAGAATAGATATGAAATAAATATTAATTCTAAAATACTC",
        #     "id": "EHICP_3230",
        #     "db_xrefs": [
        #         "SO:0000657"
        #     ]
        # },

def write_euk_utr_feature(fh, seq_id, feat, locus_counter, three=False):
    """Write a 'utr' feature."""
    start = int(feat['start'])
    stop  = int(feat['stop'])
    strand = feat['strand']

    locus = feat['locus']

    # Count occurrences for this locus
    count = locus_counter.get(locus, 0) + 1
    locus_counter[locus] = count

    # Construct ID with suffix -2, -3, etc.
    # For first entry we keep ID=locus (no -1)
    if count == 1:
        utr_id = locus
    else:
        utr_id = f"{locus}-{count}"

    # Top-level mRNA line
    attrs = {
        "ID": f"{utr_id}",
        "Parent": f"{locus}",
    }

# CAMXCT020000566.1	EMBL	three_prime_UTR	84568	84617	.	-	.	ID=id-C1SCF055_LOCUS8420;Parent=gene-C1SCF055_LOCUS8420;Note=ID:SCF055_s1507_g28601.utr3p1%3B~source:feature;gbkey=3'UTR;locus_tag=C1SCF055_LOCUS8420
# CAMXCT020000566.1	EMBL	five_prime_UTR	136251	136259	.	-	.	ID=id-C1SCF055_LOCUS8420-2;Parent=gene-C1SCF055_LOCUS8420;Note=ID:SCF055_s1507_g28601.utr5p1%3B~source:feature;gbkey=5'UTR;locus_tag=C1SCF055_LOCUS8420

    if feat.get('Note') is not None:
        attrs["Note"] = feat.get('note')


    attrs["gbkey"] = "3'UTR" if three else "5'UTR"
        
    if feat.get('Note') is not None:
        attrs["locus_tag"] = feat.get('locus')

    attr_str = ";".join(f"{k}={v}" for k, v in attrs.items())

    if three:
        gene_tag = 'three_prime_UTR'
    else:
        gene_tag = 'five_prime_UTR'

    fh.write(f"{seq_id}\tbaktfold\t{gene_tag}\t{start}\t{stop}\t.\t{strand}\t.\t{attr_str}\n")


def write_euk_trna_feature(fh, seq_id, feat):
    """
    Write a tRNA feature to GFF3 with a top-level line and single exon.
    
    Parameters
    ----------
    fh : file-like
        Open file handle to write GFF lines.
    seq_id : str
        Sequence/contig ID.
    feat : SeqFeature
        Biopython SeqFeature object of type 'tRNA'.

    Notes
    -----
    - Generates one tRNA line and one exon line.
    - Includes optional 'product' qualifier.
    """
    start = int(feat['start'])
    stop  = int(feat['stop'])
    
    strand = feat['strand']
    
    locus = feat['locus']

    trna_id = f"{locus}-T1"

    # Top-level tRNA attributes
    attrs = {
        "ID": trna_id,
        "Parent": locus
    }

    attrs = {}

    product = feat.get("product", [])

    if product:

        key = "product"         
        if isinstance(product, list):
            if len(product) == 1:
                attrs[key] = str(product[0])
            else:
                attrs[key] = ",".join(str(v) for v in product)
        else:
            attrs[key] = str(product)


    attr_str = ";".join(f"{k}={v}" for k, v in attrs.items())

    # Write top-level tRNA line
    fh.write(f"{seq_id}\tbaktfold\ttRNA\t{start}\t{stop}\t.\t{strand}\t.\t{attr_str}\n")

    # Write exon line (tRNA single-exon)
    exon_id = f"{trna_id}.exßon1"
    exon_attrs = f"ID={exon_id};Parent={trna_id}"
    fh.write(f"{seq_id}\tbaktfold\texon\t{start}\t{stop}\t.\t{strand}\t.\t{exon_attrs}\n")


def write_features(data: dict, features_by_sequence: Dict[str, dict], gff3_path: Path, prokka: bool = False, euk: bool = False):
    """Export features in GFF3 format."""
    logger.info(f'write features: path={gff3_path}')

    with gff3_path.open('wt') as fh:
        fh.write('##gff-version 3\n')  # GFF version
        fh.write('##feature-ontology https://github.com/The-Sequence-Ontology/SO-Ontologies/blob/v3.1/so.obo\n')  # SO feature version

        if(data['genome'].get('taxon', None)):  # write organism info
            fh.write(f"# organism {data['genome']['taxon']}\n")

        fh.write('# Annotated with Baktfold\n')
        fh.write(f'# Software: v{cfg.version}\n')
        fh.write(f"# Database: v{cfg.version}\n") # fix later
        #fh.write(f"# Database: v{cfg.db_info['major']}.{cfg.db_info['minor']}, {cfg.db_info['type']}\n")
        fh.write(f'# DOI: {bc.BAKTFOLD_DOI}\n')
        fh.write(f'# URL: {bc.BAKTFOLD_URL}\n')

        for seq in data['sequences']:  # write features
            if euk:
                locus_counter = {} # for UTRs
                
            fh.write(f"##sequence-region {seq['id']} 1 {seq['length']}\n")  # sequence region

            # write landmark region
            annotations = {
                'ID': seq['id'],
                'Name': seq['id']
            }
            if(seq['topology'] == bc.TOPOLOGY_CIRCULAR):
                annotations['Is_circular'] = 'true'
            annotations = encode_annotations(annotations)
            fh.write(f"{seq['id']}\tBaktfold\tregion\t1\t{str(seq['length'])}\t.\t+\t.\t{annotations}\n")

            for feat in features_by_sequence[seq['id']]:
                seq_id = feat['sequence'] if 'sequence' in feat else feat['contig']  # <1.10.0 compatibility
                start = feat['start']
                stop = feat['stop']
                if('edge' in feat):
                    stop += seq['length']

                # euks
                if euk:
                    if(feat['type'] == bc.FEATURE_REPEAT):
                        write_euk_repeat_region_feature(fh, seq_id, feat)

                    if(feat['type'] == bc.FEATURE_5UTR or feat['type'] == bc.FEATURE_3UTR):
                        if feat['type'] == bc.FEATURE_3UTR:
                            write_euk_utr_feature(fh, seq_id, feat, locus_counter, three=True)
                        elif feat['type'] == bc.FEATURE_5UTR:
                            write_euk_utr_feature(fh, seq_id, feat, locus_counter, three=False)

                if(feat['type'] == bc.FEATURE_T_RNA):

                    if euk:
                        write_euk_trna_feature(fh, seq_id, feat)
                    else:

                        trna_tool = "tRNAscan-SE"
                        if prokka:
                            trna_tool = "Aragorn"
                        annotations = {
                            'ID': feat['locus'],
                            'Name': feat['product'],
                            'locus_tag': feat['locus'],
                            'product': feat['product'],
                            'Dbxref': feat['db_xrefs']
                        }
                        if(feat.get('gene', None)):  # add gene annotation if available
                            annotations['gene'] = feat['gene']
                        if(bc.PSEUDOGENE in feat):
                            annotations[bc.INSDC_FEATURE_PSEUDOGENE] = bc.INSDC_FEATURE_PSEUDOGENE_TYPE_UNKNOWN
                        elif('truncated' in feat):
                            annotations[bc.INSDC_FEATURE_PSEUDO] = True
                        if(feat.get('anti_codon', False)):
                            annotations['anti_codon'] = feat['anti_codon']
                        if(feat.get('amino_acid', False)):
                            annotations['amino_acid'] = feat['amino_acid']
                        if(cfg.compliant):
                            gene_id = f"{feat['locus']}_gene"
                            annotations['Parent'] = gene_id
                            annotations['inference'] = 'profile:tRNAscan:2.0'
                            annotations['Dbxref'], annotations['Note'] = insdc.revise_dbxref_insdc(feat['db_xrefs'])  # remove INSDC invalid DbXrefs
                            gene_annotations = {
                                'ID': gene_id,
                                'locus_tag': feat['locus']
                            }
                            if(feat.get('gene', None)):
                                gene_annotations['gene'] = feat['gene']
                            if(bc.PSEUDOGENE in feat):
                                gene_annotations[bc.INSDC_FEATURE_PSEUDOGENE] = bc.INSDC_FEATURE_PSEUDOGENE_TYPE_UNKNOWN
                            gene_annotations = encode_annotations(gene_annotations)
                            fh.write(f"{seq_id}\t{trna_tool}\tgene\t{start}\t{stop}\t.\t{feat['strand']}\t.\t{gene_annotations}\n")
                        annotations = encode_annotations(annotations)
                        fh.write(f"{seq_id}\t{trna_tool}\t{so.SO_TRNA.name}\t{start}\t{stop}\t.\t{feat['strand']}\t.\t{annotations}\n")
                elif(feat['type'] == bc.FEATURE_TM_RNA):
                    # both prokka and bakta use Aragorn
                    annotations = {
                        'ID': feat['locus'],
                        'Name': feat['product'],
                        'locus_tag': feat['locus'],
                        'gene': feat['gene'],
                        'product': feat['product'],
                        'Dbxref': feat['db_xrefs']
                    }
                    if('tag' in feat):
                        annotations['tag_peptide'] = feat['tag']['aa']
                    if('truncated' in feat):
                        annotations[bc.INSDC_FEATURE_PSEUDO] = True
                    if(cfg.compliant):
                        gene_id = f"{feat['locus']}_gene"
                        annotations['Parent'] = gene_id
                        annotations['inference'] = 'profile:aragorn:1.2'
                        annotations['Dbxref'], annotations['Note'] = insdc.revise_dbxref_insdc(feat['db_xrefs'])  # remove INSDC invalid DbXrefs
                        if('tag' in feat):
                            annotations['tag_peptide'] = f"{feat['tag']['start']}..{feat['tag']['stop']}" if feat['strand'] == bc.STRAND_FORWARD else f"complement({feat['tag']['start']}..{feat['tag']['stop']})"
                        gene_annotations = {
                            'ID': gene_id,
                            'locus_tag': feat['locus'],
                            'gene': feat['gene']
                        }
                        if('truncated' in feat):
                            gene_annotations[bc.INSDC_FEATURE_PSEUDO] = True
                        gene_annotations = encode_annotations(gene_annotations)
                        fh.write(f"{seq_id}\tAragorn\tgene\t{start}\t{stop}\t.\t{feat['strand']}\t.\t{gene_annotations}\n")
                    annotations = encode_annotations(annotations)
                    fh.write(f"{seq_id}\tAragorn\t{so.SO_TMRNA.name}\t{start}\t{stop}\t.\t{feat['strand']}\t.\t{annotations}\n")
                elif(feat['type'] == bc.FEATURE_R_RNA):
                    rrna_tool = "Infernal"
                    if prokka:
                        rrna_tool = "barrnap"
                    annotations = {
                        'ID': feat['locus'],
                        'Name': feat['product'],
                        'locus_tag': feat['locus'],
                        'gene': feat['gene'],
                        'product': feat['product'],
                        'Dbxref': feat['db_xrefs']
                    }
                    if('truncated' in feat):
                        annotations[bc.INSDC_FEATURE_PSEUDO] = True
                    if(cfg.compliant):
                        gene_id = f"{feat['locus']}_gene"
                        annotations['Parent'] = gene_id
                        annotations['Dbxref'], annotations['Note'] = insdc.revise_dbxref_insdc(feat['db_xrefs'])  # remove INSDC invalid DbXrefs
                        for rfam_id in [dbxref.split(':')[1] for dbxref in feat['db_xrefs'] if dbxref.split(':')[0] == bc.DB_XREF_RFAM]:
                            annotations['inference'] = f'profile:Rfam:{rfam_id}'
                        gene_annotations = {
                            'ID': gene_id,
                            'locus_tag': feat['locus'],
                            'gene': feat['gene']
                        }
                        if('truncated' in feat):
                            gene_annotations[bc.INSDC_FEATURE_PSEUDO] = True
                        gene_annotations = encode_annotations(gene_annotations)
                        fh.write(f"{seq_id}\t{rrna_tool}\tgene\t{start}\t{stop}\t.\t{feat['strand']}\t.\t{gene_annotations}\n")
                    annotations = encode_annotations(annotations)
                    fh.write(f"{seq_id}\t{rrna_tool}\t{so.SO_RRNA.name}\t{start}\t{stop}\t{feat['evalue']}\t{feat['strand']}\t.\t{annotations}\n")
                elif(feat['type'] == bc.FEATURE_NC_RNA):
                    # both prokka and bakta use infernal for ncrna
                    annotations = {
                        'ID': feat['locus'],
                        'Name': feat['product'],
                        'locus_tag': feat['locus'],
                        'gene': feat['gene'],
                        'product': feat['product'],
                        'Dbxref': feat['db_xrefs']
                    }
                    if('truncated' in feat):
                        annotations[bc.INSDC_FEATURE_PSEUDO] = True
                    if(cfg.compliant):
                        gene_id = f"{feat['locus']}_gene"
                        annotations['Parent'] = gene_id
                        annotations['Dbxref'], annotations['Note'] = insdc.revise_dbxref_insdc(feat['db_xrefs'])  # remove INSDC invalid DbXrefs
                        annotations[bc.INSDC_FEATURE_NC_RNA_CLASS] = insdc.select_ncrna_class(feat)
                        for rfam_id in [dbxref.split(':')[1] for dbxref in feat['db_xrefs'] if dbxref.split(':')[0] == bc.DB_XREF_RFAM]:
                            annotations['inference'] = f'profile:Rfam:{rfam_id}'
                        gene_annotations = {
                            'ID': gene_id,
                            'locus_tag': feat['locus'],
                            'gene': feat['gene']
                        }
                        if(ba.RE_GENE_SYMBOL.fullmatch(feat['gene'])):  # discard non-standard ncRNA gene symbols
                            gene_annotations['gene'] = feat['gene']
                        else:
                            annotations.pop('gene', None)
                        if('truncated' in feat):
                            gene_annotations[bc.INSDC_FEATURE_PSEUDO] = True
                        gene_annotations = encode_annotations(gene_annotations)
                        fh.write(f"{seq_id}\tInfernal\tgene\t{start}\t{stop}\t.\t{feat['strand']}\t.\t{gene_annotations}\n")
                    annotations = encode_annotations(annotations)
                    fh.write(f"{seq_id}\tInfernal\t{so.SO_NCRNA_GENE.name}\t{start}\t{stop}\t{feat['evalue']}\t{feat['strand']}\t.\t{annotations}\n")
                elif(feat['type'] == bc.FEATURE_NC_RNA_REGION):
                    annotations = {
                        'ID': feat['id'],
                        'Name': feat['product'],
                        'product': feat['product'],
                        'Dbxref': feat['db_xrefs']
                    }
                    if('truncated' in feat):
                        annotations[bc.INSDC_FEATURE_PSEUDO] = True
                    if(cfg.compliant):
                        for rfam_id in [dbxref.split(':')[1] for dbxref in feat['db_xrefs'] if dbxref.split(':')[0] == bc.DB_XREF_RFAM]:
                            annotations['inference'] = f'profile:Rfam:{rfam_id}'
                        annotations['Dbxref'], annotations['Note'] = insdc.revise_dbxref_insdc(feat['db_xrefs'])  # remove INSDC invalid DbXrefs
                        annotations[bc.INSDC_FEATURE_REGULATORY_CLASS] = insdc.select_regulatory_class(feat)
                    annotations = encode_annotations(annotations)
                    fh.write(f"{seq_id}\tInfernal\t{so.SO_REGULATORY_REGION.name}\t{start}\t{stop}\t{feat['evalue']}\t{feat['strand']}\t.\t{annotations}\n")
                elif(feat['type'] == bc.FEATURE_CRISPR):
                    crispr_tool = "PILER-CR"
                    if prokka:
                        crispr_tool = "MinCED"
                    annotations = {
                        'ID': feat['id'],
                        'Name': feat['product'],
                        'product': feat['product']
                    }
                    feat_type = so.SO_CRISPR.name
                    if(cfg.compliant):
                        feat_type = bc.INSDC_FEATURE_REPEAT_REGION
                        annotations['inference'] = 'COORDINATES:alignment:pilercr:1.02'
                        annotations['Dbxref'], annotations['Note'] = insdc.revise_dbxref_insdc(feat['db_xrefs'])  # remove INSDC invalid DbXrefs
                        annotations[bc.INSDC_FEATURE_REPEAT_FAMILY] = 'CRISPR'
                        annotations[bc.INSDC_FEATURE_REPEAT_TYPE] = 'direct'
                        annotations[bc.INSDC_FEATURE_REPEAT_UNIT_SEQ] = feat['repeat_consensus']
                    annotations = encode_annotations(annotations)
                    fh.write(f"{seq_id}\t{crispr_tool}\t{feat_type}\t{start}\t{stop}\t.\t{feat['strand']}\t.\t{annotations}\n")
                    if(not cfg.compliant):
                        i = 0
                        # spacers and repeats wont exist if Prokka input
                        spacers = feat.get('spacers', [])
                        repeat = feat.get('repeat', [])
                        if len(spacers) > 0 and len(repeat) > 0: 
                            while i < len(feat['spacers']):
                                repeat = feat['repeats'][i]
                                annotations = {
                                    'ID': f"{feat['id']}_repeat_{i+1}",
                                    'Parent': feat['id']
                                }
                                annotations = encode_annotations(annotations)
                                # will always be PILER here as prokka won't have any
                                fh.write(f"{seq_id}\tPILER-CR\t{bc.FEATURE_CRISPR_REPEAT}\t{repeat['start']}\t{repeat['stop']}\t.\t{repeat['strand']}\t.\t{annotations}\n")
                                spacer = feat['spacers'][i]
                                annotations = {
                                    'ID': f"{feat['id']}_spacer_{i+1}",
                                    'Parent': feat['id'],
                                    'sequence': spacer['sequence']
                                }
                                annotations = encode_annotations(annotations)
                                fh.write(f"{seq_id}\tPILER-CR\t{bc.FEATURE_CRISPR_SPACER}\t{spacer['start']}\t{spacer['stop']}\t.\t{spacer['strand']}\t.\t{annotations}\n")
                                i += 1
                            if(len(feat['repeats']) - 1 == i):
                                repeat = feat['repeats'][i]
                                annotations = { 'ID': f"{feat['id']}_repeat_{i+1}" }
                                annotations = encode_annotations(annotations)
                                fh.write(f"{seq_id}\tPILER-CR\t{bc.FEATURE_CRISPR_REPEAT}\t{repeat['start']}\t{repeat['stop']}\t.\t{repeat['strand']}\t.\t{annotations}\n")
                elif feat['type'] == bc.FEATURE_CDS:
                    if euk:
                        write_euk_cds_feature(fh, seq_id, feat)
                    else:
                        annotations = {
                            'ID': feat['locus'],
                            'Name': feat['product'],
                            'locus_tag': feat['locus'],
                            'product': feat['product'],
                            'Dbxref': feat['db_xrefs']
                        }
                        if(bc.PSEUDOGENE in feat):
                            annotations[bc.INSDC_FEATURE_PSEUDOGENE] = bc.INSDC_FEATURE_PSEUDOGENE_TYPE_UNPROCESSED if feat[bc.PSEUDOGENE]['paralog'] else bc.INSDC_FEATURE_PSEUDOGENE_TYPE_UNITARY
                        elif('truncated' in feat):
                            annotations[bc.INSDC_FEATURE_PSEUDO] = True
                        if(feat.get('gene', None)):  # add gene annotation if available
                            annotations['gene'] = feat['gene']
                        source = '?' if feat.get('source', None) == bc.CDS_SOURCE_USER else 'Pyrodigal'
                        if prokka: 
                            source = 'Prodigal'
                        if(cfg.compliant):
                            gene_id = f"{feat['locus']}_gene"
                            annotations['Parent'] = gene_id
                            annotations['inference'] = 'EXISTENCE:non-experimental evidence, no additional details recorded' if feat.get('source', None) == bc.CDS_SOURCE_USER else 'ab initio prediction:Pyrodigal:3.5'
                            annotations['Dbxref'], annotations['Note'] = insdc.revise_dbxref_insdc(feat['db_xrefs'])  # remove INSDC invalid DbXrefs
                            annotations['Note'], ec_number = insdc.extract_ec_from_notes_insdc(annotations, 'Note')
                            if(ec_number is not None):
                                annotations['ec_number'] = ec_number
                            gene_annotations = {
                                'ID': gene_id,
                                'locus_tag': feat['locus']
                            }
                            if(feat.get('gene', None)):
                                gene_annotations['gene'] = feat['gene']
                            if(bc.PSEUDOGENE in feat):
                                gene_annotations[bc.INSDC_FEATURE_PSEUDOGENE] = bc.INSDC_FEATURE_PSEUDOGENE_TYPE_UNPROCESSED if feat[bc.PSEUDOGENE]['paralog'] else bc.INSDC_FEATURE_PSEUDOGENE_TYPE_UNITARY
                            gene_annotations = encode_annotations(gene_annotations)
                            fh.write(f"{seq_id}\t{source}\tgene\t{start}\t{stop}\t.\t{feat['strand']}\t.\t{gene_annotations}\n")
                        if('exception' in feat):
                            ex = feat['exception']
                            pos = f"{ex['start']}..{ex['stop']}"
                            if(feat['strand'] == bc.STRAND_REVERSE):
                                pos = f"complement({pos})"
                            annotations['transl_except']=f"(pos:{pos},aa:{ex['aa']})"
                            notes = annotations.get('Note', [])
                            notes.append(f"codon on position {ex['codon_position']} is a {ex['type']} codon")
                            if('Notes' not in annotations):
                                annotations['Note'] = notes
                        annotations = encode_annotations(annotations)
                        fh.write(f"{seq_id}\t{source}\t{so.SO_CDS.name}\t{start}\t{stop}\t.\t{feat['strand']}\t0\t{annotations}\n")
                        if(bc.FEATURE_SIGNAL_PEPTIDE in feat):
                            write_signal_peptide(fh, feat)
                elif(feat['type'] == bc.FEATURE_SORF):
                    annotations = {
                        'ID': feat['locus'],
                        'Name': feat['product'],
                        'locus_tag': feat['locus'],
                        'product': feat['product'],
                        'Dbxref': feat['db_xrefs']
                    }
                    if(feat.get('gene', None)):  # add gene annotation if available
                        annotations['gene'] = feat['gene']
                    if(cfg.compliant):
                        gene_id = f"{feat['locus']}_gene"
                        annotations['Parent'] = gene_id
                        annotations['Dbxref'], annotations['Note'] = insdc.revise_dbxref_insdc(feat['db_xrefs'])  # remove INSDC invalid DbXrefs
                        annotations['Note'], ec_number = insdc.extract_ec_from_notes_insdc(annotations, 'Note')
                        if(ec_number is not None):
                            annotations['ec_number'] = ec_number
                        gene_annotations = {
                            'ID': gene_id,
                            'locus_tag': feat['locus'],
                            'inference': 'ab initio prediction:Bakta'
                        }
                        if(feat.get('gene', None)):
                            gene_annotations['gene'] = feat['gene']
                        gene_annotations = encode_annotations(gene_annotations)
                        fh.write(f"{seq_id}\tBakta\tgene\t{start}\t{stop}\t.\t{feat['strand']}\t.\t{gene_annotations}\n")
                    annotations = encode_annotations(annotations)
                    fh.write(f"{seq_id}\tBakta\t{so.SO_CDS.name}\t{start}\t{stop}\t.\t{feat['strand']}\t0\t{annotations}\n")
                    if(bc.FEATURE_SIGNAL_PEPTIDE in feat):
                        write_signal_peptide(fh, feat)
                elif(feat['type'] == bc.FEATURE_GAP):
                    gap_tool="Bakta"
                    if prokka:
                        gap_tool="Prokka"
                    annotations = {
                        'ID': feat['id'],
                        'Name': f"gap ({feat['length']} bp)",
                        'product': f"gap ({feat['length']} bp)"
                    }
                    annotations = encode_annotations(annotations)
                    fh.write(f"{seq_id}\t{gap_tool}\t{so.SO_GAP.name}\t{start}\t{stop}\t.\t{feat['strand']}\t.\t{annotations}\n")
                elif(feat['type'] == bc.FEATURE_ORIC):
                    annotations = {
                        'ID': feat['id'],
                        'Name': feat['product']
                    }
                    if(cfg.compliant):
                        annotations['Note'] = feat['product']
                    else:
                        annotations['product'] = feat['product']
                        annotations['inference'] = 'similar to DNA sequence'
                    annotations = encode_annotations(annotations)
                    feat_type = bc.INSDC_FEATURE_ORIGIN_REPLICATION if cfg.compliant else so.SO_ORIC.name
                    fh.write(f"{seq_id}\tBLAST+\t{feat_type}\t{start}\t{stop}\t.\t{feat['strand']}\t.\t{annotations}\n")
                elif(feat['type'] == bc.FEATURE_ORIV):
                    annotations = {
                        'ID': feat['id'],
                        'Name': feat['product']
                    }
                    if(cfg.compliant):
                        annotations['Note'] = feat['product']
                    else:
                        annotations['product'] = feat['product']
                        annotations['inference'] = 'similar to DNA sequence'
                    annotations = encode_annotations(annotations)
                    feat_type = bc.INSDC_FEATURE_ORIGIN_REPLICATION if cfg.compliant else so.SO_ORIC.name
                    fh.write(f"{seq_id}\tBLAST+\t{feat_type}\t{start}\t{stop}\t.\t{feat['strand']}\t.\t{annotations}\n")
                elif(feat['type'] == bc.FEATURE_ORIT):
                    annotations = {
                        'ID': feat['id'],
                        'Name': feat['product']
                    }
                    if(cfg.compliant):
                        annotations['Note'] = feat['product']
                    else:
                        annotations['product'] = feat['product']
                        annotations['inference'] = 'similar to DNA sequence'
                    annotations = encode_annotations(annotations)
                    feat_type = bc.INSDC_FEATURE_ORIGIN_TRANSFER if cfg.compliant else so.SO_ORIT.name
                    fh.write(f"{seq_id}\tBLAST+\t{feat_type}\t{start}\t{stop}\t.\t{feat['strand']}\t.\t{annotations}\n")
                elif(feat['type'] == bc.FEATURE_GENE):
                    write_gene_feature(fh, seq_id, feat)
                elif(feat['type'] == bc.FEATURE_MRNA):
                    write_mrna_feature(fh, seq_id, feat)

        if(not cfg.compliant):
            fh.write('##FASTA\n')
            for seq in data['sequences']:  # write sequences
                fh.write(f">{seq['id']}\n")
                seq_nt = seq['nt'] if 'nt' in seq else seq['sequence']  # <1.10.0 compatibility
                fh.write(fasta.wrap_sequence(seq_nt))
    return


def encode_attribute(product: str) -> str:
    """Replace special characters forbidden in column 9 of the GFF3 format: https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md"""
    product = str(product)
    product = product.replace('%', '%25')
    product = product.replace(';', '%3B')
    product = product.replace('=', '%3D')
    product = product.replace('&', '%26')
    product = product.replace(',', '%2C')
    return product


def encode_annotations(annotations: Dict[str, Union[str, Sequence[str]]]) -> str:
    """
    Encodes annotations into a string.

    Args:
      annotations (dict): A dictionary containing the annotations.

    Returns:
      str: The encoded annotations.

    Examples:
      >>> encode_annotations({
          'ID': 'EHICP_3230_sigpep',
          'Name': 'signal peptide',
          'product': 'signal peptide',
          'score': 0.5,
          'Parent': 'EHICP_3230'
      })
      'ID=EHICP_3230_sigpep;Name=signal peptide;product=signal peptide;score=0.5;Parent=EHICP_3230'
    """
    annotation_strings = []
    for key, val in annotations.items():
        if(type(val) is list):
            if(len(val) >= 1):
                val = [encode_attribute(k) for k in val]
                annotation = f"{key}={','.join(val)}"
                annotation_strings.append(annotation)
        else:
            annotation_strings.append(f'{key}={encode_attribute(val)}')
    return ';'.join(annotation_strings)


def write_signal_peptide(fh, feat: dict):  # <1.10.0 compatibility
    """
    Writes a signal peptide feature to a file.

    Args:
      fh (file): The file handle to write to.
      feat (dict): A dictionary containing the feature information.

    Returns:
      None

    Examples:
      >>> write_signal_peptide(fh, {
          'locus': 'EHICP_3230',
          'sequence': 'DS571531.1',
          'strand': '+',
          'signal_peptide': {
              'start': 1,
              'stop': 20,
              'score': 0.5
          }
      })
    """
    sig_peptide = feat[bc.FEATURE_SIGNAL_PEPTIDE]
    annotations = {
        'ID': f"{feat['locus']}_sigpep",
        'Name': 'signal peptide',
        'product': 'signal peptide',
        'score': sig_peptide['score'],
        'Parent': feat['locus']
    }
    annotations = encode_annotations(annotations)
    seq_id = feat['sequence'] if 'sequence' in feat else feat['contig']  # <1.10.0 compatibility
    fh.write(f"{seq_id}\tDeepSig\t{so.SO_SIGNAL_PEPTIDE.name}\t{sig_peptide['start']}\t{sig_peptide['stop']}\t{sig_peptide['score']:.2f}\t{feat['strand']}\t.\t{annotations}\n")
