#!/usr/bin/env python3

import hashlib
import json
import random
import string
# from Bio import SeqIO
# import argparse
# import re
from collections import defaultdict
from datetime import datetime

from loguru import logger
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.SeqUtils import gc_fraction

import baktfold.bakta.so as so


AMINO_ACID_DICT = {
    'ala': ('A', so.SO_TRNA_ALA),
    'gln': ('Q', so.SO_TRNA_GLN),
    'glu': ('E', so.SO_TRNA_GLU),
    'gly': ('G', so.SO_TRNA_GLY),
    'pro': ('P', so.SO_TRNA_PRO),
    'met': ('M', so.SO_TRNA_MET),
    'fmet':('fM', so.SO_TRNA_MET),
    'asp': ('D', so.SO_TRNA_ASP),
    'thr': ('T', so.SO_TRNA_THR),
    'val': ('V', so.SO_TRNA_VAL),
    'tyr': ('Y', so.SO_TRNA_TYR),
    'cys': ('C', so.SO_TRNA_CYS),
    'ile': ('I', so.SO_TRNA_ILE),
    'ile2':('I', so.SO_TRNA_ILE),
    'ser': ('S', so.SO_TRNA_SER),
    'leu': ('L', so.SO_TRNA_LEU),
    'trp': ('W', so.SO_TRNA_TRP),
    'lys': ('K', so.SO_TRNA_LYS),
    'asn': ('N', so.SO_TRNA_ASN),
    'arg': ('R', so.SO_TRNA_ARG),
    'his': ('H', so.SO_TRNA_HIS),
    'phe': ('F', so.SO_TRNA_PHE),
    'sec': ('U', so.SO_TRNA_SELCYS)
}


def convert_cds_feature(feature, seq_record, translation_table, id):
    """
    Convert a Prokka CDS Biopython SeqFeature to a Bakta CDS JSON entry.
    """

    # ----------- Location info -----------

    strand = "+" if feature.location.strand == 1 else "-"

    if strand == -1:  # negative strand
        start = int(feature.location.end)     
        stop  = int(feature.location.start) - 1  
    else:  # positive strand
        start = int(feature.location.start) + 1  
        stop  = int(feature.location.end)        

    if feature.location.__class__.__name__ == "CompoundLocation":
        # Multi-exon (join)
        starts = []
        stops = []
        for part in feature.location.parts:
            if strand == -1:
                # For minus strand, 5' is end, 3' is start
                starts.append(int(part.end))
                stops.append(int(part.start) - 1)
            else:
                starts.append(int(part.start) + 1)
                stops.append(int(part.end))

    else:
        starts = None
        stops = None


    # frame: Bakta uses 1/2/3; Prokka codon_start is ["1","2","3"]
    codon_start = int(feature.qualifiers.get("codon_start", ["1"])[0])
    frame = codon_start

    # ----------- Basic qualifiers -----------
    gene = feature.qualifiers.get("gene", [None])[0]
    product = feature.qualifiers.get("product", [None])[0]

    locus_tag = feature.qualifiers.get("locus_tag", [None])[0]
    note = feature.qualifiers.get("note", [None])[0]
    locus = locus_tag

    # pseudo
    pseudo = feature.qualifiers.get("pseudo", [None])[0]

    protein_id = feature.qualifiers.get("protein_id", [None])[0]

    # ----------- Extract nucleotides -----------
    nt_seq = feature.extract(seq_record.seq)
    nt = str(nt_seq)

    # ----------- Extract amino acids -----------
    aa = feature.qualifiers.get("translation", [""])[0]

    # Compute translation if Prokka didn't provide it
    if not aa:
        try:
            aa = str(nt_seq.translate(table=translation_table, cds=True))
        except Exception:
            aa = ""

    # ----------- aa MD5 hexdigest -----------
    aa_hexdigest = hashlib.md5(aa.encode()).hexdigest()

    # ----------- Hypothetical? -----------
    hypothetical = product is None or "hypothetical protein" in product.lower()

    # ----------- Compute protein stats -----------
    seq_stats = None
    if aa:
        try:
            analysed = ProteinAnalysis(aa)
            seq_stats = {
                "molecular_weight": analysed.molecular_weight(),
                "isoelectric_point": analysed.isoelectric_point()
            }
        except Exception:
            seq_stats = None

    # Get existing db_xref list or default to [so.SO_CDS.id]
    db_xref = feature.qualifiers.get("db_xref", [so.SO_CDS.id])

    # Append so.SO_CDS.id only if it’s not already present
    if so.SO_CDS.id not in db_xref:
        db_xref.append(so.SO_CDS.id)

    # ----------- Make Bakta-format dict -----------
    bakta_cds = {
        "type": "cds",
        "sequence": seq_record.id,
        "start": start,
        "stop": stop,
        "starts": starts,
        "stops": stops,
        "strand": strand,
        "frame": frame,
        "gene": gene,
        "product": product,
        "db_xrefs": db_xref,  
        "nt": nt,
        "aa": aa,
        "aa_hexdigest": aa_hexdigest,
        "start_type": None,
        "rbs_motif": None,
        "genes": [],
        "note": note,
        "seq_stats": seq_stats,
        "id": id,
        "locus": locus,
        "protein_id": protein_id
    }

    if pseudo:
        bakta_cds["pseudo"] = True

    if hypothetical:
        bakta_cds["hypothetical"] = True

    return bakta_cds

def convert_trna_feature(feature, seq_record, id):
    """
    Convert a funannotate tRNA SeqFeature to a Bakta tRNA JSON entry.
    """

    # ------------ Location ------------

    strand = "+" if feature.location.strand == 1 else "-"


    if strand == -1:  # negative strand
        start = int(feature.location.end)     
        stop  = int(feature.location.start) - 1  
    else:  # positive strand
        start = int(feature.location.start) + 1  
        stop  = int(feature.location.end)        



    # ------------ Extract nt sequence ------------
    nt_seq = feature.extract(seq_record.seq)
    nt = str(nt_seq)

    # ------------ Basic qualifiers ------------
    product = feature.qualifiers.get("product", [None])[0]
    locus = feature.qualifiers.get("locus_tag", [None])[0]

    # ------------ amino acid ------------
    # Prokka product examples:
    #   "tRNA-Trp"
    #   "tRNA-Leu"
    amino_acid = None
    if product and product.startswith("tRNA-"):
        amino_acid = product.split("-")[1]


    # ------------ anticodon ------------
    anti_codon = None
    
    # anticodons are in notes

    notes = feature.qualifiers.get("note", [])

    # Expect a note like: "tRNA-Ser(gga)"
    for note in notes:
        # Remove spaces for safety
        n = note.replace(" ", "")

        # Extract part inside parentheses (anticodon)
        if "(" in n and ")" in n:
            anti_codon = n.split("(")[1].split(")")[0].lower()

        # Extract amino acid:
        # tRNA-Ser(gga) → "Ser"
        if "tRNA-" in n:
            try:
                # tRNA-Ser(gga) → "Ser(gga)" → split('(')[0] → "Ser"
                aa_section = n.split("tRNA-")[1]
                aa_clean = aa_section.split("(")[0]
                amino_acid = aa_clean
            except Exception:
                pass

    # ------------ Anti-codon position detection ------------
    # Prokka doesnt have it - dont include
    # anti_codon_pos = None

    # ------------ score ------------
    # nothing in prokka
    score = None

    # ------------ db_xrefs ------------
    # doesnt exist for prokka
    db_xrefs = feature.qualifiers.get("db_xref", [])
    # add so_term
    so_term = AMINO_ACID_DICT.get(amino_acid.lower(), ('', None))[1]

    if (so_term):
        db_xrefs.append(so_term.id)

    # ------------ final Bakta-form dict ------------
    bakta_trna = {
        "type": "tRNA",
        "sequence": seq_record.id,
        "start": start,
        "stop": stop,
        "strand": strand,
        "gene": "trn" + (amino_acid[0].lower() if amino_acid else "?"),
        "product": product,
        "amino_acid": amino_acid,
        "anti_codon": anti_codon,
        "score": score,
        "nt": nt,
        "db_xrefs": db_xrefs,
       #  "anti_codon_pos": anti_codon_pos,  dont include, not in output
        "id": id,
        "locus": locus
    }

    return bakta_trna

def convert_gene_feature(feature, rec, id):
    """
    Convert a Funannotate GenBank gene feature to Bakta-style JSON.
    
    Parameters:
        feature: Bio.SeqFeature
            The rRNA feature from the GBK.
        rec: str
            The record from the GBK.
    Returns:
        dict: Bakta-style rRNA feature
    """

    strand = "+" if feature.location.strand == 1 else "-"
    

    if strand == -1:  # negative strand
        start = int(feature.location.end)     
        stop  = int(feature.location.start) - 1  
    else:  # positive strand
        start = int(feature.location.start) + 1  
        stop  = int(feature.location.end)        




    qualifiers = feature.qualifiers

    locus_tag = qualifiers.get("locus_tag", [None])[0]
    gene = qualifiers.get("gene", [None])[0]

    pseudo = qualifiers.get("pseudo", [None])[0]
    

    gene_entry = {
        "type": "gene",
        "sequence": rec.id,
        "start": start,
        "stop": stop,
        "strand": strand,
        "gene": gene,
        "db_xrefs": [so.SO_GENE.id], 
        "id": id,
        "locus": locus_tag
    }

    if pseudo:
        gene_entry["pseudo"] = True
    
    return gene_entry

def convert_mrna_feature(feature, rec, id):
    """
    Convert a funannotate mrna feature to a Bakta-style feature.

    Parameters:
        feature: Bio.SeqFeature
            The mRNA feature from the GBK.
        rec: Bio.SeqRecord
            The record containing the sequence.

    Returns:
        dict: Bakta-style misc_RNA feature.
    """

    seq = str(rec.seq)

    # Coordinates (GBK is 0-based, Bakta is 1-based)
    start = int(feature.location.start) + 1
    stop = int(feature.location.end)


    # Extract location
    strand = "+" if feature.location.strand == 1 else "-"

    if strand == -1:  # negative strand
        start = int(feature.location.end)     
        stop  = int(feature.location.start) - 1  
    else:  # positive strand
        start = int(feature.location.start) + 1  
        stop  = int(feature.location.end)        

    if feature.location.__class__.__name__ == "CompoundLocation":
        # Multi-exon (join)
        starts = []
        stops = []
        for part in feature.location.parts:
            if strand == -1:
                # For minus strand, 5' is end, 3' is start
                starts.append(int(part.end))
                stops.append(int(part.start) - 1)
            else:
                starts.append(int(part.start) + 1)
                stops.append(int(part.end))

    else:
        starts = None
        stops = None



    qualifiers = feature.qualifiers

    gene = qualifiers.get("gene", [None])[0]
    product = qualifiers.get("product", [None])[0]
    locus_tag = qualifiers.get("locus_tag", [None])[0]
    # only on some
    standard_name = qualifiers.get("standard_name", [None])[0]
    pseudo = qualifiers.get("pseudo", [None])[0]


    mrna_entry = {
        "type": "mRNA", 
        "sequence": rec.id,
        "start": start,
        "stop": stop,
        "starts": starts,
        "stops": stops,
        "strand": strand,
        "gene": gene,
        "product": product,
        "db_xrefs": [so.SO_MRNA.id],        
        "id": id,
        "locus": locus_tag
    }

    if standard_name:
        mrna_entry["standard_name"] = standard_name

    if pseudo:
        mrna_entry["pseudo"] = True

    return mrna_entry


def convert_assembly_gap_feature(feature, rec, id):
    """
    Convert a GenBank assembly_gap feature to a simplified Bakta-style 'gap' feature.

    Parameters:
        feature: Bio.SeqFeature
            The assembly_gap feature from the GBK.
        rec: Bio.SeqRecord
            The full GenBank record containing the sequence.

    Returns:
        dict: Simplified Bakta-style gap feature.
    """

    # Coordinates (1-based)
    start = int(feature.location.start) + 1
    stop = int(feature.location.end)

    qualifiers = feature.qualifiers

    #  may provide estimated_length but coordinates already give an exact span
    est_len = qualifiers.get("estimated_length", [None])[0]
    if est_len is not None:
        length = int(est_len)
    else:
        length = stop - start + 1  # fallback from coordinates

    # Bakta always uses "." for strand on gaps
    strand = "."

    gap_entry = {
        "type": "gap",
        "sequence": rec.id,
        "start": start,
        "stop": stop,
        "strand": strand,
        "length": length,
        "id": id
    }

    # no need to add estimated length separately - it is covered by length in the json 

    # if est_len:
    #     gap_entry["estimated_length"] = est_len

    return gap_entry


def convert_repeat_region_feature(feature, rec, id):
    """
    Convert a Prokka GenBank repeat_region (CRISPR) feature to a simplified Bakta-style feature.

    Parameters:
        feature: Bio.SeqFeature
            The repeat_region feature (crispr) from the Prokka GBK.
        rec: Bio.SeqRecord
            The full GenBank record containing the sequence.

    Returns:
        dict: Simplified Bakta-style CRISPR feature.
    """

    # Coordinates (Bakta uses 1-based)
    start = int(feature.location.start) + 1
    stop = int(feature.location.end)

    qualifiers = feature.qualifiers
    note = qualifiers.get("note", [None])[0]
    rpt_family = qualifiers.get("rpt_family", [None])[0]
    rpt_type = qualifiers.get("rpt_type", [None])[0]
    rpt_unit_seq = qualifiers.get("rpt_unit_seq", [None])[0]

    strand = "?"

    # always just take the positive strand to get the NT seq (crispr repeat region)
    seq =  str(rec.seq)
    nt_seq = seq[start-1:stop]

    # Minimal Bakta-like CRISPR structure
    repeat_region_entry = {
        "type": "repeat_region",
        "sequence": rec.id,
        "start": start,
        "stop": stop,
        "strand": strand, # matches Bakta and is required
        "family": rpt_family,       # e.g., "LINE1" - should always be there
        "rpt_type": rpt_type,   
        "repeat_unit": rpt_unit_seq, # the actual consensus repeat if crispr
        "product": note, # won't be the same as Bakta as different lookup method used - but needed for the gff writing
        "nt": nt_seq, # needed for batka .ffn writeout
        "id": id, # bakta_id needed 
        # "locus": None, # no locus tag like Bakta
        "db_xrefs": [so.SO_REPEAT.id]
    }

    return repeat_region_entry

def convert_utr_region_feature(feature, rec, id, three):
    """
    Convert a UTR GenBank feature to a simplified Bakta-style feature.

    Parameters:
        feature: Bio.SeqFeature
            The UTR feature from the GBK.
        rec: Bio.SeqRecord
            The full GenBank record containing the sequence.

    Returns:
        dict: Simplified Bakta-style feature.
    """

    if three:
        type = "3'UTR"
        so_code =  so.SO_3UTR.id
    else:
        type = "5'UTR"
        so_code =  so.SO_5UTR.id

    # get position - can be neg or pos strand
    strand = "+" if feature.location.strand == 1 else "-"
    
    if strand == -1:  # negative strand
        start = int(feature.location.end)     
        stop  = int(feature.location.start) - 1  
    else:  # positive strand
        start = int(feature.location.start) + 1  
        stop  = int(feature.location.end)     

    qualifiers = feature.qualifiers
    note = qualifiers.get("note", [None])[0]

    locus_tag = qualifiers.get("locus_tag", [None])[0]

    # always just take the positive strand to get the NT seq (UTR region)
    seq =  str(rec.seq)
    nt_seq = seq[start-1:stop]


    # Minimal Bakta-like CRISPR structure
    utr_entry = {
        "type": type,
        "sequence": rec.id,
        "start": start,
        "stop": stop,
        "strand": strand, # matches Bakta and is required
        "product": note, 
        "nt": nt_seq, # needed for batka .ffn writeout
        "id": id, # bakta_id needed 
        "locus": locus_tag,
        "db_xrefs": [so_code]
    }

    return utr_entry

    #  3'UTR           186716..186779
    #                  /locus_tag="C1SCF055_LOCUS2256"
    #                  /note="ID:SCF055_s1087_g24103.utr3p1;
    #                  source:feature"

    #  5'UTR           218227..218248
    #                  /locus_tag="C1SCF055_LOCUS2213"
    #                  /note="ID:SCF055_s1083_g24060.utr5p1;
    #                  source:feature"

def build_bakta_sequence_entry(rec):
    """
    Convert a  SeqRecord into a Bakta-style sequence entry.
    Missing fields are filled with None.
    """

    seq = str(rec.seq)

    # -----------------------------------------
    # Extract source feature qualifiers - genbank always has source field
    # -----------------------------------------
    source_feat = next((f for f in rec.features if f.type == "source"), None)

    source_qualifiers = {}

    # Defaults (None) for all fields
    organism = None
    strain = None
    db_xref = None
    note = None

    if source_feat:
        q = source_feat.qualifiers

        # Always include mol_type if it exists
        if "mol_type" in q:
            mol_type = q.get("mol_type", [None])[0]

        # Optional qualifiers
        if "organism" in q:
            organism = q["organism"][0]

        if "strain" in q:
            strain = q["strain"][0]

        if "db_xref" in q:
            # db_xref is usually a 1-element list
            val = q["db_xref"]
            db_xref = val[0] if len(val) == 1 else val

        if "note" in q:
            note = q["note"][0]
            
    entry = {
        "id": rec.id,                                 #  contig name
        "description": None,                          #  does not have topology completeness etc 
        "nt": seq,                                    # Nucleotide sequence
        "length": len(seq),                           
        "complete": None,                             #  does not define completeness
        "type": None,                                 # plasmid/chromosome not known
        "topology": rec.annotations["topology"],                             # circular/linear unknown
        "simple_id": rec.id,                          # Same as id (placeholder)
        "orig_id": rec.id,                            # No separate original ID fr
        "orig_description": None,                     # Same as description
        "name": None,                          
    }

    # Add source qualifiers only if they exist
    if organism is not None:
        entry["organism"] = organism

    if mol_type is not None:
        entry["mol_type"] = mol_type

    if strain is not None:
        entry["strain"] = strain

    if db_xref is not None:
        entry["db_xref"] = db_xref

    if note is not None:
        entry["note"] = note

    return entry


# def parse_funannotate_version(record):
#     """
#     Extract Funannotate version from COMMENT field.
#     Example COMMENT: 'Annotated using 1.8.9'.
#     """

#     comments = record.annotations.get("comment", "") or record.annotations.get("comments", "")
#     if not comments:
#         return "unknown"

#     # Match version number like 1.8.9
#     m = re.search(r"\b(\d+\.\d+\.\d+)\b", comments)
#     if m:
#         return m.group(1)

#     return "unknown"

def random_n_letter_id(n=4):
    """
    generates a n letter id prefix 
    
    n=2 to append to   locus tag  for bakta id to make it different
    n=10 if the locus tag is somehow missing (should never happen) 
    """
    return ''.join(random.choices(string.ascii_uppercase, k=n))

def get_bakta_style_id_from_locus_tag(records):
    """
    Gets 10 char bakta-style ID tag based off the 8 char locus tag in first CDS on the first  record + 2 random chars

    Assumes all records will have the same locus tag prefix

    Will always add 2 chars to make ID unique vs locus tag
    """

    if not records:
        raise ValueError("No GenBank records found.")
    
    for record in records:
    
        for feat in record.features:
            if feat.type == "CDS":
                locus_tag_list = feat.qualifiers.get("locus_tag") # returns None if doesn't exist
                
                if locus_tag_list:
                    locus_tag = locus_tag_list[0]

                    if len(locus_tag) > 7:

                        locus_tag_prefix = locus_tag[:-7] # trims off _000001 from CDS

                        rand_two_chars = random_n_letter_id(2)

                        # by default  locus tag is 8 chars. So this returns a 10 char string (same as bakta defaults)

                        id_tag = f"{locus_tag_prefix}{rand_two_chars}"

                        return id_tag
                  

                    else:
                        return random_n_letter_id(10)

                # fallback if locus_tag missing or too short
                return random_n_letter_id(10)

    # No CDS feature found at all (shouldn't happen)
    return random_n_letter_id(10)



def calc_genome_stats(records):
    """
    Compute correct genome stats (size, GC, N-ratio, N50, N90) for records from a multi-contig
     GenBank file.
    """

    if not records:
        raise ValueError("No GenBank records found.")

    # lengths of all contigs
    contig_lengths = [len(r.seq) for r in records]
    total_length = sum(contig_lengths)

    # concatenate sequences for global GC + N calculation
    full_seq = "".join(str(r.seq) for r in records)

    # GC as fraction (Bakta wants 0–1)
    gc_perc = gc_fraction(full_seq)

    # N-ratio
    n_ratio = full_seq.count("N") / total_length

    # ---------- N50 / N90 ----------
    sorted_lengths = sorted(contig_lengths, reverse=True)

    def nx_metric(sorted_lens, total, threshold):
        """
        Generic N{threshold} function.
        threshold: 0.5 for N50, 0.9 for N90
        """
        cutoff = total * threshold
        running = 0
        for l in sorted_lens:
            running += l
            if running >= cutoff:
                return l
        return sorted_lens[-1]  # fallback (should not happen)

    n50 = nx_metric(sorted_lengths, total_length, 0.5)
    n90 = nx_metric(sorted_lengths, total_length, 0.9)

    return {
        "size": total_length,
        "gc": gc_perc,
        "n_ratio": n_ratio,
        "n50": n50,
        "n90": n90,
        "coding_ratio": None  
    }

def eukaryotic_gbk_to_json(records, output_json):
    
    if len(records) == 0:
        raise ValueError("No GenBank records found.")
    
    # elif len(records) >= 1:
    #     funannotate_version = parse_funannotate_version(records[0])

    # unique record types
    unique_feature_types = {
        feature.type
        for record in records
        for feature in record.features
    }

    ORDER = ["tRNA", "gene", "mRNA", "CDS", "assembly_gap", "gap", "repeat_region", "5'UTR", "3'UTR"]

     # source always in input - it is made in output anyway
    covered_set = set(ORDER + ["source"])

    # Compute features in records not in the covered list
    uncovered_features = unique_feature_types - covered_set

    if uncovered_features:
        logger.warning("Feature types present in your input GenBank but not convertible include:")
        for ft in sorted(uncovered_features):
            logger.warning(ft)
        logger.warning("Baktfold will always only work on CDS.")
        logger.warning("If you would like to have this feature type supported in GFF3/GenBank output, please make an issue at https://github.com/gbouras13/baktfold")
    else:
        logger.info("All feature types in your input GenBank are convertible.")


    # translation_table = get_transl_table(records)
    # not in funannotate - generally just assume blank
    translation_table = "Unknown"

    bakta_id_prefix = get_bakta_style_id_from_locus_tag(records)
    
    # ----------------------------
    # Genome block 
    # ----------------------------

    genome_block = {
        "genus": None,
        "species": None,
        "strain": None,
        "taxon": None,
        "complete": True,
        "gram": "?",
        "translation_table": translation_table
    }

    # ----------------------------
    # Stats - on whole GBK
    # ----------------------------
    stats_block = calc_genome_stats(records)

    # ----------------------------
    # Features block
    # ----------------------------

    
    # as far as I can tell these are all the relevant ones for funannotate
    # added 5'UTR and 3'UTR from genbank entry CAMXCT020000001 - assume it is common
    # can't add everything so for now this'll do

    features = []
    i = 1

    for rec in records:
        
        for ftype in ORDER:
            for feat in rec.features:
                if feat.type != ftype:
                    continue

                id = f"{bakta_id_prefix}_{i}"

                if ftype == "CDS":
                    features.append(convert_cds_feature(feat, rec, translation_table, id))
                elif ftype == "mRNA":
                    features.append(convert_mrna_feature(feat, rec, id))
                elif ftype == "gene":
                    features.append(convert_gene_feature(feat, rec, id))
                elif ftype == "tRNA":
                    features.append(convert_trna_feature(feat, rec, id))
                elif ftype == "assembly_gap" or ftype == "gap":
                    features.append(convert_assembly_gap_feature(feat, rec, id))
                elif ftype == "repeat_region":
                    features.append(convert_repeat_region_feature(feat, rec, id))
                elif ftype == "5'UTR":
                    features.append(convert_utr_region_feature(feat, rec, id, three=False))
                elif ftype == "3'UTR":
                    features.append(convert_utr_region_feature(feat, rec, id, three=True))              
                i +=1


    # ----------------------------
    # Sort features within each contig like Bakta
    # ----------------------------

    features_by_contig = defaultdict(list)
    for f in features:
        features_by_contig[f["sequence"]].append(f)

    # Sort each contig's features by start and flatten back
    sorted_features = []
    for contig_id in features_by_contig:
        sorted_features.extend(sorted(features_by_contig[contig_id], key=lambda x: x["start"]))

    # Replace the original features list
    features = sorted_features

    # ----------------------------
    # Sequences block
    # ----------------------------

    sequences = []
    for rec in records:
        sequences.append(build_bakta_sequence_entry(rec))

    # just to put in a time
    start_time = datetime.now()

    bakta_json = {
        "genome": genome_block,
        "stats": stats_block,
        "features": features,
        "sequences": sequences,
        "run": {
            "start": start_time.strftime("%Y-%m-%d %H:%M:%S"),
            "end": start_time.strftime("%Y-%m-%d %H:%M:%S"),       
            "duration": "0.00 min"   
        },
        "version": {
            "eukaryotic": "unknown",  
            "db": {
                "version": "unknown",
                "type": "unknown_dbs"
            }
        }
    }




    with open(output_json, "w") as fh:
        json.dump(bakta_json, fh, indent=4)





