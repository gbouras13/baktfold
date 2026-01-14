#!/usr/bin/env python3

import hashlib
import json
import random
import string
# from Bio import SeqIO
# import argparse
# import re
from collections import defaultdict

from tqdm import tqdm
import hashlib

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

# global variable for the genome-level random string - backup if no locus_tags in input genbank
GENOME_RANDOM_BACKUP_LOCUSTAG_STR = ''.join(random.choices(string.ascii_uppercase + string.digits, k=6))

def add_optional_qualifiers(entry, qualifiers, single_valued=None, multi_valued=None):
    """
    Add optional INSDC qualifiers to a feature entry dict in Bakta style.

    Parameters
    ----------
    entry : dict
        The feature dictionary being built.
    qualifiers : dict
        The qualifiers dictionary from Bio.SeqFeature.
    single_valued : set or list
        Qualifiers expected to be single-valued (take the first if multiple).
    multi_valued : set or list
        Qualifiers that can have multiple values (keep as list if >1, else single value).
    """

    single_valued = single_valued or set()
    multi_valued = multi_valued or set()

    # Multi-valued qualifiers
    for key in multi_valued:
        vals = qualifiers.get(key)
        if vals:
            entry[key] = vals if len(vals) > 1 else vals[0]

    # Single-valued qualifiers
    for key in single_valued:
        vals = qualifiers.get(key)
        if vals:
            if key == "locus_tag":
                entry["locus"] = vals[0] # this is what bakta needs
            else:
                entry[key] = vals[0]

def convert_cds_feature(feature, seq_record, translation_table, id):
    """
    Convert a Prokka CDS Biopython SeqFeature to a Bakta CDS JSON entry.
    """

    # ----------- Location info -----------

    strand = "+" if feature.location.strand == 1 else "-"
    start = int(feature.location.start) + 1   # Bakta uses 1-based inclusive
    stop  = int(feature.location.end)         # already inclusive after conversion

    # Handle CompoundLocation (join)
    starts = None
    stops = None

    if feature.location.__class__.__name__ == "CompoundLocation":
        starts = []
        stops = []
        for part in feature.location.parts:
            starts.append(int(part.start) + 1)
            stops.append(int(part.end))


    # frame: Bakta uses 1/2/3; Prokka codon_start is ["1","2","3"]
    codon_start = int(feature.qualifiers.get("codon_start", ["1"])[0])
    frame = codon_start

    qualifiers = feature.qualifiers

    # ----------- Basic qualifiers -----------
    gene = qualifiers.get("gene", [None])[0]
    product = qualifiers.get("product", [None])[0]


    # fall back to start_stop_strand if there is no locus tag
    if 'locus_tag' in qualifiers and qualifiers['locus_tag']:
        locus_tag = qualifiers['locus_tag'][0]
    else:
        logger.warning(f"No locus_tag found for feature {id}")
        locus_tag = f"{GENOME_RANDOM_BACKUP_LOCUSTAG_STR}_{start}_{stop}"
        logger.warning(f"Generating a locus_tag: {locus_tag}")

    note = qualifiers.get("note", [None])[0]
    locus = locus_tag

    # pseudo

    protein_id = qualifiers.get("protein_id", [None])[0]

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
    db_xref = qualifiers.get("db_xref", [so.SO_CDS.id])

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

# Feature Key           CDS

# Definition            coding sequence; sequence of nucleotides that
#                       corresponds with the sequence of amino acids in a
#                       protein (location includes stop codon); 
#                       feature includes amino acid conceptual translation.

# Optional qualifiers   /allele="text"
#                       /artificial_location="[artificial_location_value]"
#                       /circular_RNA
#                       /codon_start=<1 or 2 or 3>
#                       /db_xref="<database>:<identifier>"
#                       /EC_number="text"
#                       /exception="[exception_value]"
#                       /experiment="[CATEGORY:]text"
#                       /function="text"
#                       /gene="text"
#                       /gene_synonym="text"
#                       /inference="[CATEGORY:]TYPE[ (same species)][:EVIDENCE_BASIS]"
#                       /locus_tag="text" (single token)
#                       /map="text"
#                       /note="text"
#                       /number=unquoted text (single token)
#                       /old_locus_tag="text" (single token)
#                       /operon="text"
#                       /product="text"
#                       /protein_id="<identifier>"
#                       /pseudo
#                       /pseudogene="TYPE"
#                       /ribosomal_slippage
#                       /standard_name="text"
#                       /translation="text"
#                       /transl_except=(pos:<location>,aa:<amino_acid>)
#                       /transl_table =<integer>
#                       /trans_splicing

    multi_valued = {"EC_number", "exception", "experiment", "function",  "gene_synonym",  "inference", }
    single_valued = {"allele", "artificial_location",  "map", "number",  "old_locus_tag", "operon", "phenotype", "pseudogene", "standard_name", "transl_except", "transl_table"}

    add_optional_qualifiers(bakta_cds, qualifiers, single_valued, multi_valued)

    # Flags (boolean-like qualifiers)
    for flag in ["circular_RNA", "pseudo", "ribosomal_slippage", "trans_splicing"]:
        if flag in qualifiers:
            bakta_cds[flag] = flag in qualifiers

    if hypothetical:
        bakta_cds["hypothetical"] = True

    return bakta_cds

def convert_trna_feature(feature, seq_record, id):
    """
    Convert a funannotate tRNA SeqFeature to a Bakta tRNA JSON entry.
    """

    # ------------ Location ------------

    strand = "+" if feature.location.strand == 1 else "-"
    start = int(feature.location.start) + 1   # Bakta uses 1-based inclusive
    stop  = int(feature.location.end)         # already inclusive after conversion

    # Handle CompoundLocation (join)
    starts = None
    stops = None

    if feature.location.__class__.__name__ == "CompoundLocation":
        starts = []
        stops = []
        for part in feature.location.parts:
            starts.append(int(part.start) + 1)
            stops.append(int(part.end))



    # ------------ Extract nt sequence ------------
    nt_seq = feature.extract(seq_record.seq)
    nt = str(nt_seq)

    # ------------ Basic qualifiers ------------
    product = feature.qualifiers.get("product", [None])[0]

    # fall back to start_stop_strand if there is no locus tag
    if 'locus_tag' in qualifiers and qualifiers['locus_tag']:
        locus_tag = qualifiers['locus_tag'][0]
    else:
        logger.warning(f"No locus_tag found for feature {id}")
        locus_tag = f"{GENOME_RANDOM_BACKUP_LOCUSTAG_STR}_{start}_{stop}"
        logger.warning(f"Generating a locus_tag: {locus_tag}")

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
    bakta_trna_entry = {
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
        "locus": locus_tag,
        "id": id,
    }

# Feature Key           tRNA


# Definition            mature transfer RNA, a small RNA molecule (75-85 bases
#                       long) that mediates the translation of a nucleic acid
#                       sequence into an amino acid sequence;

# Optional qualifiers   /allele="text"
#                       /anticodon=(pos:<location>,aa:<amino_acid>,seq:<text>)
#                       /circular_RNA
#                       /db_xref="<database>:<identifier>"
#                       /experiment="[CATEGORY:]text"
#                       /function="text"
#                       /gene="text"
#                       /gene_synonym="text"
#                       /inference="[CATEGORY:]TYPE[ (same species)][:EVIDENCE_BASIS]"
#                       /locus_tag="text" (single token)
#                       /map="text"
#                       /note="text"
#                       /old_locus_tag="text" (single token)
#                       /operon="text"
#                       /product="text"
#                       /pseudo
#                       /pseudogene="TYPE"
#                       /standard_name="text"
#                       /trans_splicing

    multi_valued = {"experiment", "function",  "gene_synonym",  "inference", "note" }
    single_valued = {"allele",  "map",    "old_locus_tag", "standard_name"}

    qualifiers = feature.qualifiers

    add_optional_qualifiers(bakta_trna_entry, qualifiers, single_valued, multi_valued)

    # Flags (boolean-like qualifiers)
    for flag in ["circular_RNA", "pseudo", "trans_splicing"]:
        if flag in qualifiers:
            bakta_trna_entry[flag] = flag in qualifiers

    return bakta_trna_entry

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

    # Extract location
    strand = "+" if feature.location.strand == 1 else "-"
    start = int(feature.location.start) + 1   # Bakta uses 1-based inclusive
    stop  = int(feature.location.end)         # already inclusive after conversion


    qualifiers = feature.qualifiers

    # fall back to start_stop_strand if there is no locus tag
    if 'locus_tag' in qualifiers and qualifiers['locus_tag']:
        locus_tag = qualifiers['locus_tag'][0]
    else:
        logger.warning(f"No locus_tag found for feature {id}")
        locus_tag = f"{GENOME_RANDOM_BACKUP_LOCUSTAG_STR}_{start}_{stop}"
        logger.warning(f"Generating a locus_tag: {locus_tag}")
        
    

    gene_entry = {
        "type": "gene",
        "sequence": rec.id,
        "start": start,
        "stop": stop,
        "strand": strand,
        "db_xrefs": [so.SO_GENE.id], 
        "id": id,
        "locus": locus_tag
    }


# Feature Key           gene 


# Definition            region of biological interest identified as a gene 
#                       and for which a name has been assigned;

# Optional qualifiers   /allele="text"
#                       /db_xref="<database>:<identifier>"
#                       /experiment="[CATEGORY:]text"
#                       /function="text"
#                       /gene="text"
#                       /gene_synonym="text"
#                       /inference="[CATEGORY:]TYPE[ (same species)][:EVIDENCE_BASIS]"
#                       /locus_tag="text" (single token)
#                       /map="text"
#                       /note="text"
#                       /old_locus_tag="text" (single token)
#                       /operon="text"
#                       /product="text"
#                       /pseudo
#                       /pseudogene="TYPE"
#                       /phenotype="text"
#                       /standard_name="text"
#                       /trans_splicing

        
# Comment               the gene feature describes the interval of DNA that 
#                       corresponds to a genetic trait or phenotype; the feature is,
#                       by definition, not strictly bound to it's positions at the 
#                       ends;  it is meant to represent a region where the gene is 
#                       located.
    

    multi_valued = {"experiment", "function",  "gene_synonym",  "inference", "note" }
    single_valued = {"allele",  "map",  "old_locus_tag", "operon", "phenotype", "standard_name"}

    qualifiers = feature.qualifiers

    add_optional_qualifiers(gene_entry, qualifiers, single_valued, multi_valued)

    # Flags (boolean-like qualifiers)
    for flag in ["pseudo", "trans_splicing"]:
        if flag in qualifiers:
            gene_entry[flag] = flag in qualifiers

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
        dict: Bakta-style mRNA feature.
    """

    # seq = str(rec.seq)

    # Extract location
    strand = "+" if feature.location.strand == 1 else "-"
    start = int(feature.location.start) + 1   # Bakta uses 1-based inclusive
    stop  = int(feature.location.end)         # already inclusive after conversion

    # Handle CompoundLocation (join)
    starts = None
    stops = None

    if feature.location.__class__.__name__ == "CompoundLocation":
        starts = []
        stops = []
        for part in feature.location.parts:
            starts.append(int(part.start) + 1)
            stops.append(int(part.end))

    else:
        starts = None
        stops = None


    qualifiers = feature.qualifiers


    # fall back to start_stop_strand if there is no locus tag
    if 'locus_tag' in qualifiers and qualifiers['locus_tag']:
        locus_tag = qualifiers['locus_tag'][0]
    else:
        logger.warning(f"No locus_tag found for feature {id}")
        locus_tag = f"{GENOME_RANDOM_BACKUP_LOCUSTAG_STR}_{start}_{stop}"
        logger.warning(f"Generating a locus_tag: {locus_tag}")


    mrna_entry = {
        "type": "mRNA", 
        "sequence": rec.id,
        "start": start,
        "stop": stop,
        "starts": starts,
        "stops": stops,
        "strand": strand,
        "db_xrefs": [so.SO_MRNA.id],        
        "id": id,
        "locus": locus_tag
    }



# Feature Key           mRNA


# Definition            messenger RNA; includes 5'untranslated region (5'UTR),
#                       coding sequences (CDS, exon) and 3'untranslated region
#                       (3'UTR);

# Optional qualifiers   /allele="text"
#                       /artificial_location="[artificial_location_value]"
#                       /circular_RNA
#                       /db_xref="<database>:<identifier>"
#                       /experiment="[CATEGORY:]text"
#                       /function="text"
#                       /gene="text"
#                       /gene_synonym="text"
#                       /inference="[CATEGORY:]TYPE[ (same species)][:EVIDENCE_BASIS]"
#                       /locus_tag="text" (single token)
#                       /map="text"
#                       /note="text"
#                       /old_locus_tag="text" (single token)
#                       /operon="text"
#                       /product="text"
#                       /pseudo
#                       /pseudogene="TYPE"
#                       /standard_name="text"
#                       /trans_splicing


    multi_valued = {"experiment", "function",  "gene_synonym",  "inference", "note" }
    single_valued = {"allele", "artificial_location", "gene",  "locus_tag", "map",  "old_locus_tag", "operon", "phenotype", "product", "standard_name"}

    qualifiers = feature.qualifiers

    add_optional_qualifiers(mrna_entry, qualifiers, single_valued, multi_valued)

    # Flags (boolean-like qualifiers)
    for flag in ["circular_RNA", "pseudo", "trans_splicing"]:
        if flag in qualifiers:
            mrna_entry[flag] = flag in qualifiers

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
    strand = "." # bakta uses "." for strand on gaps
    start = int(feature.location.start) + 1   # Bakta uses 1-based inclusive
    stop  = int(feature.location.end)         # already inclusive after conversion

    qualifiers = feature.qualifiers

    #  may provide estimated_length but coordinates already give an exact span
    est_len = qualifiers.get("estimated_length", [None])[0]
    if est_len is not None:
        length = int(est_len)
    else:
        length = stop - start + 1  # fallback from coordinates


    gap_entry = {
        "type": "gap",
        "sequence": rec.id,
        "start": start,
        "stop": stop,
        "strand": strand,
        "length": length,
        "id": id,
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
    strand = "."
    start = int(feature.location.start) + 1   # Bakta uses 1-based inclusive
    stop  = int(feature.location.end)         # already inclusive after conversion


    qualifiers = feature.qualifiers
    note = qualifiers.get("note", [None])[0]
    rpt_family = qualifiers.get("rpt_family", [None])[0]
    rpt_type = qualifiers.get("rpt_type", [None])[0]
    rpt_unit_seq = qualifiers.get("rpt_unit_seq", [None])[0]

    # always just take the positive strand to get the NT seq (crispr repeat region)
    seq =  str(rec.seq)
    nt_seq = seq[start-1:stop]


# Feature Key           repeat_region


# Definition            region of genome containing repeating units;

# Optional qualifiers   /allele="text"
#                       /db_xref="<database>:<identifier>" 
#                       /experiment="[CATEGORY:]text"
#                       /function="text"
#                       /gene="text"
#                       /gene_synonym="text"
#                       /inference="[CATEGORY:]TYPE[ (same species)][:EVIDENCE_BASIS]"
#                       /locus_tag="text" (single token)
#                       /map="text"
#                       /note="text"
#                       /old_locus_tag="text" (single token)
#                       /rpt_family="text"
#                       /rpt_type=<repeat_type>
#                       /rpt_unit_range=<base_range>
#                       /rpt_unit_seq="text"
#                       /satellite="<satellite_type>[:<class>][ <identifier>]"
#                       /standard_name="text"

    so_code =  so.SO_REPEAT.id

    # Get existing db_xref list or default to [so.SO_CDS.id]
    db_xrefs = feature.qualifiers.get("db_xref", [so_code])

    # Append only if it’s not already present
    if so_code not in db_xrefs:
        db_xrefs.append(so_code)

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
        "db_xrefs": db_xrefs
    }

    multi_valued = {"experiment", "function",  "gene_synonym",  "inference", "note" }
    single_valued = {"satellite", "gene",  "locus_tag", "map",  "old_locus_tag", "operon", "phenotype", "product", "standard_name"}

    qualifiers = feature.qualifiers

    add_optional_qualifiers(repeat_region_entry, qualifiers, single_valued, multi_valued)


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

    # Extract location
    strand = "+" if feature.location.strand == 1 else "-"
    start = int(feature.location.start) + 1   # Bakta uses 1-based inclusive
    stop  = int(feature.location.end)         # already inclusive after conversion


    qualifiers = feature.qualifiers
    note = qualifiers.get("note", [None])[0]


    # fall back to start_stop_strand if there is no locus tag
    if 'locus_tag' in qualifiers and qualifiers['locus_tag']:
        locus_tag = qualifiers['locus_tag'][0]
    else:
        logger.warning(f"No locus_tag found for feature {id}")
        locus_tag = f"{GENOME_RANDOM_BACKUP_LOCUSTAG_STR}_{start}_{stop}"
        logger.warning(f"Generating a locus_tag: {locus_tag}")

    # always just take the positive strand to get the NT seq (UTR region)
    seq =  str(rec.seq)
    nt_seq = seq[start-1:stop]


# Feature Key           3'UTR


# Definition            1) region at the 3' end of a mature transcript (following 
#                       the stop codon) that is not translated into a protein;
#                       2) region at the 3' end of an RNA virus (following the last stop
#                       codon) that is not translated into a protein;

# Optional qualifiers   /allele="text"
#                       /db_xref="<database>:<identifier>"
#                       /experiment="[CATEGORY:]text"
#                       /function="text"
#                       /gene="text"
#                       /gene_synonym="text"
#                       /inference="[CATEGORY:]TYPE[ (same species)][:EVIDENCE_BASIS]"
#                       /locus_tag="text" (single token)
#                       /map="text"
#                       /note="text"
#                       /old_locus_tag="text" (single token)
#                       /standard_name="text"
#                       /trans_splicing



# Feature Key           5'UTR


# Definition            1) region at the 5' end of a mature transcript (preceding 
#                       the initiation codon) that is not translated into a protein;
#                       2) region at the 5' end of an RNA virus genome (preceding the first 
#                       initiation codon) that is not translated into a protein;

# Optional qualifiers   /allele="text"
#                       /db_xref="<database>:<identifier>"
#                       /experiment="[CATEGORY:]text"
#                       /function="text"
#                       /gene="text"
#                       /gene_synonym="text"
#                       /inference="[CATEGORY:]TYPE[ (same species)][:EVIDENCE_BASIS]"
#                       /locus_tag="text" (single token)
#                       /map="text"
#                       /note="text"
#                       /old_locus_tag="text" (single token)
#                       /standard_name="text"
#                       /trans_splicing

    so_code =  so.SO_REPEAT.id

    # Get existing db_xref list or default to [so.SO_CDS.id]
    db_xrefs = feature.qualifiers.get("db_xref", [so_code])

    # Append only if it’s not already present
    if so_code not in db_xrefs:
        db_xrefs.append(so_code)


    # Minimal Bakta-like structure
    utr_entry = {
        "type": type,
        "sequence": rec.id,
        "start": start,
        "stop": stop,
        "strand": strand, # matches Bakta and is required
        "product": note, 
        "nt": nt_seq, # needed for batka .ffn writeout
        "id": id, # bakta_id needed 
        "db_xrefs": db_xrefs,
        "locus": locus_tag
    }

    multi_valued = {"experiment", "function",  "gene_synonym",  "inference", "note" }
    single_valued = {"allele", "gene",   "map",  "old_locus_tag", "operon", "phenotype", "standard_name"}

    qualifiers = feature.qualifiers

    add_optional_qualifiers(utr_entry, qualifiers, single_valued, multi_valued)

    # Flags (boolean-like qualifiers)
    for flag in ["pseudo", "trans_splicing"]:
        if flag in qualifiers:
            utr_entry[flag] = flag in qualifiers

    return utr_entry

    #  3'UTR           186716..186779
    #                  /locus_tag="C1SCF055_LOCUS2256"
    #                  /note="ID:SCF055_s1087_g24103.utr3p1;
    #                  source:feature"

    #  5'UTR           218227..218248
    #                  /locus_tag="C1SCF055_LOCUS2213"
    #                  /note="ID:SCF055_s1083_g24060.utr5p1;
    #                  source:feature"

def convert_misc_rna_feature(feature, rec, id):
    """
    Convert a GenBank misc_rna feature to a simplified Bakta-style 'misc_rna' feature.

    Parameters:
        feature: Bio.SeqFeature
            The assembly_gap feature from the GBK.
        rec: Bio.SeqRecord
            The full GenBank record containing the sequence.

    Returns:
        dict: Simplified Bakta-style gap feature.

    """

        # from ensemble genomes
        # misc_RNA        complement(437333..442742)
        #             /gene="YPL060C-A"
        #             /note="transposable_element"
        #             /standard_name="YPL060C-A"

    # Coordinates (1-based)

    # Extract location
    strand = "+" if feature.location.strand == 1 else "-"
    start = int(feature.location.start) + 1   # Bakta uses 1-based inclusive
    stop  = int(feature.location.end)         # already inclusive after conversion


    qualifiers = feature.qualifiers
    gene = qualifiers.get("gene", [None])[0]

# Feature Key           misc_RNA


# Definition            any transcript or RNA product that cannot be defined by
#                       other RNA keys (prim_transcript, precursor_RNA, mRNA,
#                       5'UTR, 3'UTR, exon, CDS, sig_peptide, transit_peptide,
#                       mat_peptide, intron, polyA_site, ncRNA, rRNA and tRNA);

# Optional qualifiers   /allele="text"
#                       /db_xref="<database>:<identifier>"
#                       /experiment="[CATEGORY:]text"
#                       /function="text"
#                       /gene="text"
#                       /gene_synonym="text"
#                       /inference="[CATEGORY:]TYPE[ (same species)][:EVIDENCE_BASIS]"
#                       /locus_tag="text" (single token)
#                       /map="text"
#                       /note="text"
#                       /old_locus_tag="text" (single token)
#                       /operon="text"
#                       /product="text"
#                       /pseudo
#                       /pseudogene="TYPE"
#                       /standard_name="text"
#                       /trans_splicing

    misc_rna_entry = {
        "type": "misc_RNA", # expects lowercase 
        "sequence": rec.id,
        "start": start,
        "stop": stop,
        "strand": strand, # matches Bakta and is required
        "gene": gene,
        "id": id
    }


    multi_valued = {"experiment", "function",  "gene_synonym",  "inference", "note" }
    single_valued = {"allele", "gene",  "locus_tag", "map",  "old_locus_tag", "operon", "product", "phenotype", "standard_name"}

    qualifiers = feature.qualifiers

    add_optional_qualifiers(misc_rna_entry, qualifiers, single_valued, multi_valued)

    # Flags (boolean-like qualifiers)
    for flag in ["pseudo", "trans_splicing"]:
        if flag in qualifiers:
            misc_rna_entry[flag] = flag in qualifiers

    return misc_rna_entry

def convert_exon_feature(feature, rec, id):
    """
    Convert a GenBank exon feature to a simplified Bakta-style 'exon' feature.

    Parameters
    ----------
    feature : Bio.SeqFeature
        The exon feature from the GenBank record.
    rec : Bio.SeqRecord
        The full GenBank record.
    id : str
        Unique feature ID.

    Returns
    -------
    dict
        Bakta-style exon feature.
    """

    # Extract location
    strand = "+" if feature.location.strand == 1 else "-"
    start = int(feature.location.start) + 1   # Bakta uses 1-based inclusive
    stop  = int(feature.location.end)         # already inclusive after conversion

    qualifiers = feature.qualifiers

    db_xrefs = qualifiers.get("db_xref", [])

# Optional qualifiers   /allele="text"
#                       /db_xref="<database>:<identifier>"
#                       /EC_number="text"
#                       /experiment="[CATEGORY:]text"
#                       /function="text"
#                       /gene="text"
#                       /gene_synonym="text"
#                       /inference="[CATEGORY:]TYPE[ (same species)][:EVIDENCE_BASIS]"
#                       /locus_tag="text" (single token)
#                       /map="text"
#                       /note="text"
#                       /number=unquoted text (single token)
#                       /old_locus_tag="text" (single token)
#                       /product="text"
#                       /pseudo
#                       /pseudogene="TYPE"
#                       /standard_name="text"
#                       /trans_splicing


    # Extract commonly used INSDC qualifiers
    exon_entry = {
            "type": "exon",
            "sequence": rec.id,
            "start": start,
            "stop": stop,
            "strand": strand,
            "id": id,
            "db_xrefs": db_xrefs
        }
    
    multi_valued = {"EC_number","experiment","function",  "gene_synonym",  "inference","note" }
    single_valued = {"allele", "gene", "locus_tag", "map", "number",   "old_locus_tag", "operon", "pseudogene", "standard_name"   }

    add_optional_qualifiers(exon_entry, qualifiers, single_valued, multi_valued)

    # Flags (boolean-like qualifiers)
    for flag in ["pseudo", "trans_splicing"]:
        if flag in qualifiers:
            exon_entry[flag] = True

    return exon_entry


def convert_mat_peptide_feature(feature, rec, id):
    """
    Convert a mat_peptide feature to a Bakta-style feature.

    mus musculus chrom 1 NC_000067

    Parameters:
        feature: Bio.SeqFeature
            The mRNA feature from the GBK.
        rec: Bio.SeqRecord
            The record containing the sequence.

    Returns:
        dict: Bakta-style misc_RNA feature.
    """

    seq = str(rec.seq)


    # Extract location
    strand = "+" if feature.location.strand == 1 else "-"
    start = int(feature.location.start) + 1   # Bakta uses 1-based inclusive
    stop  = int(feature.location.end)         # already inclusive after conversion

    # Handle CompoundLocation (join)
    starts = None
    stops = None

    if feature.location.__class__.__name__ == "CompoundLocation":
        starts = []
        stops = []
        for part in feature.location.parts:
            starts.append(int(part.start) + 1)
            stops.append(int(part.end))

    so_code =  so.SO_MAT_PEPTIDE.id

    # Get existing db_xref list or default to [so.SO_CDS.id]
    db_xrefs = feature.qualifiers.get("db_xref", [so_code])

    # Append only if it’s not already present
    if so_code not in db_xrefs:
        db_xrefs.append(so_code)

    qualifiers = feature.qualifiers


    # Extract commonly used INSDC qualifiers
    mat_peptide_entry = {
            "type": "mat_peptide",
            "sequence": rec.id,
            "start": start,
            "stop": stop,
            # Join support
            "starts": starts,
            "stops": stops,
            "strand": strand,
            "id": id,
            "db_xrefs": db_xrefs
        }
   

# Feature Key           mat_peptide


# Definition            mature peptide or protein coding sequence; coding
#                       sequence for the mature or final peptide or protein
#                       product following post-translational modification; the
#                       location does not include the stop codon (unlike the
#                       corresponding CDS);

# Optional qualifiers   /allele="text"
#                       /db_xref="<database>:<identifier>"
#                       /EC_number="text"
#                       /experiment="[CATEGORY:]text"
#                       /function="text"
#                       /gene="text"
#                       /gene_synonym="text"
#                       /inference="[CATEGORY:]TYPE[ (same species)][:EVIDENCE_BASIS]"
#                       /locus_tag="text" (single token)
#                       /map="text"
#                       /note="text"
#                       /old_locus_tag="text" (single token)
#                       /product="text"
#                       /pseudo
#                       /pseudogene="TYPE"
#                       /standard_name="text"

    multi_valued = {"EC_number","experiment", "function",  "gene_synonym",  "inference","note" }
    single_valued = {"allele", "gene", "locus_tag", "map", "number",   "old_locus_tag", "operon", "pseudogene", "standard_name"}

    add_optional_qualifiers(mat_peptide_entry, qualifiers, single_valued, multi_valued)

    # Flags (boolean-like qualifiers) - no flags
    # for flag in ["pseudo"]:
    #     if flag in qualifiers:
    #         mat_peptide_entry[flag] = True


    #  mat_peptide     complement(join(194724303..194724321,194744661..194744721,
    #                  194746996..194747031,194750435..194750476,
    #                  194757818..194757865,194759962..194760144,
    #                  194764890..194765087,194765856..194765944,
    #                  194767641..194767743,194768400..194768583))
    #                  /gene="Cd46"
    #                  /gene_synonym="Mcp"
    #                  /product="Membrane cofactor protein. /id=PRO_0000238971"
    #                  /note="propagated from UniProtKB/Swiss-Prot (O88174.1)"

    return mat_peptide_entry

def convert_mobile_element_feature(feature, rec, id):
    """
    Convert a GenBank mobile_element feature to a Bakta-style feature.
    """

    # Extract location
    strand = "+" if feature.location.strand == 1 else "-"
    start = int(feature.location.start) + 1   # Bakta uses 1-based inclusive
    stop  = int(feature.location.end)         # already inclusive after conversion

    qualifiers = feature.qualifiers

    # Mandatory qualifier check (INSDC requirement)
    mobile_element_type = qualifiers.get("mobile_element_type", [None])[0]
    if mobile_element_type is None:
        raise ValueError(
            f"mobile_element feature {id} is missing mandatory "
            "/mobile_element_type qualifier"
        )

    so_code =  so.SO_MOBILE_ELEMENT.id

    # Get existing db_xref list or default to [so.SO_CDS.id]
    db_xrefs = feature.qualifiers.get("db_xref", [so_code])

    # Append only if it’s not already present
    if so_code not in db_xrefs:
        db_xrefs.append(so_code)


# Feature Key           mobile_element


# Definition            region of genome containing mobile elements;

# Mandatory qualifiers  /mobile_element_type="<mobile_element_type>
#                       [:<mobile_element_name>]"

# Optional qualifiers   /allele="text"
#                       /db_xref="<database>:<identifier>" 
#                       /experiment="[CATEGORY:]text"
#                       /function="text"
#                       /gene="text"
#                       /gene_synonym="text"
#                       /inference="[CATEGORY:]TYPE[ (same species)][:EVIDENCE_BASIS]"
#                       /locus_tag="text" (single token)
#                       /map="text"
#                       /note="text"
#                       /old_locus_tag="text" (single token)
#                       /rpt_family="text"
#                       /rpt_type=<repeat_type>
#                       /standard_name="text"


    # Extract commonly used INSDC qualifiers
    mobile_element_entry = {
            "type": "mobile_element",
            "sequence": rec.id,
            "start": start,
            "stop": stop,
            "strand": strand,
            "id": id,
            "db_xrefs": db_xrefs,
                    # Mandatory
            "mobile_element_type": mobile_element_type,
        }
    
    multi_valued = {"experiment", "function",  "gene_synonym",  "inference", "note" }
    single_valued = {"allele", "gene", "locus_tag", "map",    "old_locus_tag", "standard_name", "rpt_family", "rpt_type"}

    add_optional_qualifiers(mobile_element_entry, qualifiers, single_valued, multi_valued)

    # Flags (boolean-like qualifiers)
    # for flag in ["pseudo"]:
    #   if flag in qualifiers:
    #     mobile_element_entry[flag] = True


    #  mobile_element  57369551..57369723
    #                  /note="Derived by automated computational analysis using
    #                  gene prediction method: RefSeqFE."
    #                  /mobile_element_type="SINE:AmnSINE1"
    #                  /db_xref="GeneID:106707176"

    return mobile_element_entry


def convert_ncrna_feature(feature, rec, id):
    """
    Convert a ncrna feature to a Bakta-style feature.

    Parameters:
        feature: Bio.SeqFeature
            The mRNA feature from the GBK.
        rec: Bio.SeqRecord
            The record containing the sequence.

    Returns:
        dict: Bakta-style misc_RNA feature.
    """

    # seq = str(rec.seq)

    # Extract location
    strand = "+" if feature.location.strand == 1 else "-"
    start = int(feature.location.start) + 1   # Bakta uses 1-based inclusive
    stop  = int(feature.location.end)         # already inclusive after conversion

    # Handle CompoundLocation (join)
    starts = None
    stops = None

    if feature.location.__class__.__name__ == "CompoundLocation":
        starts = []
        stops = []
        for part in feature.location.parts:
            starts.append(int(part.start) + 1)
            stops.append(int(part.end))

    qualifiers = feature.qualifiers

    so_code =  so.SO_NCRNA_GENE.id

    # Get existing db_xref list or default to [so.SO_CDS.id]
    db_xrefs = feature.qualifiers.get("db_xref", [so_code])

    # Append so only if it’s not already present
    if so_code not in db_xrefs:
        db_xrefs.append(so_code)

    # Mandatory qualifier (INSDC requirement)
    ncrna_class = qualifiers.get("ncRNA_class", [None])[0]
    if ncrna_class is None:
        raise ValueError(
            f"ncRNA feature {id} is missing mandatory /ncRNA_class qualifier"
        )

    ncrna_entry = {
        "type": "ncRNA",
        "sequence": rec.id,
        "start": start,
        "stop": stop,
        "strand": strand,
        "id": id,

        # Join support
        "starts": starts,
        "stops": stops,

        # Mandatory
        "ncRNA_class": ncrna_class,

        # Multi-valued qualifiers
        "db_xrefs": db_xrefs,


    }

# Feature Key           ncRNA

# Definition            a non-protein-coding gene, other than ribosomal RNA and
#                       transfer RNA, the functional molecule of which is the RNA
#                       transcript;

# Mandatory qualifiers  /ncRNA_class="TYPE"
                      
# Optional qualifiers   /allele="text"
#                       /db_xref="<database>:<identifier>"
#                       /experiment="[CATEGORY:]text"
#                       /function="text"
#                       /gene="text"
#                       /gene_synonym="text"
#                       /inference="[CATEGORY:]TYPE[ (same species)][:EVIDENCE_BASIS]"
#                       /locus_tag="text" (single token)
#                       /map="text"
#                       /note="text"
#                       /old_locus_tag="text" (single token)
#                       /operon="text"
#                       /product="text"
#                       /pseudo
#                       /pseudogene="TYPE"
#                       /standard_name="text"
#                       /trans_splicing

# Example               /ncRNA_class="miRNA"
#                       /ncRNA_class="siRNA"
#                       /ncRNA_class="scRNA"       

# Comment               the ncRNA feature is not used for ribosomal and transfer
#                       RNA annotation, for which the rRNA and tRNA feature keys
#                       should be used, respectively;

    multi_valued = {"experiment", "function",  "gene_synonym",  "inference", "note" }
    single_valued = {"allele", "gene", "locus_tag", "map",    "old_locus_tag", "operon", "product", "standard_name", "pseudogene"}

    add_optional_qualifiers(ncrna_entry, qualifiers, single_valued, multi_valued)

    # Flags (boolean-like qualifiers)
    for flag in ["pseudo", "trans_splicing"]:
        if flag in qualifiers:
            ncrna_entry[flag] = flag in qualifiers

    #  ncRNA           join(189791085..189791793,189798997..189799081,
    #                  189819873..189820364,189821703..189822337)
    #                  /ncRNA_class="lncRNA"
    #                  /gene="Gm30446"
    #                  /product="predicted gene, 30446, transcript variant X6"
    #                  /note="Derived by automated computational analysis using
    #                  gene prediction method: Gnomon. Supporting evidence
    #                  includes similarity to: 100% coverage of the annotated
    #                  genomic feature by RNAseq alignments, including 2 samples
    #                  with support for all annotated introns"
    #                  /transcript_id="XR_001779629.1"
    #                  /db_xref="GeneID:102632350"
    #                  /db_xref="MGI:MGI:5589605"

    return ncrna_entry



def convert_misc_feature(feature, rec, id):
    """
    Convert a misc feature to a Bakta-style feature.

    Parameters:
        feature: Bio.SeqFeature
            The mRNA feature from the GBK.
        rec: Bio.SeqRecord
            The record containing the sequence.

    Returns:
        dict: Bakta-style misc_feature feature.
    """

    seq = str(rec.seq)

    # Extract location
    strand = "+" if feature.location.strand == 1 else "-"
    start = int(feature.location.start) + 1   # Bakta uses 1-based inclusive
    stop  = int(feature.location.end)         # already inclusive after conversion

    # Handle CompoundLocation (join)
    starts = None
    stops = None

    if feature.location.__class__.__name__ == "CompoundLocation":
        starts = []
        stops = []
        for part in feature.location.parts:
            starts.append(int(part.start) + 1)
            stops.append(int(part.end))

    qualifiers = feature.qualifiers

    so_code =  so.SO_MISC_REGION.id

    # Get existing db_xref list or default to [so.SO_CDS.id]
    db_xrefs = feature.qualifiers.get("db_xref", [so_code])

    # Append so.SO_CDS.id only if it’s not already present
    if so_code not in db_xrefs:
        db_xrefs.append(so_code)


    misc_feature_entry = {
            "type": "misc_feature",
            "sequence": rec.id,
            "start": start,
            "stop": stop,
            "strand": strand,
            "id": id,

            # Join support
            "starts": starts,
            "stops": stops,

            # Multi-valued
            "db_xrefs": db_xrefs,


        }

    multi_valued = {"experiment", "function",  "gene_synonym",  "inference", "note", "phenotype"}
    single_valued = {"allele", "gene", "locus_tag", "map",    "old_locus_tag", "operon", "product", "standard_name",  "pseudogene"}

    add_optional_qualifiers(misc_feature_entry, qualifiers, single_valued, multi_valued)

    # Flags (boolean-like qualifiers)
    for flag in ["pseudo",]:
        if flag in qualifiers:
            misc_feature_entry[flag] = True

# Feature Key           misc_feature


# Definition            region of biological interest which cannot be described
#                       by any other feature key; a new or rare feature;

# Optional qualifiers   /allele="text"
#                       /db_xref="<database>:<identifier>"
#                       /experiment="[CATEGORY:]text"
#                       /function="text"
#                       /gene="text"
#                       /gene_synonym="text"
#                       /inference="[CATEGORY:]TYPE[ (same species)][:EVIDENCE_BASIS]"
#                       /locus_tag="text" (single token)
#                       /map="text"
#                       /note="text"
#                       /number=unquoted text (single token)
#                       /old_locus_tag="text" (single token)
#                       /phenotype="text"
#                       /product="text"
#                       /pseudo
#                       /pseudogene="TYPE"
#                       /standard_name="text"

# Comment               this key should not be used when the need is merely to 
#                       mark a region in order to comment on it or to use it in 
#                       another feature's location

    #  misc_feature    join(78488668..78488692,78499322..78499359)
    #                  /gene="Mogat1"
    #                  /gene_synonym="0610030A14Rik; 1110064N14Rik; Dgat2l;
    #                  Dgat2l1; mDC2; MGAT1; WI1-2612I11.1"
    #                  /note="propagated from UniProtKB/Swiss-Prot (Q91ZV4.2);
    #                  transmembrane region"

    #  misc_feature    78179419..78180585
    #                  /standard_name="Pax3 upstream hypaxial enhancer"
    #                  /note="Region: biological region; Derived by automated
    #                  computational analysis using gene prediction method:
    #                  RefSeqFE."
    #                  /function="regulatory_interactions: LOC107980439 | Pax3"
    #                  /db_xref="GeneID:107980442"    

    return misc_feature_entry


def convert_proprotein_propeptide_feature(feature, rec, id):
    """
    Convert a proprotein or propeptide feature to a Bakta-style feature.

    Parameters:
        feature: Bio.SeqFeature
            The mRNA feature from the GBK.
        rec: Bio.SeqRecord
            The record containing the sequence.

    Returns:
        dict: Bakta-style proprotein feature.
    """

    # seq = str(rec.seq)

    # Extract location
    strand = "+" if feature.location.strand == 1 else "-"
    start = int(feature.location.start) + 1   # Bakta uses 1-based inclusive
    stop  = int(feature.location.end)         # already inclusive after conversion

    # Handle CompoundLocation (join)
    starts = None
    stops = None

    if feature.location.__class__.__name__ == "CompoundLocation":
        starts = []
        stops = []
        for part in feature.location.parts:
            starts.append(int(part.start) + 1)
            stops.append(int(part.end))


    qualifiers = feature.qualifiers

    so_code =  so.SO_PROPEPTIDE.id

    # Get existing db_xref list or default to [so.SO_CDS.id]
    db_xrefs = feature.qualifiers.get("db_xref", [so_code])

    # Append so.SO_CDS.id only if it’s not already present
    if so_code not in db_xrefs:
        db_xrefs.append(so_code)


    propeptide_entry = {
        "type": "propeptide",
        "sequence": rec.id,
        "start": start,
        "stop": stop,
        "strand": strand,
        "id": id,

        # Join support
        "starts": starts,
        "stops": stops,

        # Multi-valued
        "db_xrefs": qualifiers.get("db_xref", []),

    }


# Feature Key           propeptide


# Definition            propeptide coding sequence; coding sequence for the domain of a 
#                       proprotein that is cleaved to form the mature protein product.

# Optional qualifiers   /allele="text"
#                       /db_xref="<database>:<identifier>"
#                       /experiment="[CATEGORY:]text"
#                       /function="text"
#                       /gene="text"
#                       /gene_synonym="text"
#                       /inference="[CATEGORY:]TYPE[ (same species)][:EVIDENCE_BASIS]"
#                       /locus_tag="text" (single token)
#                       /map="text"
#                       /note="text"
#                       /old_locus_tag="text" (single token)
#                       /product="text"
#                       /pseudo
#                       /pseudogene="TYPE"
#                       /standard_name="text"


    multi_valued = {"experiment", "function",  "gene_synonym",  "inference", "note"}
    single_valued = {"allele", "gene", "locus_tag", "map",    "old_locus_tag", "product", "standard_name", "pseudogene"}

    add_optional_qualifiers(propeptide_entry, qualifiers, single_valued, multi_valued)

    # Flags (boolean-like qualifiers)
    for flag in ["pseudo",]:
        if flag in qualifiers:
            propeptide_entry[flag] = True

    #  proprotein      join(171053237..171053367,171053712..171053832)
    #                  /gene="Apoa2"
    #                  /gene_synonym="Alp-2; Apo-AII; Apoa-2; ApoA-II; ApoAII;
    #                  Hdl-1"
    #                  /product="apolipoprotein A-II proprotein"  

    return propeptide_entry



def convert_precursor_rna_feature(feature, rec, id):
    """
    Convert a GenBank precursor_RNA feature to a Bakta-style feature.
    """

    # Extract location
    strand = "+" if feature.location.strand == 1 else "-"
    start = int(feature.location.start) + 1   # Bakta uses 1-based inclusive
    stop  = int(feature.location.end)         # already inclusive after conversion

    qualifiers = feature.qualifiers

    so_code =  so.SO_PRECURSOR_RNA.id

    # Get existing db_xref list or default to [so.SO_CDS.id]
    db_xrefs = feature.qualifiers.get("db_xref", [so_code])

    # Append only if it’s not already present
    if so_code not in db_xrefs:
        db_xrefs.append(so_code)

    precursor_rna_entry = {
            "type": "precursor_RNA",
            "sequence": rec.id,
            "start": start,
            "stop": stop,
            "strand": strand,
            "db_xrefs": db_xrefs,
            "id": id,
        }


    multi_valued = {"experiment", "function",  "gene_synonym",  "inference", "note"}
    single_valued = {"allele", "gene", "locus_tag", "map",  "operon",  "old_locus_tag", "product", "standard_name"}

    add_optional_qualifiers(precursor_rna_entry, qualifiers, single_valued, multi_valued)

    # Flags (boolean-like qualifiers)
    for flag in ["trans_splicing"]:
        if flag in qualifiers:
            precursor_rna_entry[flag] = True

#     Feature Key           precursor_RNA


# Definition            any RNA species that is not yet the mature RNA product;
#                       may include ncRNA, rRNA, tRNA, 5' untranslated region
#                       (5'UTR), coding sequences (CDS, exon), intervening
#                       sequences (intron) and 3' untranslated region (3'UTR);

# Optional qualifiers   /allele="text"
#                       /db_xref="<database>:<identifier>"  
#                       /experiment="[CATEGORY:]text"
#                       /function="text"
#                       /gene="text"
#                       /gene_synonym="text"
#                       /inference="[CATEGORY:]TYPE[ (same species)][:EVIDENCE_BASIS]"
#                       /locus_tag="text" (single token)
#                       /map="text"
#                       /note="text"
#                       /old_locus_tag="text" (single token)
#                       /operon="text"
#                       /product="text"
#                       /standard_name="text"
#                       /trans_splicing


    #  precursor_RNA   194719348..194719428
    #                  /gene="Mir29b-2"
    #                  /gene_synonym="mir-29b-2; Mirn29b-2"
    #                  /product="microRNA 29b-2"
    #                  /note="Derived by automated computational analysis using
    #                  gene prediction method: BestRefSeq."
    #                  /transcript_id="NR_029809.1"
    #                  /db_xref="GeneID:723963"
    #                  /db_xref="MGI:MGI:3619047"
    #                  /db_xref="miRBase:MI0000712"

    return precursor_rna_entry


def convert_protein_bind_feature(feature, rec, id):
    """
    Convert a GenBank protein_bind feature to a Bakta-style feature.
    """

    # Extract location
    strand = "+" if feature.location.strand == 1 else "-"
    start = int(feature.location.start) + 1   # Bakta uses 1-based inclusive
    stop  = int(feature.location.end)         # already inclusive after conversion

    qualifiers = feature.qualifiers

    # Mandatory qualifier
    bound_moiety = qualifiers.get("bound_moiety", [None])[0]
    if bound_moiety is None:
        raise ValueError(
            f"protein_bind feature {id} is missing mandatory /bound_moiety qualifier"
        )

    so_code =  so.SO_PROTEINBIND.id

    # Get existing db_xref list or default to [so.SO_CDS.id]
    db_xrefs = feature.qualifiers.get("db_xref", [so_code])

    # Append only if it’s not already present
    if so_code not in db_xrefs:
        db_xrefs.append(so_code)

    protein_bind_entry = {
        "type": "protein_bind",
        "sequence": rec.id,
        "start": start,
        "stop": stop,
        "strand": strand,
        "bound_moiety": bound_moiety,
        "db_xrefs": db_xrefs,
        "id": id,
    }


# Feature Key           protein_bind


# Definition            non-covalent protein binding site on nucleic acid;

# Mandatory qualifiers  /bound_moiety="text"

# Optional qualifiers   /allele="text"
#                       /db_xref="<database>:<identifier>"
#                       /experiment="[CATEGORY:]text"
#                       /function="text"
#                       /gene="text"
#                       /gene_synonym="text"
#                       /inference="[CATEGORY:]TYPE[ (same species)][:EVIDENCE_BASIS]"
#                       /locus_tag="text" (single token)
#                       /map="text"
#                       /note="text"
#                       /old_locus_tag="text" (single token)
#                       /operon="text"
#                       /standard_name="text"

# Comment               note that feature key regulatory with /regulatory_class="ribosome_binding_site"
#                       should be used for ribosome binding sites.


    multi_valued = {"experiment", "function",  "gene_synonym",  "inference", "note"}
    single_valued = {"allele", "gene", "locus_tag", "map",  "operon",  "old_locus_tag", "product", "standard_name"}

    add_optional_qualifiers(protein_bind_entry, qualifiers, single_valued, multi_valued)

    # Flags (boolean-like qualifiers)
    # for flag in ["trans_splicing"]:
    #     if flag in qualifiers:
    #         protein_bind_entry[flag] = True

    return protein_bind_entry

def convert_rrna_feature(feature, rec, id):
    """
    Convert a GenBank rRNA feature to a Bakta-style feature.
    """

    # Extract location
    strand = "+" if feature.location.strand == 1 else "-"
    start = int(feature.location.start) + 1   # Bakta uses 1-based inclusive
    stop  = int(feature.location.end)         # already inclusive after conversion

    qualifiers = feature.qualifiers

    so_code =  so.SO_RRNA.id

    # Get existing db_xref list or default to [so.SO_CDS.id]
    db_xrefs = feature.qualifiers.get("db_xref", [so_code])

    # Append only if it’s not already present
    if so_code not in db_xrefs:
        db_xrefs.append(so_code)

# Feature Key           rRNA


# Definition            mature ribosomal RNA; RNA component of the
#                       ribonucleoprotein particle (ribosome) which assembles
#                       amino acids into proteins.

# Optional qualifiers   /allele="text"
#                       /db_xref="<database>:<identifier>"
#                       /experiment="[CATEGORY:]text"
#                       /function="text"
#                       /gene="text"
#                       /gene_synonym="text"
#                       /inference="[CATEGORY:]TYPE[ (same species)][:EVIDENCE_BASIS]"
#                       /locus_tag="text" (single token)
#                       /map="text"
#                       /note="text"
#                       /old_locus_tag="text" (single token)
#                       /operon="text"
#                       /product="text"
#                       /pseudo
#                       /pseudogene="TYPE"
#                       /standard_name="text"

# Comment               rRNA sizes should be annotated with the /product
#                       qualifier.  


    rrna_entry = {
            "type": "rRNA",
            "sequence": rec.id,
            "start": start,
            "stop": stop,
            "strand": strand,
            "db_xrefs": db_xrefs,
            "id": id,
        }
    
    multi_valued = {"experiment", "function",  "gene_synonym",  "inference", "note"}
    single_valued = {"allele", "gene", "locus_tag", "map",  "operon",  "old_locus_tag", "product", "pseudogene", "standard_name"}

    add_optional_qualifiers(rrna_entry, qualifiers, single_valued, multi_valued)

    # Flags (boolean-like qualifiers)
    for flag in ["pseudo"]:
        if flag in qualifiers:
            rrna_entry[flag] = True

    
    #  rRNA            46413357..46413475
    #                  /gene="n-R5s211"
    #                  /product="5S ribosomal RNA"
    #                  /inference="COORDINATES: nucleotide
    #                  motif:Rfam:12.0:RF00001"
    #                  /inference="COORDINATES: profile:INFERNAL:1.1.1"
    #                  /note="Derived by automated computational analysis using
    #                  gene prediction method: cmsearch."
    #                  /transcript_id="XR_004936691.1"
    #                  /db_xref="GeneID:115487577"
    #                  /db_xref="RFAM:RF00001"
    #                  /db_xref="MGI:MGI:4422076"

    return rrna_entry


def convert_regulatory_feature(feature, rec, id):
    """
    Convert a GenBank regulatory feature to a Bakta-style feature.
    """

    # Extract location
    strand = "+" if feature.location.strand == 1 else "-"
    start = int(feature.location.start) + 1   # Bakta uses 1-based inclusive
    stop  = int(feature.location.end)         # already inclusive after conversion

    qualifiers = feature.qualifiers

    # Mandatory qualifier
    regulatory_class = qualifiers.get("regulatory_class", [None])[0]
    if regulatory_class is None:
        raise ValueError(
            f"regulatory feature {id} is missing mandatory /regulatory_class qualifier"
        )

    so_code =  so.SO_REGULATORY_REGION.id

    # Get existing db_xref list or default to [so.SO_CDS.id]
    db_xrefs = feature.qualifiers.get("db_xref", [so_code])

    # Append only if it’s not already present
    if so_code not in db_xrefs:
        db_xrefs.append(so_code)


    regulatory_entry = {
            "type": "regulatory",
            "sequence": rec.id,
            "start": start,
            "stop": stop,
            "strand": strand,
            "regulatory_class": regulatory_class,
            "db_xrefs": db_xrefs,
            "id": id,
        }
    

# Feature Key           regulatory


# Definition            any region of sequence that functions in the regulation of
#                       transcription, translation, replication, recombination, or chromatin structure;

# Mandatory qualifiers  /regulatory_class="TYPE"

# Optional qualifiers   /allele="text"
#                       /bound_moiety="text"
#                       /db_xref="<database>:<identifier>"
#                       /experiment="[CATEGORY:]text"
#                       /function="text"
#                       /gene="text"
#                       /gene_synonym="text"
#                       /inference="[CATEGORY:]TYPE[ (same species)][:EVIDENCE_BASIS]"
#                       /locus_tag="text" (single token)
#                       /map="text"
#                       /note="text"
#                       /old_locus_tag="text" (single token)
#                       /operon="text"
#                       /phenotype="text"
#                       /pseudo
#                       /pseudogene="TYPE"
#                       /standard_name="text"

# Comment	              This feature has replaced the following Feature Keys on 15-DEC-2014:
#                       enhancer, promoter, CAAT_signal, TATA_signal, -35_signal, -10_signal,
#                       RBS, GC_signal, polyA_signal, attenuator, terminator, misc_signal.

    multi_valued = {"experiment", "function",  "gene_synonym",  "inference", "note"}
    single_valued = {"allele", "bound_moiety", "gene", "locus_tag", "map",  "operon",  "old_locus_tag", "phenotype", "product", "pseudogene", "standard_name"}

    add_optional_qualifiers(regulatory_entry, qualifiers, single_valued, multi_valued)

    # Flags (boolean-like qualifiers)
    for flag in ["pseudo"]:
        if flag in qualifiers:
            regulatory_entry[flag] = True

    #  regulatory      195030925..195032349
    #                  /regulatory_class="enhancer"
    #                  /experiment="EXISTENCE:reporter gene assay evidence
    #                  [ECO:0000049][PMID:32912294]"
    #                  /note="C2 STARR-seq-only enhancer starr_03508"
    #                  /function="activates a minimal SCP1 promoter by STARR-seq
    #                  in ground-state (2iL) and metastable (SL) mouse embryonic
    #                  stem cells {active_cell/tissue: mESC(E14 +2i+LIF or
    #                  +serum+LIF)}"
    #                  /db_xref="GeneID:131296982"

    return regulatory_entry



def convert_sig_peptide_feature(feature, rec, id):
    """
    Convert a sig_peptide feature to a Bakta-style feature.

    mus musculus chrom 1 NC_000067

    Parameters:
        feature: Bio.SeqFeature
            The mRNA feature from the GBK.
        rec: Bio.SeqRecord
            The record containing the sequence.

    Returns:
        dict: Bakta-style sig_peptide feature.
    """

    # Extract location
    strand = "+" if feature.location.strand == 1 else "-"
    start = int(feature.location.start) + 1   # Bakta uses 1-based inclusive
    stop  = int(feature.location.end)         # already inclusive after conversion

    # Handle CompoundLocation (join)
    starts = None
    stops = None

    if feature.location.__class__.__name__ == "CompoundLocation":
        starts = []
        stops = []
        for part in feature.location.parts:
            starts.append(int(part.start) + 1)
            stops.append(int(part.end))


    qualifiers = feature.qualifiers

    so_code =  so.SO_SIGNAL_PEPTIDE.id

    # Get existing db_xref list or default to [so.SO_CDS.id]
    db_xrefs = feature.qualifiers.get("db_xref", [so_code])

    # Append only if it’s not already present
    if so_code not in db_xrefs:
        db_xrefs.append(so_code)


    sig_peptide_entry = {
        "type": "sig_peptide",
        "sequence": rec.id,
        "start": start,
        "stop": stop,
        "strand": strand,
        "id": id,

        # Join support
        "starts": starts,
        "stops": stops,

        # Multi-valued
        "db_xrefs": qualifiers.get("db_xref", []),

    }

  


# Feature Key           sig_peptide


# Definition            signal peptide coding sequence; coding sequence for an
#                       N-terminal domain of a secreted protein; this domain is
#                       involved in attaching nascent polypeptide to the
#                       membrane leader sequence;

# Optional qualifiers   /allele="text"
#                       /db_xref="<database>:<identifier>"
#                       /experiment="[CATEGORY:]text"
#                       /function="text"
#                       /gene="text"
#                       /gene_synonym="text"
#                       /inference="[CATEGORY:]TYPE[ (same species)][:EVIDENCE_BASIS]"
#                       /locus_tag="text" (single token)
#                       /map="text"
#                       /note="text"
#                       /old_locus_tag="text" (single token)
#                       /product="text"
#                       /pseudo
#                       /pseudogene="TYPE"
#                       /standard_name="text"


    multi_valued = {"experiment", "function",  "gene_synonym",  "inference", "note"}
    single_valued = {"allele", "gene", "locus_tag", "map",  "operon",  "old_locus_tag", "phenotype", "product", "pseudogene", "standard_name"}

    add_optional_qualifiers(sig_peptide_entry, qualifiers, single_valued, multi_valued)

    # Flags (boolean-like qualifiers)
    for flag in ["pseudo"]:
        if flag in qualifiers:
            sig_peptide_entry[flag] = True


    #  sig_peptide     complement(join(194768584..194768588,
    #                  194774407..194774533))
    #                  /gene="Cd46"
    #                  /gene_synonym="Mcp"
    #                  /inference="COORDINATES: ab initio prediction:SignalP:6.0"

    return sig_peptide_entry


def convert_transit_peptide_feature(feature, rec, id):
    """
    Convert a transit_peptide feature to a Bakta-style feature.

    mus musculus chrom 1 NC_000067

    Parameters:
        feature: Bio.SeqFeature
            The mRNA feature from the GBK.
        rec: Bio.SeqRecord
            The record containing the sequence.

    Returns:
        dict: Bakta-style transit_peptide feature.
    """

    # seq = str(rec.seq)

    # Extract location
    strand = "+" if feature.location.strand == 1 else "-"
    start = int(feature.location.start) + 1   # Bakta uses 1-based inclusive
    stop  = int(feature.location.end)         # already inclusive after conversion

    # Handle CompoundLocation (join)
    starts = None
    stops = None

    if feature.location.__class__.__name__ == "CompoundLocation":
        starts = []
        stops = []
        for part in feature.location.parts:
            starts.append(int(part.start) + 1)
            stops.append(int(part.end))


    qualifiers = feature.qualifiers

    so_code =  so.SO_TRANSIT_PEPTIDE.id

    # Get existing db_xref list or default to [so.SO_CDS.id]
    db_xrefs = feature.qualifiers.get("db_xref", [so_code])

    # Append only if it’s not already present
    if so_code not in db_xrefs:
        db_xrefs.append(so_code)


    transit_peptide_entry = {
            "type": "transit_peptide",
            "sequence": rec.id,
            "start": start,
            "stop": stop,
            "strand": strand,
            "id": id,

            # Join support
            "starts": starts,
            "stops": stops,

            # Multi-valued
            "db_xrefs": qualifiers.get("db_xref", []),

        }



# Feature Key           transit_peptide


# Definition            transit peptide coding sequence; coding sequence for an
#                       N-terminal domain of a nuclear-encoded organellar
#                       protein; this domain is involved in post-translational
#                       import of the protein into the organelle;

# Optional qualifiers   /allele="text"
#                       /db_xref="<database>:<identifier>"
#                       /experiment="[CATEGORY:]text"
#                       /function="text"
#                       /gene="text"
#                       /gene_synonym="text"
#                       /inference="[CATEGORY:]TYPE[ (same species)][:EVIDENCE_BASIS]"
#                       /locus_tag="text" (single token)
#                       /map="text"
#                       /note="text"
#                       /old_locus_tag="text" (single token)
#                       /product="text"
#                       /pseudo
#                       /pseudogene="TYPE"
#                       /standard_name="text"

    multi_valued = {"experiment", "function",  "gene_synonym",  "inference", "note"}
    single_valued = {"allele", "gene", "locus_tag", "map",  "operon",  "old_locus_tag", "phenotype", "product", "pseudogene", "standard_name"}

    add_optional_qualifiers(transit_peptide_entry, qualifiers, single_valued, multi_valued)

    # Flags (boolean-like qualifiers)
    for flag in ["pseudo"]:
        if flag in qualifiers:
            transit_peptide_entry[flag] = True

    #  transit_peptide complement(join(180006550..180006849,
    #                  180009627..180009803))
    #                  /gene="Coq8a"
    #                  /gene_synonym="4632432J16Rik; Adck3; Cabc1; mKIAA0451"
    #                  /note="Mitochondrion.
    #                  /evidence=ECO:0000250|UniProtKB:Q8NI60; propagated from
    #                  UniProtKB/Swiss-Prot (Q60936.2)"

    return transit_peptide_entry

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
    mol_type = None
    organism = None
    strain = None
    db_xref = None
    note = None

    plasmid = None
    chromosome = None
    completeness_hint = None

    if source_feat:
        q = source_feat.qualifiers

        mol_type = q.get("mol_type", [None])[0]
        organism = q.get("organism", [None])[0]
        strain = q.get("strain", [None])[0]
        note = q.get("note", [None])[0]

        if "db_xref" in q:
            val = q["db_xref"]
            db_xref = val[0] if len(val) == 1 else val

        plasmid = q.get("plasmid", [None])[0]
        chromosome = q.get("chromosome", [None])[0]
        completeness_hint = q.get("completeness", [None])[0]
            
    # -----------------------------------------
    # Infer topology
    # -----------------------------------------
    topology = rec.annotations.get("topology")
    if topology not in {"linear", "circular"}:
        topology = "linear"

    # -----------------------------------------
    # Infer type
    # -----------------------------------------
    if plasmid is not None or "plasmid" in rec.annotations:
        seq_type = "plasmid"
    elif chromosome is not None or "chromosome" in rec.annotations:
        seq_type = "chromosome"
    else:
        seq_type = "contig"

    # -----------------------------------------
    # Infer completeness (conservative)
    # -----------------------------------------
    complete = False

    if topology == "circular":
        complete = True
    elif completeness_hint is not None and completeness_hint.lower() == "complete":
        complete = True
    elif note and "complete genome" in note.lower():
        complete = True

    # -----------------------------------------
    # Infer genetic codefor description
    # -----------------------------------------
    gcode = None

    if "genetic_code" in rec.annotations:
        gcode = rec.annotations["genetic_code"]
    elif "gcode" in rec.annotations:
        gcode = rec.annotations["gcode"]
    elif source_feat and "transl_table" in source_feat.qualifiers:
        gcode = source_feat.qualifiers["transl_table"][0]

    # Conservative fallback to 1 for euks
    if gcode is None:
        gcode = 1 

    description_parts = [
        f"[gcode={gcode}]",
        f"[topology={topology}]",
    ]

    description = " ".join(description_parts)

    # -----------------------------------------
    # Build entry
    # -----------------------------------------
    entry = {
        "id": rec.id,
        "description": description,
        "nt": seq,
        "length": len(seq),
        "complete": complete,
        "type": seq_type,
        "topology": topology,
        "simple_id": rec.id,
        "orig_id": rec.id,
        "orig_description": None,
    }

    # -----------------------------------------
    # Add source qualifiers if present
    # -----------------------------------------
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


    # this is from bakta
    # "id": "contig_1",
    # "description": "[gcode=11] [topology=linear]",
    # "nt": "AT"
    # "length": 5165988,
    # "complete": false,
    # "type": "contig",
    # "topology": "linear",
    # "simple_id": "contig_1",
    # "orig_id": "GCF_002368115_000000000001",
    # "orig_description": ""

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

def eukaryotic_gbk_to_json(records, output_json, verbose):
    
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

# https://www.insdc.org/submitting-standards/feature-table/
# eventually add them
    ORDER = ["tRNA", "gene", "mRNA", "CDS", "assembly_gap", "gap", "repeat_region", "5'UTR", "3'UTR", "misc_RNA", "exon",
             "mat_peptide", "mobile_element", "ncRNA", "misc_feature", "precursor_RNA", "proprotein",  # proprotein is actually not indsc compliant but is in Refseq
             "propeptide",
            "protein_bind", "rRNA",
             "regulatory", "sig_peptide", "transit_peptide"] 

     # source always in input - it is made in output anyway
    covered_set = set(ORDER + ["source"])

    # Compute features in records not in the covered list
    uncovered_features = unique_feature_types - covered_set

    if unique_feature_types:
        logger.info("Feature types present in your input GenBank and convertible include:")
        for ft in sorted(unique_feature_types):
            logger.info(ft)

    if uncovered_features:
        logger.warning("Feature types present in your input GenBank but NOT convertible include:")
        for ft in sorted(uncovered_features):
            logger.warning(ft)
        logger.warning("Baktfold will always only work on annotating CDS - however, these not converted feature types will not make it into the Baktfold output")
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

    # get total number of features

    total_features = 0

    for rec in records:

        total_features += len(rec.features)

    logger.info(f"Total {total_features} found")  

    # map feature type → list of features

    for rec in records:

        logger.info(f"Converting record: {rec.id} with {len(rec.features)} features")        

        # get a ftype, feature list
        features_by_type = defaultdict(list)
        for feat in rec.features:
            features_by_type[feat.type].append(feat)
        
        for ftype in ORDER:
            feats_by_type = features_by_type.get(ftype, [])
            num_feats = len(feats_by_type)
            if num_feats > 0:
                logger.info(f"Converting feature type: {ftype} with {num_feats} features") 

                # choose iterator based on verbose or not
                feat_iter = (
                    tqdm(
                        feats_by_type,
                        desc=f"{ftype}",
                        unit="feat",
                        leave=False,
                    )
                    if verbose
                    else feats_by_type
                )

                for feat in feat_iter:

                    id = f"{bakta_id_prefix}_{i}"

                    if verbose:
                        logger.info(f"Converting feature {id}")

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
                    elif ftype == "misc_RNA":
                        features.append(convert_misc_rna_feature(feat, rec, id))  
                    elif ftype == "exon":
                        features.append(convert_exon_feature(feat, rec, id))       
                    elif ftype == "mat_peptide":
                        features.append(convert_mat_peptide_feature(feat, rec, id))  
                    elif ftype == "mobile_element":
                        features.append(convert_mobile_element_feature(feat, rec, id))  
                    elif ftype == "ncRNA":
                        features.append(convert_ncrna_feature(feat, rec, id))  
                    elif ftype == "misc_feature":
                        features.append(convert_misc_feature(feat, rec, id))  
                    elif ftype == "precursor_RNA":
                        features.append(convert_precursor_rna_feature(feat, rec, id)) 
                    elif ftype == "proprotein":
                        features.append(convert_proprotein_propeptide_feature(feat, rec, id)) 
                    elif ftype == "propeptide":
                        features.append(convert_proprotein_propeptide_feature(feat, rec, id)) 
                    elif ftype == "protein_bind":
                        features.append(convert_protein_bind_feature(feat, rec, id)) 
                    elif ftype == "rRNA":
                        features.append(convert_rrna_feature(feat, rec, id)) 
                    elif ftype == "regulatory":
                        features.append(convert_regulatory_feature(feat, rec, id))
                    elif ftype == "sig_peptide":
                        features.append(convert_sig_peptide_feature(feat, rec, id))
                    elif ftype == "transit_peptide":
                        features.append(convert_transit_peptide_feature(feat, rec, id))
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

    logger.info(f"{total_features} total features converted successfully") 
    logger.info(f"Output saved as {output_json}")





