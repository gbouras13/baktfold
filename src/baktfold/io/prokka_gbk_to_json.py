#!/usr/bin/env python3
# import argparse
import hashlib
import json
import random
import re
import string
from collections import defaultdict
from datetime import datetime

# from Bio import SeqIO
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
    start = int(feature.location.start) + 1     # Bakta uses 1-based inclusive
    stop  = int(feature.location.end)           # already one past in BioPython
    strand = "+" if feature.location.strand == 1 else "-"

    # frame: Bakta uses 1/2/3; Prokka codon_start is ["1","2","3"]
    codon_start = int(feature.qualifiers.get("codon_start", ["1"])[0])
    frame = codon_start

    # ----------- Basic qualifiers -----------
    gene = feature.qualifiers.get("gene", [None])[0]
    product = feature.qualifiers.get("product", [None])[0]

    locus_tag = feature.qualifiers.get("locus_tag", [None])[0]
    locus = locus_tag

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

    # ----------- Make Bakta-format dict -----------
    bakta_cds = {
        "type": "cds",
        "sequence": seq_record.id,
        "start": start,
        "stop": stop,
        "strand": strand,
        "frame": frame,
        "gene": gene,
        "product": product,
        "db_xrefs": feature.qualifiers.get("db_xref", [so.SO_CDS.id]),  # there will be no other db_xref 
        "nt": nt,
        "aa": aa,
        "aa_hexdigest": aa_hexdigest,
        "start_type": None,
        "rbs_motif": None,
        "genes": [],
        "seq_stats": seq_stats,
        "id": id,
        "locus": locus,
    }

    if hypothetical:
        bakta_cds["hypothetical"] = True

    return bakta_cds

def convert_trna_feature(feature, seq_record, id):
    """
    Convert a Prokka tRNA SeqFeature to a Bakta tRNA JSON entry.
    """

    # ------------ Location ------------
    start = int(feature.location.start) + 1
    stop  = int(feature.location.end)
    strand = "+" if feature.location.strand == 1 else "-"

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
       #  "anti_codon_pos": anti_codon_pos,  dont include, not in prokka output
        "id": id,
        "locus": locus
    }

    return bakta_trna

def convert_rrna_feature(feature, rec, id):
    """
    Convert a Prokka GenBank rRNA feature to Bakta-style JSON.
    
    Parameters:
        feature: Bio.SeqFeature
            The rRNA feature from the Prokka GBK.
        rec: str
            The record from the GBK.
    Returns:
        dict: Bakta-style rRNA feature
    """
    start = int(feature.location.start) + 1  # GBK is 0-based
    stop = int(feature.location.end)
    strand = "+" if feature.location.strand == 1 else "-"
    
    qualifiers = feature.qualifiers
    product = qualifiers.get("product", [None])[0]
    locus_tag = qualifiers.get("locus_tag", [None])[0]
    
    # Infer gene type from product if possible
    gene_map = {
        "16S ribosomal RNA": "rrs",
        "23S ribosomal RNA": "rrl",
        "5S ribosomal RNA": "rrf"
    }
    gene = gene_map.get(product, None)

    so_map = {
        "16S ribosomal RNA": so.SO_RRNA_16S.id,
        "23S ribosomal RNA": so.SO_RRNA_23S.id,
        "5S ribosomal RNA": so.SO_RRNA_5S.id
    }

    specific_so = so_map.get(product, None)

    contig_seq = str(rec.seq)
    
    nt_seq = contig_seq[start-1:stop]
    if strand == "-":
        comp = str.maketrans("ACGTacgt", "TGCAtgca")
        nt_seq = nt_seq.translate(comp)[::-1]
    
    rrna_entry = {
        "type": "rRNA",
        "sequence": rec.id,
        "start": start,
        "stop": stop,
        "strand": strand,
        "gene": gene,
        "product": product,
        "coverage": None,  # Prokka does not provide
        "score": None,     # Prokka does not provide
        "evalue": None,    # Prokka does not provide
        "db_xrefs": [so.SO_RRNA.id, specific_so], 
        "nt": nt_seq,
        "id": id,
        "locus": locus_tag
    }
    
    return rrna_entry

def convert_misc_rna_feature(feature, rec, id):
    """
    Convert a Prokka GenBank misc_RNA (nc_rna) feature to a Bakta-style feature.

    Parameters:
        feature: Bio.SeqFeature
            The misc_RNA feature from the Prokka GBK.
        rec: Bio.SeqRecord
            The record containing the sequence.

    Returns:
        dict: Bakta-style misc_RNA feature.
    """

    seq = str(rec.seq)

    # Coordinates (GBK is 0-based, Bakta is 1-based)
    start = int(feature.location.start) + 1
    stop = int(feature.location.end)
    strand = "+" if feature.location.strand == 1 else "-"

    qualifiers = feature.qualifiers

    # Fields that Prokka may or may not include
    gene = qualifiers.get("gene", [None])[0]
    product = qualifiers.get("product", [None])[0]
    locus_tag = qualifiers.get("locus_tag", [None])[0]

    # If gene missing, use product (Bakta often fills 'gene' for sRNA/tmRNA/etc.)
    if gene is None:
        gene = product

    # Extract nucleotide sequence
    nt_seq = seq[start-1:stop]
    if strand == "-":
        comp = str.maketrans("ACGTacgt", "TGCAtgca")
        nt_seq = nt_seq.translate(comp)[::-1]

    misc_entry = {
        "type": "ncRNA", # bakta uses ncRNA -> will be misc_rna in Prokka
        "class": None,             # bakta's classes are not in Prokka
        "sequence": rec.id,
        "start": start,
        "stop": stop,
        "strand": strand,
        "gene": gene,
        "product": product,
        "score": None,
        "evalue": None,
        "db_xrefs": [so.SO_NCRNA_GENE.id],        
        "nt": nt_seq,
        "id": id,
        "locus": locus_tag
    }

    return misc_entry

def convert_tmrna_feature(feature, rec, id):
    """
    Convert a Prokka GenBank tmRNA feature to Bakta-style feature.
    
    Parameters:
        feature: Bio.SeqFeature
            The tmRNA feature from the Prokka GBK.
        rec: str
            The record from the GBK
    
    Returns:
        dict: Bakta-style tmRNA feature
    """

    seq =  str(rec.seq)

    start = int(feature.location.start) + 1  # GBK is 0-based
    stop = int(feature.location.end)
    strand = "+" if feature.location.strand == 1 else "-"
    
    qualifiers = feature.qualifiers
    gene = qualifiers.get("gene", [None])[0]
    locus_tag = qualifiers.get("locus_tag", [None])[0]
    product = qualifiers.get("product", [None])[0]
    
    # Extract the nucleotide sequence of the tmRNA
    nt_seq = seq[start-1:stop]
    if strand == "-":
        comp = str.maketrans("ACGTacgt", "TGCAtgca")
        nt_seq = nt_seq.translate(comp)[::-1]
    

    
    tmrna_entry = {
        "type": "tmRNA",
        "sequence": rec.id,
        "start": start,
        "stop": stop,
        "strand": strand,
        "gene": gene,
        "product": product,
        "db_xrefs": [so.SO_TMRNA.id], 
        # "tag": tag_info,  no tag in tmrna for prokka - no information on it in the output
        "nt": nt_seq,
        "id": id,
        "locus": locus_tag
    }
    
    return tmrna_entry

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
    crispr_entry = {
        "type": "crispr",
        "sequence": rec.id,
        "start": start,
        "stop": stop,
        "strand": strand, # matches Bakta and is required
        "family": rpt_family,       # e.g., "CRISPR"
        "rpt_type": rpt_type,       # e.g., "direct"
        "repeat_unit": rpt_unit_seq, # the actual consensus repeat
        "product": note, # won't be the same as Bakta as different lookup method used - but needed for the gff writing
        "nt": nt_seq, # needed for batka .ffn writeout
        "id": id, # bakta_id needed 
        # "locus": None, # no locus tag like Bakta
        "db_xrefs": [so.SO_CRISPR.id]
    }

    return crispr_entry

def convert_assembly_gap_feature(feature, rec, id):
    """
    Convert a Prokka GenBank assembly_gap feature to a simplified Bakta-style 'gap' feature.

    Parameters:
        feature: Bio.SeqFeature
            The assembly_gap feature from the Prokka GBK.
        rec: Bio.SeqRecord
            The full GenBank record containing the sequence.

    Returns:
        dict: Simplified Bakta-style gap feature.
    """

    # Coordinates (1-based)
    start = int(feature.location.start) + 1
    stop = int(feature.location.end)

    qualifiers = feature.qualifiers

    # Prokka may provide estimated_length but coordinates already give an exact span
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

    return gap_entry


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

    # Conservative fallback to 11 for prokka
    if gcode is None:
        gcode = 11 

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



def parse_prokka_version(record):
    """
    Extract Prokka version from COMMENT field:
    Example COMMENT:
    'Annotated using prokka 1.14.6 ...'
    """
    
    comments = record.annotations.get("comment", "") or record.annotations.get("comments", "")
    if not comments:
        return "unknown"

    m = re.search(r"[Pp]rokka[\s_]?v?(\d+\.\d+\.\d+)", comments)
    if m:
        return m.group(1)

    # fallback pattern
    m = re.search(r"prokka[^0-9]*(\d+\.\d+\.\d+)", comments)
    if m:
        return m.group(1)

    return "unknown"


def get_transl_table(records):
    """
    Gets translation table based off the first CDS on the first record
    """

    if not records:
        raise ValueError("No GenBank records found.")
    
    record_1 = records[0]

    for feat in record_1.features:
        if feat.type == "CDS":
            # Translation table may be string → convert to int
            transl = feat.qualifiers.get("transl_table", ["11"])[0]
            try:
                return int(transl)
            except ValueError:
                return 11

        # If no CDS found, default to 11
        return 11

def random_n_letter_id(n=4):
    """
    generates a n letter id prefix 
    
    n=2 to append to  Prokka locus tag  for bakta id to make it different
    n=10 if the locus tag is somehow missing (should never happen) 
    """
    return ''.join(random.choices(string.ascii_uppercase, k=n))

def get_bakta_style_id_from_locus_tag(records):
    """
    Gets 10 char bakta-style ID tag based off the 8 char locus tag in first CDS on the first Prokka record + 2 random chars

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

                    if len(locus_tag) > 6:

                        locus_tag_prefix = locus_tag[:-6] # trims off _00001 from CDS

                        rand_two_chars = random_n_letter_id(2)

                        # by default prokka locus tag is 8 chars. So this returns a 10 char string (same as bakta defaults)

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
    Prokka GenBank file.
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

def prokka_gbk_to_json(records, output_json):
    """
    Convert Prokka-generated GenBank SeqRecord objects into a Bakta-style JSON
    annotation file.

    This function takes one or more Biopython SeqRecord objects (typically parsed
    from a Prokka GenBank file) and reconstructs a JSON structure following the
    Bakta output schema. It extracts genome metadata, statistics, annotated features,
    nucleotide sequences, and Prokka version information. Features are converted
    to Bakta-compatible dictionaries and sorted in the same order Bakta expects.

    The JSON file contains:
      • `genome` block – high-level organism metadata (genus, species, strain, etc.)
      • `stats` block – genome statistics derived from all records
      • `features` block – all features converted to Bakta-style objects, sorted
        by feature type and genomic position
      • `sequences` block – contig/sequence entries in Bakta format
      • `run` block – timestamps and duration placeholder
      • `version` block – Prokka version and database metadata

    Parameters
    ----------
    records : list of SeqRecord
        A list of Biopython SeqRecord objects already parsed from a Prokka GenBank
        file. Must contain at least one record. The COMMENT field is expected to
        contain Prokka metadata.

    output_json : str
        Path to the output JSON file to be written.

    Returns
    -------
    bool
        True if the JSON was successfully written.

    Raises
    ------
    ValueError
        If `records` is empty.

    Notes
    -----
    • Features are processed in a fixed Bakta-like order:
        ["tRNA", "tmRNA", "rRNA", "misc_RNA", "repeat_region", "CDS", "assembly_gap"]
    • Feature IDs are generated in Bakta-style using the locus tag prefix.
    • Per-contig sorting is performed by genomic start coordinate.
    • Runtime values in the `run` block are placeholders (duration = "0.00 min").
    • The function does not validate that the GenBank file truly originates from
      Prokka; that should be checked beforehand.
    """
    
    complete = False

    # records = list(SeqIO.parse(genbank_path, "genbank"))
    
    if len(records) == 0:
        raise ValueError("No GenBank records found.")
    elif len(records) >= 1:
        prokka_version = parse_prokka_version(records[0])

    translation_table = get_transl_table(records)

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

    # order by id creation block to match how ids are generated bakta
    ORDER = ["tRNA", "tmRNA", "rRNA", "misc_RNA", "repeat_region", "CDS", "assembly_gap"]
    # bakta has oriC detection too - prokka doesn't I think so leaving it out 

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
                elif ftype == "tRNA":
                    features.append(convert_trna_feature(feat, rec, id))
                elif ftype == "tmRNA":
                    features.append(convert_tmrna_feature(feat, rec, id))
                elif ftype == "rRNA":
                    features.append(convert_rrna_feature(feat, rec, id))
                elif ftype == "misc_RNA":
                    features.append(convert_misc_rna_feature(feat, rec, id))
                elif ftype == "repeat_region":
                    features.append(convert_repeat_region_feature(feat, rec, id))
                elif ftype == "assembly_gap":
                    features.append(convert_assembly_gap_feature(feat, rec, id))
                
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
            "prokka": prokka_version,  
            "db": {
                "version": prokka_version,
                "type": "prokka_dbs"
            }
        }
    }



    with open(output_json, "w") as fh:
        json.dump(bakta_json, fh, indent=4)
        complete = True

    return complete




