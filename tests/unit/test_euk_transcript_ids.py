"""Eukaryotic GFF3 hierarchy (audit N9): transcript IDs, shared CDS IDs, tRNAs.

Before: every transcript was '<locus>-T1' (duplicate IDs for alternative
transcripts), each CDS segment had its own ID ('.cds-1', '.cds-2', so GFF3
readers and table2asn saw one CDS per exon), every CDS pointed at '-T1', and
the tRNA writer discarded its own ID/Parent.

Run::

    pytest tests/unit/test_euk_transcript_ids.py -v
"""
from __future__ import annotations

import io as _io

import pytest

from baktfold.io import gff
from baktfold.io.euk_transcripts import gff_ids, transcript_id_of


@pytest.mark.parametrize(
    "feat, expected",
    [
        ({"type": "mRNA", "transcript_id": "XM_001234.1"}, "XM_001234.1"),
        ({"type": "mRNA", "standard_name": "ENST0001"}, "ENST0001"),
        ({"type": "mRNA", "note": ["transcript EHI_111330A"]}, "EHI_111330A"),
        ({"type": "cds", "note": "encoded by transcript EHI_111330A; predicted by EVM"}, "EHI_111330A"),
        ({"type": "cds", "protein_id": "ncbi:I1I15_012813-T2"}, "I1I15_012813-T2"),
        ({"type": "cds", "protein_id": "EDS88578.1"}, None),     # plain protein accession
        ({"type": "mRNA"}, None),
    ],
)
def test_transcript_id_of(feat, expected):
    assert transcript_id_of(feat) == expected


def _f(t, start, stop, parts=None, **kw):
    parts = parts or [(start, stop)]
    f = {"type": t, "locus": "g", "sequence": "chr", "strand": "+", "start": start, "stop": stop,
         "starts": [p[0] for p in parts], "stops": [p[1] for p in parts], "product": "x"}
    f.update(kw)
    return f


def test_alternative_transcripts_unique_ids_and_parents():
    gene = {"type": "gene", "locus": "g"}
    m3 = _f("mRNA", 100, 1000, [(100, 300), (400, 600), (700, 1000)])
    m2 = _f("mRNA", 100, 1000, [(100, 600), (700, 1000)])
    c3 = _f("cds", 150, 900, [(150, 300), (400, 600), (700, 900)])
    c2 = _f("cds", 150, 900, [(150, 600), (700, 900)])
    mrna_ids, cds_ids, _, _ = gff_ids([gene, m3, m2, c2, c3])
    assert len({mrna_ids[id(m3)], mrna_ids[id(m2)]}) == 2          # unique
    assert cds_ids[id(c3)][1] == mrna_ids[id(m3)]                   # own transcript
    assert cds_ids[id(c2)][1] == mrna_ids[id(m2)]
    assert cds_ids[id(c3)][0] == f"{mrna_ids[id(m3)]}.cds"


def test_input_transcript_ids_used_and_kept_unique():
    m1 = _f("mRNA", 100, 500, note="transcript TA001")
    m2 = _f("mRNA", 100, 500, note="transcript TA001")   # duplicated in the input
    mrna_ids, _, _, _ = gff_ids([m1, m2])
    assert mrna_ids[id(m1)] == "TA001"
    assert mrna_ids[id(m2)] == "TA001_2"


def test_cds_without_mrna_hangs_off_gene_only_if_it_exists():
    c = _f("cds", 100, 500)
    assert gff_ids([{"type": "gene", "locus": "g"}, c])[1][id(c)] == ("g.cds", "g")
    assert gff_ids([c])[1][id(c)] == ("g.cds", None)


def test_spliced_cds_written_with_one_id():
    c = _f("cds", 150, 900, [(150, 300), (400, 600), (700, 900)])
    m = _f("mRNA", 100, 1000, [(100, 300), (400, 600), (700, 1000)])
    ids = gff_ids([{"type": "gene", "locus": "g"}, m, c])
    fh = _io.StringIO()
    gff.write_euk_cds_feature(fh, "chr", c, ids)
    attrs = [l.split("\t")[8] for l in fh.getvalue().strip().split("\n")]
    assert len(attrs) == 3
    assert len({a.split(";")[0] for a in attrs}) == 1   # same ID on every segment
    assert all(f"Parent={ids[0][id(m)]}" in a for a in attrs)


def test_trna_has_id_parent_and_exon():
    t = {"type": "tRNA", "locus": "t1", "sequence": "chr", "strand": "-", "start": 10, "stop": 80,
         "product": "tRNA-Ala"}
    ids = gff_ids([{"type": "gene", "locus": "t1"}, t])
    fh = _io.StringIO()
    gff.write_euk_trna_feature(fh, "chr", t, ids)
    trna, exon = fh.getvalue().strip().split("\n")
    assert "ID=t1-T1;Parent=t1;product=tRNA-Ala" in trna
    assert "ID=t1-T1.exon1;Parent=t1-T1" in exon
    assert "ß" not in exon


def test_ensembl_gene_link_via_gene_qualifier():
    """Ensembl main/Plants: the gene has /locus_tag (a symbol) and /gene; the
    mRNA and CDS have only /gene, which the converter now uses as their locus."""
    gene = {"type": "gene", "locus": "aap-1", "gene": "WBGene00000001",
            "sequence": "chr", "strand": "+", "start": 100, "stop": 1000}
    m = _f("mRNA", 100, 1000, locus="WBGene00000001", standard_name="Y110A7A.10.1")
    c = _f("cds", 150, 900, locus="WBGene00000001", note="transcript_id=Y110A7A.10.1")
    mrna_ids, cds_ids, parents, _ = gff_ids([gene, m, c])
    assert mrna_ids[id(m)] == "Y110A7A.10.1"
    assert parents[id(m)] == "aap-1"                       # Parent = the gene's ID
    assert cds_ids[id(c)] == ("Y110A7A.10.1.cds", "Y110A7A.10.1")


def test_resolve_locus_prefers_locus_tag_then_gene():
    from baktfold.io.eukaryotic_to_json import resolve_locus
    assert resolve_locus({"locus_tag": ["L1"], "gene": ["G1"]}, "id", 1, 2) == "L1"
    assert resolve_locus({"gene": ["G1"]}, "id", 1, 2) == "G1"
    assert resolve_locus({}, "id", 1, 2).endswith("_1_2")


def test_two_genes_sharing_a_gene_name_get_unique_ids():
    """E. cuniculi (Genoscope) gives a protein gene and a tRNA gene the same
    /gene; each gene gets its own ID and each child the gene containing it."""
    g1 = {"type": "gene", "locus": "ECU11_0375", "sequence": "chr", "strand": "+", "start": 100, "stop": 300}
    g2 = {"type": "gene", "locus": "ECU11_0375", "sequence": "chr", "strand": "-", "start": 400, "stop": 480}
    t = {"type": "tRNA", "locus": "ECU11_0375", "sequence": "chr", "strand": "-", "start": 400, "stop": 480}
    _, _, parents, gene_ids = gff_ids([g1, g2, t])
    assert len({gene_ids[id(g1)], gene_ids[id(g2)]}) == 2
    assert parents[id(t)] == gene_ids[id(g2)]


def test_refseq_drosophila_pairs_via_note_alias():
    """mRNA /transcript_id is the RefSeq accession; the CDS names the FlyBase
    transcript, which the mRNA lists only in /note."""
    from baktfold.io.euk_transcripts import pair_transcripts
    ma = _f("mRNA", 100, 1000, transcript_id="NM_1.1", note="PlexB-RA; CG17245-RA")
    mb = _f("mRNA", 100, 1000, transcript_id="NM_2.1", note="PlexB-RB; CG17245-RB")
    cb = _f("cds", 150, 900, note="CG17245 gene product from transcript CG17245-RB")
    ca = _f("cds", 150, 900, note="CG17245 gene product from transcript CG17245-RA")
    cds_of, _, sure, _ = pair_transcripts([ma, mb, cb, ca])
    assert cds_of[id(ma)] is ca and cds_of[id(mb)] is cb      # identical CDS structure
    assert {id(ma), id(mb)} <= sure                           # paired by ID, not guessed
