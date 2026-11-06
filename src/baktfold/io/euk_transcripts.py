"""Link eukaryotic CDS features to their mRNA (transcript), and name transcripts.

GenBank/EMBL input does not link a CDS to its mRNA explicitly, and different
annotation sources record transcript identifiers in different places:

* ``/transcript_id`` on the mRNA (NCBI RefSeq/GenBank) and sometimes the CDS;
* ``/standard_name`` on the mRNA;
* ``/note="transcript X"`` on the mRNA and ``/note="encoded by transcript X"``
  on the CDS (INSDC records produced by EVM-style pipelines, e.g. Entamoeba);
* the CDS ``/protein_id`` carrying the transcript name with a ``-T<n>`` suffix
  (funannotate: ``ncbi:LOCUS-T1``).

A CDS is paired with an mRNA by identical transcript identifier when both have
one; otherwise by structure: same locus tag, sequence and strand, the CDS lies
within the mRNA, and their splice junctions agree inside the CDS. Unambiguous
pairs are assigned first, repeatedly. Anything still ambiguous (e.g. nested
single-exon transcripts) is paired in feature order for the GFF3 hierarchy,
but flagged so labels are not copied across on a guess.
"""
from __future__ import annotations

import re
from typing import Dict, List, Optional, Sequence, Tuple

import baktfold.bakta.constants as bc

# "transcript EHI_111330A", "encoded by transcript X", "transcript_id=Y110A7A.10.1".
# The identifier must contain a digit, so RefSeq prose such as "transcript
# variant A" is not mistaken for one.
RE_NOTE_TRANSCRIPT = re.compile(r'\btranscript(?:_id)?[ :=]+(?=[A-Za-z0-9_.:\-]*\d)([A-Za-z0-9_.:\-]+?)[;,.]?(?:\s|$)')
RE_PROTEIN_ID_TRANSCRIPT = re.compile(r'^(?:[A-Za-z]+:)?(.+-T\d+)$')


def _as_text(value) -> str:
    if value is None:
        return ''
    if isinstance(value, (list, tuple)):
        return ' '.join(str(v) for v in value)
    return str(value)


def transcript_id_of(feat: dict) -> Optional[str]:
    """The input's own identifier for this feature's transcript, if any."""
    tid = _as_text(feat.get('transcript_id')).strip()
    if tid:
        return tid.split(':', 1)[1] if tid.startswith('ncbi:') else tid
    if feat.get('type') == bc.FEATURE_MRNA:
        std = _as_text(feat.get('standard_name')).strip()
        if std:
            return std
    m = RE_NOTE_TRANSCRIPT.search(_as_text(feat.get('note')))
    if m:
        return m.group(1)
    if feat.get('type') == bc.FEATURE_CDS:
        m = RE_PROTEIN_ID_TRANSCRIPT.match(_as_text(feat.get('protein_id')).strip())
        if m:
            return m.group(1)
    return None


RE_ID_TOKEN = re.compile(r'[A-Za-z0-9_.:\-]*\d[A-Za-z0-9_.:\-]*')


def mrna_aliases(feat: dict) -> set:
    """Every identifier an mRNA is known by: its transcript ID plus
    identifier-like tokens in its /note (RefSeq Drosophila mRNAs carry the
    RefSeq accession as /transcript_id but list the FlyBase transcript name,
    which the CDS refers to, only in /note: "PlexB-RB; ...; CG17245-RB")."""
    names = set(RE_ID_TOKEN.findall(_as_text(feat.get('note'))))
    tid = transcript_id_of(feat)
    if tid:
        names.add(tid)
    return names


def _parts(feat: dict) -> list:
    """Exon/CDS parts as ascending (start, stop) pairs."""
    starts, stops = feat.get('starts'), feat.get('stops')
    if isinstance(starts, list) and isinstance(stops, list) and starts and len(starts) == len(stops):
        return sorted(zip(starts, stops))
    return [(feat['start'], feat['stop'])]


def _junctions(parts: list) -> set:
    """Introns as (end of left part, start of right part), genomic order."""
    return {(parts[k][1], parts[k + 1][0]) for k in range(len(parts) - 1)}


def cds_fits_mrna(cds: dict, mrna: dict) -> bool:
    """The CDS lies within the mRNA and their splice structure agrees inside
    the CDS: every CDS intron is an mRNA intron, and the mRNA has no intron
    within the CDS span that the CDS lacks."""
    if not (mrna['start'] <= cds['start'] and cds['stop'] <= mrna['stop']):
        return False
    cds_j = _junctions(_parts(cds))
    mrna_j = {j for j in _junctions(_parts(mrna)) if cds['start'] < j[0] and j[1] < cds['stop']}
    return cds_j == mrna_j


def pair_transcripts(features: Sequence[dict]) -> Tuple[Dict[int, dict], Dict[int, dict], set, Dict[int, List[dict]]]:
    """Pair CDS with mRNA features.

    Returns ``(cds_of_mrna, mrna_of_cds, confident, candidates)`` keyed by
    ``id(feature)``. ``confident`` holds the ``id`` of mRNAs paired by
    transcript identifier or by a unique structural match; ``candidates``
    maps each mRNA to every structurally compatible CDS.
    """
    groups: Dict[tuple, Dict[str, list]] = {}
    for feat in features:
        if feat.get('type') in (bc.FEATURE_CDS, bc.FEATURE_MRNA):
            key = (feat.get('locus'), feat.get('sequence'), feat.get('strand'))
            slot = 'cds' if feat['type'] == bc.FEATURE_CDS else 'mrna'
            groups.setdefault(key, {'cds': [], 'mrna': []})[slot].append(feat)

    cds_of_mrna: Dict[int, dict] = {}
    mrna_of_cds: Dict[int, dict] = {}
    confident: set = set()
    candidates: Dict[int, List[dict]] = {}

    def pair(m, c, sure):
        cds_of_mrna[id(m)] = c
        mrna_of_cds[id(c)] = m
        if sure:
            confident.add(id(m))

    for group in groups.values():
        mrnas, cdss = group['mrna'], group['cds']
        for m in mrnas:
            candidates[id(m)] = [c for c in cdss if cds_fits_mrna(c, m)]

        # 1. identical transcript identifiers
        cds_by_tid: Dict[str, list] = {}
        for c in cdss:
            tid = transcript_id_of(c)
            if tid:
                cds_by_tid.setdefault(tid, []).append(c)
        for m in mrnas:
            matches = [c for alias in mrna_aliases(m) for c in cds_by_tid.get(alias, [])]
            unique_matches = {id(c): c for c in matches}
            if len(unique_matches) == 1:
                c = next(iter(unique_matches.values()))
                if id(c) not in mrna_of_cds:
                    pair(m, c, True)

        # 2. unique structural matches, repeated as pairs are removed
        progress = True
        while progress:
            progress = False
            for m in mrnas:
                if id(m) in cds_of_mrna:
                    continue
                free = [c for c in candidates[id(m)] if id(c) not in mrna_of_cds]
                if len(free) == 1:
                    pair(m, free[0], True)
                    progress = True

        # 3. still ambiguous: pair in feature order for the hierarchy only
        for m in mrnas:
            if id(m) in cds_of_mrna:
                continue
            free = [c for c in candidates[id(m)] if id(c) not in mrna_of_cds]
            if free:
                pair(m, free[0], False)

    return cds_of_mrna, mrna_of_cds, confident, candidates


def gff_ids(features: Sequence[dict]) -> Tuple[Dict[int, str], Dict[int, Tuple[str, str]], Dict[int, Optional[str]], Dict[int, str]]:
    """GFF3 identifiers for eukaryotic mRNA and CDS features.

    Returns ``(mrna_ids, cds_ids, parents, gene_ids)``; ``parents[id(transcript)]`` is
    the GFF3 ID of the transcript's gene, or None if the input has no gene for
    it. A transcript belongs to the gene whose locus tag or /gene equals its
    locus (Ensembl links by /gene). Also: ``mrna_ids[id(mRNA)]`` is the transcript
    ID, and ``cds_ids[id(CDS)]`` is ``(cds_id, parent_id)``; ``parent_id`` is None when there is neither an mRNA nor a gene to hang the CDS on. Transcripts use
    the input's own identifier where it has one, else ``<locus>-T<n>``
    numbered per locus in feature order. Every ID is unique in the file, and
    all segments of one CDS share one ID (GFF3: a spliced CDS is a single
    feature written over several lines).
    """
    _, mrna_of_cds, _, _ = pair_transcripts(features)
    used: set = set()

    def unique(name: str) -> str:
        if name not in used:
            used.add(name)
            return name
        k = 2
        while f"{name}_{k}" in used:
            k += 1
        used.add(f"{name}_{k}")
        return f"{name}_{k}"

    # gene IDs are their locus tags (see gff.write_gene_feature); a gene can
    # be referred to by that or by its /gene value
    # (inputs occasionally give two different genes the same /gene, so gene
    # IDs are made unique too and a child is attached to the gene containing it)
    genes_of_key: Dict[str, list] = {}
    gene_ids: Dict[int, str] = {}
    for f in features:
        if f.get('type') == bc.FEATURE_GENE:
            base = f.get('locus') or f"{f.get('sequence')}_{f['start']}_{f['stop']}_{f.get('strand')}"
            gene_ids[id(f)] = unique(base)
            for key in {f.get('locus'), f.get('gene')} - {None, ''}:
                genes_of_key.setdefault(key, []).append(f)

    def gene_for(child: dict) -> Optional[str]:
        cands = genes_of_key.get(child.get('locus'), [])
        inside = [g for g in cands if g.get('sequence') == child.get('sequence')
                  and g['start'] <= child['start'] and child['stop'] <= g['stop']]
        pick = (inside or cands or [None])[0]
        return gene_ids[id(pick)] if pick is not None else None

    cds_by_id = {id(f): f for f in features if f.get('type') == bc.FEATURE_CDS}
    paired_cds_of_mrna = {id(m): cds_by_id[cid] for cid, m in mrna_of_cds.items()}

    mrna_ids: Dict[int, str] = {}
    parents: Dict[int, Optional[str]] = {}
    per_locus: Dict[str, int] = {}
    for feat in features:
        # transcripts: mRNAs, and eukaryotic tRNA/rRNA/ncRNA (which share the
        # per-locus -T<n> numbering so their IDs cannot collide with each
        # other or with the gene line, which uses the bare locus tag)
        if feat.get('type') not in (bc.FEATURE_MRNA, bc.FEATURE_T_RNA, bc.FEATURE_R_RNA, bc.FEATURE_NC_RNA):
            continue
        locus = feat.get('locus') or f"{feat.get('sequence')}_{feat['start']}_{feat['stop']}_{feat.get('strand')}"
        tid = transcript_id_of(feat)
        if not tid and id(feat) in paired_cds_of_mrna:
            tid = transcript_id_of(paired_cds_of_mrna[id(feat)])
        if not tid:
            per_locus[locus] = per_locus.get(locus, 0) + 1
            tid = f"{locus}-T{per_locus[locus]}"
        mrna_ids[id(feat)] = unique(tid)
        parents[id(feat)] = gene_for(feat)

    cds_ids: Dict[int, Tuple[str, str]] = {}
    for feat in features:
        if feat.get('type') != bc.FEATURE_CDS:
            continue
        locus = feat.get('locus') or f"{feat.get('sequence')}_{feat['start']}_{feat['stop']}_{feat.get('strand')}"
        mrna = mrna_of_cds.get(id(feat))
        if mrna is not None:
            parent = mrna_ids[id(mrna)]
            cds_ids[id(feat)] = (unique(f"{parent}.cds"), parent)
        else:
            # no mRNA in the input: hang the CDS off its gene, if there is one
            base = transcript_id_of(feat) or locus
            cds_ids[id(feat)] = (unique(f"{base}.cds"), gene_for(feat))
    return mrna_ids, cds_ids, parents, gene_ids
