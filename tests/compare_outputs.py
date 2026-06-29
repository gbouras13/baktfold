#!/usr/bin/env python3
"""
Compare baktfold outputs between two runs (dev pholdlib-refactored vs bioconda reference).
Ignores timestamp-dependent content (log files, dates in annotation JSON run key).

Usage:
    python tests/compare_outputs.py <dir_dev> <dir_ref> [--cpu]

Exit code 0 = identical (modulo timestamps), non-zero = differences found.
"""
import json
import math
import re
import sys
from pathlib import Path

# ── patterns that are timestamp/run-specific and should be ignored ──────────
SKIP_LINE_PATTERNS = [
    re.compile(r"^\d{4}-\d{2}-\d{2}"),          # log lines: 2026-05-26 ...
    re.compile(r"#.*baktfold.*run"),              # any baktfold run pragma
]

# file extensions to skip entirely
SKIP_EXTENSIONS = {".log"}

# files to skip by name
SKIP_FILENAMES = set()

# directory components to skip entirely (any file under these dirs is ignored)
SKIP_DIRS = {"logs", "logdir"}


# Foldseek hit accessions are non-deterministic: for near-tied structural hits
# (same function, different accession) GPU score wobble flips which accession
# wins the tophit, and that flip propagates into the db_xref / hit columns of
# every output (gff/gbff/embl/tsv/inference.tsv/json). In tolerant (non-strict)
# mode the *category* of hit (swissprot/afdbclusters/pdb/cath/custom) is kept but
# the specific accession is normalised away, so the comparison still flags a
# changed product, a missing DB category, a coordinate change, etc. — just not
# which of two equivalent accessions won the tie. Bakta's own db_xrefs (SO:,
# UniRef:, ...) are deterministic and are NOT touched.
_FOLDSEEK_ACCESSION_RE = re.compile(
    r"\b(swissprot|afdbclusters|pdb|cath|custom)_[A-Za-z0-9._-]+"
)


def _normalize_accessions(text: str) -> str:
    return _FOLDSEEK_ACCESSION_RE.sub(r"\1_<X>", text)


# Columns of the per-protein TSV (inference.tsv / proteins .tsv) derived from
# the non-deterministic ProstT5/Foldseek numerics — dropped before comparison.
_CONFIDENCE_COLUMNS = {"Annotation_Confidence", "TMscore", "LDDT"}


def _drop_confidence_columns(lines: list) -> list:
    """Header-aware drop of the confidence/TM-score/LDDT columns from a TSV.

    Finds the header row (first non-``#`` line); if it carries any of
    ``_CONFIDENCE_COLUMNS``, drops those columns from the header and every data
    row. Files without those columns (e.g. the human-readable feature TSV, or a
    reference produced before the feature existed) are returned unchanged, so
    the two sides still align.
    """
    header_idx = next((i for i, l in enumerate(lines) if not l.startswith("#")), None)
    if header_idx is None:
        return lines
    header = lines[header_idx].split("\t")
    drop = {i for i, col in enumerate(header) if col in _CONFIDENCE_COLUMNS}
    if not drop:
        return lines
    out = []
    for i, line in enumerate(lines):
        if i < header_idx or line.startswith("#"):
            out.append(line)
        else:
            out.append("\t".join(p for j, p in enumerate(line.split("\t")) if j not in drop))
    return out


def filter_lines(path: Path) -> list:
    """Read a file and return lines with timestamp-like content removed."""
    try:
        lines = path.read_text(errors="replace").splitlines()
    except Exception as e:
        return [f"<ERROR reading {path}: {e}>"]
    return [l for l in lines if not any(p.search(l) for p in SKIP_LINE_PATTERNS)]


def _csv_float_differ(lines_dev: list, lines_ref: list, tol: float = 0.01) -> list:
    """Return diff messages for two sorted lists of CSV lines where the second
    column is a float (e.g. mean_probabilities.csv: 'seq_id,mean_prob').
    Lines with matching seq_ids are compared numerically; count mismatches
    are flagged."""
    row_diffs = []
    if len(lines_dev) != len(lines_ref):
        row_diffs.append(f"    line count: dev={len(lines_dev)} ref={len(lines_ref)}")
    for i, (a, b) in enumerate(zip(lines_dev, lines_ref)):
        if a == b:
            continue
        pa, pb = a.split(",", 1), b.split(",", 1)
        if pa[0] != pb[0]:
            row_diffs.append(f"    seq_id mismatch dev[{i}]: {a[:120]} | ref: {b[:120]}")
            continue
        try:
            if not math.isclose(float(pa[1]), float(pb[1]), abs_tol=tol):
                row_diffs.append(f"    value mismatch dev[{i}]: {a[:120]}")
                row_diffs.append(f"                   ref[{i}]: {b[:120]}")
        except ValueError:
            row_diffs.append(f"    parse error dev[{i}]: {a[:120]}")
    return row_diffs


def _jsonl_float_differ(lines_dev: list, lines_ref: list, tol: float = 0.01) -> list:
    """Return diff messages for two sorted lists of JSONL lines.

    Each line is a JSON object {"seq_id": str, "probability": [float, ...]}.
    Float values in the probability list are compared with tolerance *tol*
    to absorb any floating-point differences between torch versions.
    The seq_id is compared exactly; the two lists must be in the same order.
    """
    row_diffs = []
    if len(lines_dev) != len(lines_ref):
        row_diffs.append(f"    line count: dev={len(lines_dev)} ref={len(lines_ref)}")
    for i, (a, b) in enumerate(zip(lines_dev, lines_ref)):
        if a == b:
            continue
        try:
            da, db = json.loads(a), json.loads(b)
        except Exception:
            row_diffs.append(f"    JSON parse error dev[{i}]: {a[:120]}")
            continue
        if da.get("seq_id") != db.get("seq_id"):
            row_diffs.append(
                f"    seq_id mismatch dev[{i}]: {da.get('seq_id')} | ref: {db.get('seq_id')}"
            )
            continue
        pa, pb = da.get("probability", []), db.get("probability", [])
        if len(pa) != len(pb):
            row_diffs.append(
                f"    prob length mismatch for {da.get('seq_id')}: "
                f"dev={len(pa)} ref={len(pb)}"
            )
            continue
        bad_pos = [
            j for j, (fa, fb) in enumerate(zip(pa, pb))
            if not math.isclose(fa, fb, abs_tol=tol)
        ]
        if bad_pos:
            row_diffs.append(
                f"    prob mismatch for {da.get('seq_id')} at positions {bad_pos[:10]}: "
                f"dev={[pa[j] for j in bad_pos[:3]]} ref={[pb[j] for j in bad_pos[:3]]}"
            )
    return row_diffs


def _parse_fasta(lines: list) -> dict:
    """Parse FASTA lines into an ordered {id: sequence} dict (id = first token)."""
    records = {}
    current = None
    for line in lines:
        if line.startswith(">"):
            current = line[1:].split()[0]
            records[current] = []
        elif current is not None:
            records[current].append(line.strip())
    return {k: "".join(v) for k, v in records.items()}


def _fasta_3di_differ(lines_dev: list, lines_ref: list, min_identity: float) -> list:
    """Compare two 3Di FASTA files tolerantly.

    ProstT5 3Di predictions are not bit-identical across hardware / driver /
    torch versions, so each sequence is compared by per-residue identity and
    only flagged when it falls below *min_identity*. Missing/extra sequences and
    length mismatches are always reported (those indicate a real change, since
    the 3Di length equals the protein length).
    """
    dev = _parse_fasta(lines_dev)
    ref = _parse_fasta(lines_ref)
    diffs = []

    only_dev = sorted(set(dev) - set(ref))
    only_ref = sorted(set(ref) - set(dev))
    if only_dev:
        diffs.append(f"    seq ids only in dev: {only_dev[:10]}")
    if only_ref:
        diffs.append(f"    seq ids only in ref: {only_ref[:10]}")

    low = []
    for sid in sorted(set(dev) & set(ref)):
        a, b = dev[sid], ref[sid]
        if len(a) != len(b):
            diffs.append(f"    length mismatch {sid}: dev={len(a)} ref={len(b)}")
            continue
        if not a:
            continue
        identity = sum(1 for x, y in zip(a, b) if x == y) / len(a)
        if identity < min_identity:
            low.append((sid, identity))

    if low:
        low.sort(key=lambda t: t[1])
        worst = ", ".join(f"{sid}={ident:.0%}" for sid, ident in low[:10])
        diffs.append(f"    {len(low)} seq(s) below {min_identity:.0%} 3Di identity (worst: {worst})")
    return diffs


def _tophit_differ(lines_dev: list, lines_ref: list, strict: bool = False) -> list:
    """Compare two Foldseek *_tophit.tsv files.

    Columns: query target bitscore fident evalue qStart qEnd qLen qCov tStart
    tEnd tLen tCov. Because the 3Di input is itself non-deterministic, every
    per-alignment number wobbles run-to-run on GPU (scores, identity, alignment
    extent, coverage), and the *target accession* itself flips between near-tied
    structural homologs. So in tolerant (non-strict) mode only the deterministic
    facts are kept — which query proteins got a hit, and their length — as a
    deduplicated set; a protein gaining/losing a hit in a DB is still caught.
    In strict mode (run_comparison, deterministic CPU) only the three wobbly
    quality scores are dropped.
    """
    if strict:
        keep, dedupe = (0, 1, 5, 6, 7, 8, 9, 10, 11, 12), False  # drop bitscore/fident/evalue
    else:
        keep, dedupe = (0, 7), True  # query, qLen

    def _scrub(lines):
        out = []
        for line in lines:
            parts = line.split("\t")
            if len(parts) >= 13:
                parts = [parts[i] for i in keep]
            out.append("\t".join(parts))
        return sorted(set(out)) if dedupe else sorted(out)

    sd, sr = _scrub(lines_dev), _scrub(lines_ref)
    diffs = []
    if sd != sr:
        for i, (a, b) in enumerate(zip(sd, sr)):
            if a != b:
                diffs.append(f"    dev[{i}]: {a[:140]}")
                diffs.append(f"    ref[{i}]: {b[:140]}")
                if i > 10:
                    diffs.append("    ... (truncated)")
                    break
        if len(sd) != len(sr):
            diffs.append(f"    line count: dev={len(sd)} ref={len(sr)}")
    return diffs


# Keys dropped from the annotation JSON before comparison because they are not
# reproducible run-to-run on GPU. Hit *identity* is still checked (via the
# feature 'db_xrefs' and the pstc 'source'/'description'); only the raw ProstT5
# 3Di string and the Foldseek alignment numerics are dropped.
_VOLATILE_JSON_KEYS = (
    "id",            # feature ids carry a random 2-char suffix (locus is stable)
    "3di",           # raw ProstT5 3Di prediction (checked tolerantly in _3di.fasta)
    "score",         # Foldseek bitscore
    "evalue",        # Foldseek e-value
    "query_cov",     # Foldseek query coverage
    "subject_cov",   # Foldseek subject coverage
    "identity",      # Foldseek fraction-identity
    "prostt5_confidence",     # mean ProstT5 confidence (wobbles run-to-run)
    "annotation_confidence",  # derived from the above non-deterministic metrics
    "tmscore",       # Foldseek TM-score (structure input)
    "lddt",          # Foldseek LDDT (structure input)
)


def _strip_volatile_fields(obj) -> None:
    """Recursively strip non-reproducible fields from bakta annotation objects.

    See ``_VOLATILE_JSON_KEYS``. ProstT5 3Di prediction and Foldseek alignment
    scores are not bit-identical across runs/hardware, so they are removed here
    so the comparison only flags real annotation differences.
    """
    if isinstance(obj, list):
        for item in obj:
            _strip_volatile_fields(item)
    elif isinstance(obj, dict):
        for key in _VOLATILE_JSON_KEYS:
            obj.pop(key, None)
        for v in obj.values():
            _strip_volatile_fields(v)


def _compare_annotation_json(fd: Path, fr: Path, normalize: bool = False) -> list:
    """Compare two bakta annotation JSON files, ignoring the 'run' key (timestamps),
    randomly-generated feature 'id' values and the non-reproducible ProstT5/Foldseek
    fields. With *normalize*, Foldseek hit accessions are also normalised away.

    Returns a list of diff messages (empty = identical modulo timestamps).
    """
    diffs = []
    try:
        dev_obj = json.loads(fd.read_text(errors="replace"))
        ref_obj = json.loads(fr.read_text(errors="replace"))
    except json.JSONDecodeError as e:
        diffs.append(f"    JSON parse error: {e}")
        return diffs

    # strip run-specific keys
    for obj in (dev_obj, ref_obj):
        obj.pop("run", None)
        _strip_volatile_fields(obj)

    dev_str = json.dumps(dev_obj, sort_keys=True, separators=(",", ":"))
    ref_str = json.dumps(ref_obj, sort_keys=True, separators=(",", ":"))
    if normalize:
        dev_str = _normalize_accessions(dev_str)
        ref_str = _normalize_accessions(ref_str)
    if dev_str != ref_str:
        for i, (a, b) in enumerate(zip(dev_str, ref_str)):
            if a != b:
                start = max(0, i - 50)
                diffs.append(f"    first diff at char {i}:")
                diffs.append(f"      dev: ...{dev_str[start:i+40]!r}...")
                diffs.append(f"      ref: ...{ref_str[start:i+40]!r}...")
                break
        if len(dev_str) != len(ref_str):
            diffs.append(f"    length: dev={len(dev_str)} ref={len(ref_str)}")
    return diffs


def compare_dirs(dir_dev: Path, dir_ref: Path, strict: bool = False) -> list:
    """Recursively compare two directories. Returns list of diff messages."""
    diffs = []

    dev_files = {f.relative_to(dir_dev) for f in dir_dev.rglob("*") if f.is_file()}
    ref_files = {f.relative_to(dir_ref) for f in dir_ref.rglob("*") if f.is_file()}

    def should_skip(rel: Path) -> bool:
        return (
            rel.suffix in SKIP_EXTENSIONS
            or rel.name in SKIP_FILENAMES
            or bool(SKIP_DIRS.intersection(rel.parts))
        )

    for f in sorted(dev_files - ref_files):
        if not should_skip(f):
            diffs.append(f"  ONLY IN DEV : {f}")

    for f in sorted(ref_files - dev_files):
        if not should_skip(f):
            diffs.append(f"  ONLY IN REF : {f}")

    for rel in sorted(dev_files & ref_files):
        if should_skip(rel):
            continue

        fd = dir_dev / rel
        fr = dir_ref / rel

        ld = filter_lines(fd)
        lr = filter_lines(fr)
        if not strict:  # ignore non-deterministic Foldseek hit accessions
            ld = [_normalize_accessions(line) for line in ld]
            lr = [_normalize_accessions(line) for line in lr]

        # ── mean_probabilities.csv ─────────────────────────────────────────
        if "mean_probabilities" in rel.name and rel.suffix == ".csv":
            prob_tol = 0.01 if strict else 0.5
            sd, sr = sorted(ld), sorted(lr)
            row_diffs = _csv_float_differ(sd, sr, tol=prob_tol)
            if row_diffs:
                diffs.append(f"  DIFFER (sorted, tol={prob_tol}) : {rel}")
                diffs.extend(row_diffs[:22])

        # ── all_probabilities.json (JSONL probability file) ────────────────
        elif rel.suffix == ".json" and "_all_probabilities" in rel.name:
            json_tol = 0.01 if strict else 0.5
            def _sort_key(line):
                try:
                    return json.loads(line).get("seq_id", line)
                except Exception:
                    return line
            sd = sorted(ld, key=_sort_key)
            sr = sorted(lr, key=_sort_key)
            row_diffs = _jsonl_float_differ(sd, sr, tol=json_tol)
            if row_diffs:
                diffs.append(f"  DIFFER (sorted by seq_id, tol={json_tol}) : {rel}")
                diffs.extend(row_diffs[:22])
                if len(sd) != len(sr):
                    diffs.append(f"    line count: dev={len(sd)} ref={len(sr)}")

        # ── bakta annotation JSON (single large JSON, not JSONL) ───────────
        # Strip run.start/run.end timestamps before comparing.
        elif rel.suffix == ".json":
            ann_diffs = _compare_annotation_json(fd, fr, normalize=not strict)
            if ann_diffs:
                diffs.append(f"  DIFFER (annotation JSON, 'run' key stripped) : {rel}")
                diffs.extend(ann_diffs[:22])

        # ── Foldseek tophit TSV (non-deterministic accession/scores ignored) ─
        elif rel.suffix == ".tsv" and "_tophit" in rel.name:
            row_diffs = _tophit_differ(ld, lr, strict)
            if row_diffs:
                diffs.append(f"  DIFFER (tophit, hit set) : {rel}")
                diffs.extend(row_diffs[:22])

        # ── TSV/CSV/TXT (exact, sorted) ────────────────────────────────────
        elif rel.suffix in {".tsv", ".csv", ".txt"}:
            # drop the confidence/TM-score/LDDT columns (derived from the
            # non-deterministic ProstT5/Foldseek metrics) where present
            ld = _drop_confidence_columns(ld)
            lr = _drop_confidence_columns(lr)
            sd, sr = sorted(ld), sorted(lr)
            if sd != sr:
                diffs.append(f"  DIFFER (sorted) : {rel}")
                for i, (a, b) in enumerate(zip(sd, sr)):
                    if a != b:
                        diffs.append(f"    dev[{i}]: {a[:140]}")
                        diffs.append(f"    ref[{i}]: {b[:140]}")
                        if i > 10:
                            diffs.append("    ... (truncated)")
                            break
                if len(sd) != len(sr):
                    diffs.append(f"    line count: dev={len(sd)} ref={len(sr)}")

        # ── 3Di FASTA (ProstT5 prediction; per-residue identity tolerance) ──
        elif rel.suffix == ".fasta" and "_3di" in rel.name:
            min_identity = 1.0 if strict else 0.95
            row_diffs = _fasta_3di_differ(ld, lr, min_identity)
            if row_diffs:
                diffs.append(f"  DIFFER (3Di identity < {min_identity:.0%}) : {rel}")
                diffs.extend(row_diffs[:22])

        # ── FASTA and everything else (exact) ──────────────────────────────
        else:
            if ld != lr:
                diffs.append(f"  DIFFER : {rel}")
                for i, (a, b) in enumerate(zip(ld, lr)):
                    if a != b:
                        diffs.append(f"    dev[{i}]: {a[:140]}")
                        diffs.append(f"    ref[{i}]: {b[:140]}")
                        if i > 10:
                            diffs.append("    ... (truncated)")
                            break
                if len(ld) != len(lr):
                    diffs.append(f"    line count: dev={len(ld)} ref={len(lr)}")

    return diffs


def main():
    import argparse
    parser = argparse.ArgumentParser(
        description="Compare baktfold outputs between two runs."
    )
    parser.add_argument("dir_dev", type=Path, help="Dev output directory")
    parser.add_argument("dir_ref", type=Path, help="Ref output directory")
    parser.add_argument(
        "--cpu",
        action="store_true",
        help=(
            "Both runs used --cpu (deterministic). "
            "Probability CSVs compared with abs_tol=0.01; all other files exact."
        ),
    )
    args = parser.parse_args()

    dir_dev = args.dir_dev
    dir_ref = args.dir_ref
    strict = args.cpu

    for d, label in [(dir_dev, "dev"), (dir_ref, "ref")]:
        if not d.is_dir():
            print(f"ERROR: {label} directory does not exist: {d}")
            sys.exit(2)

    print(f"Comparing:\n  DEV: {dir_dev}\n  REF: {dir_ref}\n")
    if strict:
        print("Mode: --cpu (prob CSVs/JSON: abs_tol=0.01; all else exact)\n")
    else:
        print("Mode: MPS/GPU (prob CSVs/JSON: abs_tol=0.5)\n")

    diffs = compare_dirs(dir_dev, dir_ref, strict=strict)

    if diffs:
        print(f"DIFFERENCES FOUND ({len(diffs)} issues):")
        for d in diffs:
            print(d)
        sys.exit(1)
    else:
        print("ALL OUTPUTS MATCH (modulo timestamps).")
        sys.exit(0)


if __name__ == "__main__":
    main()
