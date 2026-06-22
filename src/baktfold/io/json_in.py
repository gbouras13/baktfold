import json
import sys
from datetime import datetime
from pathlib import Path

from loguru import logger
from xopen import xopen

import baktfold.bakta.constants as bc
import baktfold.bakta.config as cfg
# import baktfold.utils as bu
# import baktfold.io.fasta as fasta
# import baktfold.io.tsv as tsv
# import baktfold.io.gff as gff
# import baktfold.io.insdc as insdc
# import baktfold.plot as plot


def parse_json_input(input_path, faa_path, all_proteins, protein_json_flag):
    """
    Parses genome annotations from input JSON file.

    Args:
      input_path (str): Path to input JSON file.
      faa_path (str): Path to output file for hypothetical proteins.
      all_proteins (bool): Whether to keep all proteins or only hypothetical ones.
      protein_json_flag (bool): Whether input is protein JSON

    Returns:
      tuple: A tuple containing the data, features, and whether there are duplicate locus tags.

    Examples:
      >>> parse_json_input('input.json', 'hypotheticals.faa', False, False)
      (data, features, False, False)
    """

    

    ############################################################################
    # Checks and configurations
    # - check parameters and setup global configuration
    # - test database
    # - test binary dependencies
    ############################################################################

    try:
        if input_path == '':
            raise ValueError('File path argument must be non-empty')
        annotation_path = Path(input_path).resolve()
        cfg.check_readability('annotation', annotation_path)
        cfg.check_content_size('annotation', annotation_path)
    except Exception as e:
        logger.error(f'ERROR: annotation file {annotation_path} not valid! {e}')
        sys.exit(1)
    
    #print(f'baktfold v{cfg.version}')

    logger.info(f'Parsing annotations from input: {annotation_path}')
    with xopen(str(annotation_path), threads=0) as fh:
        data = json.load(fh)


    features = data['features']

    # features_by_sequence = {seq['id']: [] for seq in data['sequences']}
    # for feature in data['features']:
    #     seq_id = feature['sequence'] if 'sequence' in feature else feature['contig']  # <1.10.0 compatibility
    #     sequence_features = features_by_sequence.get(seq_id)
    #     sequence_features.append(feature)

    # keep all proteins
    if all_proteins:
        hypotheticals = [feat for feat in features if feat['type'] == bc.FEATURE_CDS ]
    else:
        hypotheticals = [feat for feat in features if feat['type'] == bc.FEATURE_CDS and 'hypothetical' in feat]


    if protein_json_flag: # this will also be only hypotheticals if protein mode (or else why not just run with the FASTA)
        version = data.get("version", {})
        return features, hypotheticals, version


    # check if dupe locus tags (euks can have multiple CDS same locus tag e.g. Cladocopium goreaui CAMXCT020000001.1)
    seen_loci = set()
    has_duplicate_locus = False

    for feat in hypotheticals:
        locus = feat['locus']
        if locus in seen_loci:
            has_duplicate_locus = True
            logger.warning("Multiple CDS per locus tag were detected in your input JSON.")
            logger.warning("CDS id (which is unique) rather than locus tag will be used for ProstT5+Foldseek searches.")
            break
        seen_loci.add(locus)

    # this is done after getting all the sequences into the dict for baktfold proteins

    if has_duplicate_locus:
        # write hypothetical proteins to file with id (not locus) as guaranteed exists and unique
        with faa_path.open('wt') as fh:
            for feat in hypotheticals:
                fh.write(f">{feat['id']}\n{feat['aa']}\n")

    else:
        # write hypothetical proteins to file - almost always
        with faa_path.open('wt') as fh:
            for feat in hypotheticals:
                fh.write(f">{feat['locus']}\n{feat['aa']}\n")

    # none of this is relevant for proteins
    try:
        genome_block = data.get("genome")

        if genome_block is None:
            logger.error("No 'genome' block found in input JSON. Please check.")
            translation_table = None
        else:
            if "translation_table" not in genome_block:
                logger.error("No translation table found in input JSON. Please check your input.")
            else:
                raw_value = genome_block["translation_table"]

                try:
                    translation_table = int(raw_value)
                    logger.info(
                        f"Translation table {translation_table} detected from input JSON"
                    )

                except (ValueError, TypeError):
                    translation_table = str(raw_value)
                    logger.warning(
                        f"Translation table '{raw_value}' is not an integer. "
                        f"Parsing it as a string."
                    )

    except Exception as e:
        logger.exception(
            f"Unexpected error while parsing translation table: {e}"
        )
        translation_table = None

    # input detection

    version = data.get("version", {})

    prokka = False
    other_genbank = False

    if "prokka" in version:
        prokka = True
        logger.info("Prokka input detected")
    if  "prokka"  not in version and "bakta" not in version:
        other_genbank = True

    logger.info('Parsing complete')

    return data, features, has_duplicate_locus, translation_table, prokka, other_genbank, version


def log_for_other_genbank_tools(cds_program,trna_program, rrna_program, tmrna_program, ncrna_program):

    logger.warning("Neither bakta nor prokka input detected")
    logger.info("If you would like to specify consituent inference tools for CDS, tRNA, rRNA, tmRNA and ncRNA")
    logger.info("Reminder: please use --cds-tool --trna-program --rrna-program  --tmrna-program --ncrna-program to modify them if you haven't already")
    logger.info(f"For this genome, they are --cds_program {cds_program} --trna-program {trna_program} --rrna_program {rrna_program} --tmrna-program {tmrna_program} --ncrna-program  {ncrna_program}")


############################################################################
# Reconstruction of all (non-Foldseek) outputs from a baktfold/bakta JSON
# - used by the ``baktfold json`` subcommand
############################################################################

# Default constituent-tool inference strings. Mirror the ``compare`` CLI
# defaults so reconstruction is identical for non-Bakta/non-Prokka input when
# the JSON carries no provenance block and the user passes no overrides.
_DEFAULT_PROGRAMS = {
    'cds_program': 'Prodigal:2.6',
    'trna_program': 'tRNAscan-SE:2.0.12',
    'rrna_program': 'INFERNAL:1.1.5',
    'tmrna_program': 'INFERNAL:1.1.5',
    'ncrna_program': 'INFERNAL:1.1.5',
}

# Feature types that only ever appear in eukaryotic annotations.
_EUK_FEATURE_TYPES = frozenset({
    bc.FEATURE_GENE, bc.FEATURE_MRNA, bc.FEATURE_5UTR, bc.FEATURE_3UTR, bc.FEATURE_REPEAT,
})


def _first_set(*values):
    """Return the first value that is not None (used for override > provenance > inference)."""
    for value in values:
        if value is not None:
            return value
    return None


def _infer_euk(features) -> bool:
    """Eukaryotic input is unambiguous from feature types (gene/mRNA/UTR/repeat_region)."""
    return any(feat.get('type') in _EUK_FEATURE_TYPES for feat in features)


def _infer_custom_db(features) -> bool:
    """A custom-DB hit leaves a ``custom:`` db_xref on at least one CDS."""
    for feat in features:
        for xref in feat.get('db_xrefs', []):
            if isinstance(xref, str) and xref.startswith('custom:'):
                return True
    return False


def _detect_duplicate_locus(features) -> bool:
    """Reproduce parse_json_input's duplicate-locus scan over hypothetical CDS."""
    seen = set()
    for feat in features:
        if feat.get('type') != bc.FEATURE_CDS or 'hypothetical' not in feat:
            continue
        locus = feat.get('locus')
        if locus in seen:
            return True
        seen.add(locus)
    return False


def _coerce_translation_table(raw, default: int = 11):
    if raw is None:
        return default
    try:
        return int(raw)
    except (ValueError, TypeError):
        logger.warning(f"Translation table '{raw}' is not an integer; using it as-is.")
        return raw


def _restore_run_timing(data: dict) -> None:
    """Best-effort restoration of cfg.run_start/run_end from the JSON 'run' block."""
    run_block = data.get('run') or {}
    for key, attr in (('start', 'run_start'), ('end', 'run_end')):
        value = run_block.get(key)
        if not value:
            continue
        try:
            setattr(cfg, attr, datetime.strptime(value, '%Y-%m-%d %H:%M:%S'))
        except (ValueError, TypeError):
            pass  # leave the module default in place


def parse_baktfold_json_for_reconstruction(
    input_path,
    euk_override=None,
    custom_db_override=None,
    fast_override=None,
    program_overrides=None,
):
    """Parse a baktfold (or bakta) JSON and reconstitute every argument needed to
    re-run the output writers (genome -> write_bakta_outputs, proteins ->
    write_bakta_proteins_outputs). No database or Foldseek run is required.

    Resolution order for the runtime flags that are not otherwise recoverable
    from feature data is: CLI override > JSON ``baktfold_run`` provenance block >
    inference from features > hard default.

    Returns:
        dict: keyed by 'mode' ('genome' | 'proteins') plus everything the
        matching writer needs.
    """
    program_overrides = program_overrides or {}

    # ---- validate & load -------------------------------------------------
    try:
        if input_path == '':
            raise ValueError('File path argument must be non-empty')
        annotation_path = Path(input_path).resolve()
        cfg.check_readability('annotation', annotation_path)
        cfg.check_content_size('annotation', annotation_path)
    except Exception as e:
        logger.error(f'ERROR: annotation file {input_path} not valid! {e}')
        sys.exit(1)

    logger.info(f'Parsing baktfold JSON for reconstruction: {annotation_path}')
    with xopen(str(annotation_path), threads=0) as fh:
        data = json.load(fh)

    if 'features' not in data:
        logger.error("Input JSON has no 'features' key - not a valid baktfold/bakta JSON.")
        sys.exit(1)

    features = data['features']
    version = data.get('version', {})
    provenance = data.get('baktfold_run', {})

    # ---- restore faithful cfg state (headers, version strings, timing) ---
    cfg.version = version.get('baktfold', cfg.version)
    cfg.db_version = version.get('baktfold_db', cfg.db_version)
    _restore_run_timing(data)

    # ---- detect genome vs proteins mode ----------------------------------
    mode = provenance.get('mode')
    if mode not in ('genome', 'proteins'):
        mode = 'genome' if 'sequences' in data else 'proteins'

    # ---- flags resolvable for both modes ---------------------------------
    custom_db = _first_set(custom_db_override, provenance.get('custom_db'), _infer_custom_db(features))
    if custom_db is None:
        custom_db = False
    # ``fast`` is NOT reliably inferable (absence of AFDB hits is ambiguous):
    # rely on provenance, else default off.
    fast = _first_set(fast_override, provenance.get('fast'))
    if fast is None:
        fast = False

    if mode == 'proteins':
        logger.info('Proteins-mode JSON detected.')
        return {
            'mode': 'proteins',
            'data': data,
            'aas': features,
            'features': features,
            'custom_db': bool(custom_db),
            'fast': bool(fast),
            'bakta_version': version,
        }

    # ---- genome mode -----------------------------------------------------
    if 'sequences' not in data:
        logger.error("Genome-mode JSON is missing the 'sequences' block; cannot reconstruct.")
        sys.exit(1)

    prokka = provenance.get('prokka')
    if prokka is None:
        prokka = 'prokka' in version
    other_genbank = provenance.get('other_genbank')
    if other_genbank is None:
        other_genbank = ('prokka' not in version) and ('bakta' not in version)

    euk = _first_set(euk_override, provenance.get('euk'), _infer_euk(features))
    if euk is None:
        euk = False

    translation_table = provenance.get('translation_table')
    if translation_table is None:
        translation_table = _coerce_translation_table((data.get('genome') or {}).get('translation_table'))

    has_duplicate_locus = provenance.get('has_duplicate_locus')
    if has_duplicate_locus is None:
        has_duplicate_locus = _detect_duplicate_locus(features)

    programs = {}
    for name, default in _DEFAULT_PROGRAMS.items():
        programs[name] = _first_set(program_overrides.get(name), provenance.get(name), default)

    # ---- rebuild features_by_sequence (preserving the JSON feature order,
    #      which is already start-sorted per sequence) and a flattened list
    #      consistent with it. Mirrors run/compare (skip discarded, <1.10.0
    #      'contig' fallback). -------------------------------------------------
    features_by_sequence = {seq['id']: [] for seq in data['sequences']}
    for feature in features:
        if 'discarded' in feature:
            continue
        seq_id = feature['sequence'] if 'sequence' in feature else feature.get('contig')  # <1.10.0 compat
        if seq_id is None:
            logger.warning(f"Feature missing 'sequence', skipping: id={feature.get('id')}")
            continue
        bucket = features_by_sequence.get(seq_id)
        if bucket is None:
            logger.warning(f"Feature references unknown sequence '{seq_id}', skipping: id={feature.get('id')}")
            continue
        bucket.append(feature)

    flattened = []
    for seq in data['sequences']:
        flattened.extend(features_by_sequence[seq['id']])

    if other_genbank:
        log_for_other_genbank_tools(
            programs['cds_program'], programs['trna_program'], programs['rrna_program'],
            programs['tmrna_program'], programs['ncrna_program'],
        )

    logger.info(
        f"Reconstruction settings: euk={euk}, custom_db={bool(custom_db)}, fast={bool(fast)}, "
        f"prokka={prokka}, other_genbank={other_genbank}, translation_table={translation_table}, "
        f"has_duplicate_locus={has_duplicate_locus}"
    )

    return {
        'mode': 'genome',
        'data': data,
        'features': flattened,
        'features_by_sequence': features_by_sequence,
        'has_duplicate_locus': bool(has_duplicate_locus),
        'translation_table': translation_table,
        'prokka': bool(prokka),
        'other_genbank': bool(other_genbank),
        'euk': bool(euk),
        'custom_db': bool(custom_db),
        'fast': bool(fast),
        'cds_program': programs['cds_program'],
        'trna_program': programs['trna_program'],
        'rrna_program': programs['rrna_program'],
        'tmrna_program': programs['tmrna_program'],
        'ncrna_program': programs['ncrna_program'],
        'bakta_version': version,
    }


    