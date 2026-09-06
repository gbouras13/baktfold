import json
from collections import OrderedDict
from pathlib import Path
from typing import Sequence

from loguru import logger

import baktfold.bakta
import baktfold.bakta.constants as bc
import baktfold.bakta.config as cfg


def write_json(data: dict, features: Sequence[dict], json_path: Path, bakta_version: dict, baktfold_run: dict = None):
    logger.info(f'write JSON: path={json_path}' )

    # clean feature attributes
    for feat in features:
        if(feat['type'] == bc.FEATURE_CDS or feat['type'] == bc.FEATURE_SORF):

            if isinstance(feat, dict) and 'aa_digest' in feat:
                feat.pop('aa_digest')  # remove binary aa digest before JSON serialization
            # remove redundant IPS Dbxrefs
            ips = feat.get('ips', None)
            if isinstance(ips, dict):
                ips.pop('db_xrefs', None)

            # remove redundant PSC Dbxrefs
            psc = feat.get('psc', None)
            if isinstance(psc, dict):
                psc.pop('db_xrefs', None)

    version = bakta_version
    version['baktfold'] = cfg.version
    version['baktfold_db'] = cfg.db_version

    # version['db'] = {
    #     'version': f"{cfg.db_info['major']}.{cfg.db_info['minor']}",
    #     'type': cfg.db_info['type']
    # }
    data['version'] = version

    # Persist a self-describing provenance block so that ``baktfold json`` can
    # later reconstitute every non-Foldseek output without the user having to
    # re-supply runtime flags (euk / custom_db / fast / inference tool strings).
    # These flags are NOT otherwise recoverable from the feature data alone.
    if baktfold_run is not None:
        data['baktfold_run'] = baktfold_run

    with json_path.open('wt') as fh:
        json.dump(data, fh, indent=4)
