import json
from loguru import logger

from collections import OrderedDict
from pathlib import Path
from typing import Sequence

import baktfold.bakta
import baktfold.bakta.constants as bc
import baktfold.bakta.config as cfg





def write_json(data: dict, features: Sequence[dict], json_path: Path):
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

    version = OrderedDict()
    version['bakta'] = cfg.version
    version['db'] = cfg.version
    # version['db'] = {
    #     'version': f"{cfg.db_info['major']}.{cfg.db_info['minor']}",
    #     'type': cfg.db_info['type']
    # }
    data['version'] = version

    with json_path.open('wt') as fh:
        json.dump(data, fh, indent=4)
