"""Query against Simbad to acquire aliases."""
import re

import numpy as np
from astropy.coordinates import SkyCoord as coord
import astropy.units as un
from astroquery.simbad import Simbad
from astropy.table import vstack
import warnings

from astrocats.catalog.photometry import PHOTOMETRY
from astrocats.catalog.utils import is_number, pbar, single_spaces, uniq_cdl
from ..faststars import FASTSTARS
from ..utils import name_clean


def do_simbad(catalog):
    task_str = catalog.get_current_task_str()
    keys = list(catalog.entries.keys())
    
    customSimbad = Simbad()
    customSimbad.ROW_LIMIT = -1
    customSimbad.TIMEOUT = 120
    customSimbad.add_votable_fields('otype', 'sptype', 'sp_bibcode', 'id')

    Mnames, Mradec = [], []
    for oname in pbar(keys, task_str):
        # Some events may be merged in cleanup process, skip them if
        # non-existent.
        try:
            name = catalog.add_entry(oname)
        except Exception:
            catalog.log.warning(
                '"{}" was not found, suggests merge occurred in cleanup '
                'process.'.format(oname))
            continue

        if (FASTSTARS.RA not in catalog.entries[name] or
                FASTSTARS.DEC not in catalog.entries[name]):
            continue
        else:
            Mnames.append(name)
            radec= str(catalog.entries[name][FASTSTARS.RA][0]['value'])+str(catalog.entries[name][FASTSTARS.DEC][0]['value'])
            c=coord(radec,unit=(un.hourangle, un.deg),frame='icrs')

            cnttry = 0
            failedtofindstar = False
            with warnings.catch_warnings():
                warnings.filterwarnings('error')
                while cnttry >-0.5 and cnttry < 10:
                    try:
                        cnttry += 1
                        result = customSimbad.query_region(c,radius='0d0m10s')
                        aliases = re.sub(r'b\'(.*)\'', r'\1',str(result['ID'].tolist()[0])).split(',')
                    except TypeError: 
                        continue
                    except UserWarning:
                        failedtofindstar = True
                        catalog.log.warning(
                            '"{}" was not found in Simbad.'.format(name)) 
                    cnttry = -1.0
                    
            if failedtofindstar: continue
            source = (catalog.entries[name]
                      .add_source(name='SIMBAD astronomical database',
                                  bibcode="2000A&AS..143....9W",
                                  url="http://simbad.u-strasbg.fr/",
                                  secondary=True))
            for alias in aliases:
                ali = single_spaces(re.sub(r'\[[^)]*\]', '', alias).strip())
                if is_number(ali.replace(' ', '')):
                    continue
                #if ali in simbadbannednames:
                #    continue
                ali = name_clean(ali)
                catalog.entries[name].add_quantity(FASTSTARS.ALIAS,
                                                   ali, source)
                
    catalog.journal_entries()
    return
