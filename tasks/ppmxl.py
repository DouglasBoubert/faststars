"""Query against PPMXL to acquire proper motions."""
import re

import numpy as np
from astropy.coordinates import SkyCoord as coord
import astropy.units as un
from astroquery.irsa import Irsa
from astropy.table import vstack
import warnings
import time

from astrocats.catalog.photometry import PHOTOMETRY
from astrocats.catalog.utils import is_number, pbar, single_spaces, uniq_cdl
from ..faststars import FASTSTARS
from ..utils import name_clean

degperyrtomasperyear = 1e3*3600.

def do_ppmxl(catalog):
    task_str = catalog.get_current_task_str()
    keys = list(catalog.entries.keys())
    
    warnings.filterwarnings("ignore")
    
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
            radec= str(catalog.entries[name][FASTSTARS.RA][0]['value'])+str(catalog.entries[name][FASTSTARS.DEC][0]['value'])
            c=coord(radec,unit=(un.hourangle, un.deg),frame='icrs')

            cnttry = 0
            foundstar = False
            while foundstar == False and cnttry < 1:
                try:
                    cnttry += 1
                    time.sleep(0.1)
                    result = Irsa.query_region(c,catalog='ppmxl',radius='0d0m10s')
                except TypeError:
                    #print(radec,cnttry)
                    continue
                if len(result) > 1:
                    foundstar = True
                    
            if foundstar == True:
                source = (catalog.entries[name]
                          .add_source(name='The PPMXL Catalog',
                                      bibcode="2010AJ....139.2440R",
                                      url="https://irsa.ipac.caltech.edu/Missions/ppmxl.html",
                                      secondary=True))
                catalog.entries[name].add_quantity(FASTSTARS.PROPER_MOTION_RA,  str(result['pmra'][0]*degperyrtomasperyear), source, e_value=str(result['e_pmra'][0]*degperyrtomasperyear), u_value='mas/yr')
                catalog.entries[name].add_quantity(FASTSTARS.PROPER_MOTION_DEC, str(result['pmde'][0]*degperyrtomasperyear), source, e_value=str(result['e_pmde'][0]*degperyrtomasperyear), u_value='mas/yr')
                
    catalog.journal_entries()
    return
