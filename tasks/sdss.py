"""Query against the SDSS database.
"""
import re

import numpy as np
from astropy.coordinates import SkyCoord as coord
import astropy.units as un
from astroquery.sdss import SDSS
from astropy.table import vstack

from astrocats.catalog.photometry import PHOTOMETRY
from astrocats.catalog.utils import is_number, pbar, single_spaces, uniq_cdl
from ..faststars import FASTSTARS
from ..utils import name_clean


def do_sdss(catalog):
    task_str = catalog.get_current_task_str()
    keys = list(catalog.entries.keys())

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
            Mradec.append(str(catalog.entries[name][FASTSTARS.RA][0]['value'])+str(catalog.entries[name][FASTSTARS.DEC][0]['value']))
    
    c=coord(Mradec,unit=(un.hourangle, un.deg),frame='icrs')

    # We must step through the data to query it, >100 is too many
    maxstep = 100 # stepsize
    Nentries = len(Mradec)
    roundindex = np.zeros(0) # the stepping round that this star's data was acquired in, needed to connect the obj_id in the query to the name of the star
    for i in range(int(Nentries/maxstep)):
        result_tmp = SDSS.query_crossid(c[maxstep*i:min(maxstep*(i+1),Nentries)], timeout=200.,photoobj_fields=['u','err_u','g','err_g','r', 'err_r','i','err_i','z','err_z','MJD'])
        roundindex = np.concatenate([roundindex,i*np.ones(len(result_tmp['obj_id']))])
        if i == 0:
            result = result_tmp
        else:
            result = vstack([result,result_tmp])
            
    flagsuccess = result['obj_id']
    listfilter = ['u','g','r','i','z']
    for i in range(len(flagsuccess)):
        Mi = int(flagsuccess[i].strip('obj_'))+maxstep*int(roundindex[i])
        name = Mnames[Mi]
        source = catalog.entries[name].add_source(bibcode='2015ApJS..219...12A')
        for j in range(5):
            catalog.entries[name].add_photometry(
                time=str(result[i]['MJD']),
                u_time='MJD',
                telescope='SDSS',
                band=listfilter[j],
                magnitude=str(result[i][listfilter[j]]),
                e_magnitude=str(result[i]['err_'+listfilter[j]]),
                source=source)
                
    catalog.journal_entries()
    return
