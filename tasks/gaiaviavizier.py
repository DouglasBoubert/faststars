"""Query against the Gaia database.
"""
import re

import numpy as np
from astropy.coordinates import SkyCoord as coord
import astropy.units as un
from astroquery.gaia import Gaia
from astroquery.vizier import Vizier
from astropy.table import vstack
import warnings

from astrocats.catalog.photometry import PHOTOMETRY
from astrocats.catalog.utils import is_number, pbar, single_spaces, uniq_cdl
from ..faststars import FASTSTARS
from ..utils import name_clean, parallax_to_distance

# astroquery.gaia is very verbose
import sys
from io import StringIO 

class NullIO(StringIO):
    def write(self, txt):
        pass


def silent(fn):
    """Decorator to silence functions."""
    def silent_fn(*args, **kwargs):
        saved_stdout = sys.stdout
        sys.stdout = NullIO()
        result = fn(*args, **kwargs)
        sys.stdout = saved_stdout
        return result
    return silent_fn
    
#silentgaiaquery = silent(Gaia.query_object)
v=Vizier(columns=['RA_ICRS','DE_ICRS','Plx','pmRA','pmDE','<Gmag>','e_RA_ICRS','e_DE_ICRS','e_Plx','e_pmRA','e_pmDE','RADEcor','RAPlxcor','RApmRAcor','RApmDEcor','DEPlxcor','DEpmRAcor','DEpmDEcor','PlxpmRAcor','PlxpmDEcor','pmRApmDEcor'])
silentgaiaquery = silent(v.query_region)


def do_gaiaviavizier(catalog):
    # Disable warnings
    warnings.filterwarnings("ignore")
    catalog.log.warning(
                'Some warnings are thrown by the Gaia module which do not affect the result of the Gaia queries.')
    
    task_str = catalog.get_current_task_str()
    keys = list(catalog.entries.keys())

    cntgphot = 0
    cntgast = 0
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
            Mradec=str(catalog.entries[name][FASTSTARS.RA][0]['value'])+str(catalog.entries[name][FASTSTARS.DEC][0]['value'])
            #print(name,Mradec)
            c=coord(Mradec,unit=(un.hourangle, un.deg),frame='icrs')

            #result = silentgaiaquery(c, radius=un.Quantity(2.0, un.arcsecond))
            result = silentgaiaquery(c, radius=un.Quantity(2.0, un.arcsecond), catalog='I/337/gaia')
            if len(result)>0:
                result=result[0]
                cntgphot += 1
                source = catalog.entries[name].add_source(bibcode='2016A&A...595A...2G')
                catalog.entries[name].add_photometry(
                        time=57023,
                        u_time='MJD',
                        telescope='Gaia',
                        band='G',
                        magnitude=str(result['__Gmag_'][0]),
                        source=source)
                if result['Plx'][0] != '--':
                    cntgast += 1
                    catalog.log.warning(
                        '"{}" has Gaia astrometry.'.format(name)) 
                    catalog.entries[name].add_quantity(FASTSTARS.PROPER_MOTION_RA, str(result['pmRA'][0]), source, e_value=str(result['e_pmRA'][0]), u_value='mas/yr')
                    catalog.entries[name].add_quantity(FASTSTARS.PROPER_MOTION_DEC, str(result['pmDE'][0]), source, e_value=str(result['e_pmDE'][0]), u_value='mas/yr')
                    catalog.entries[name].add_quantity(FASTSTARS.PARALLAX, str(result['Plx'][0]), source, e_value=str(result['e_Plx'][0]), u_value='mas')
                   
                    # Convert parallax to distance
                    if (FASTSTARS.LUM_DIST in catalog.entries[name]):
                        catalog.log.warning(
                        '"{}" has distance from photmetric distance prior.'.format(name)) 
                        distance, distance_error = parallax_to_distance(result['Plx'][0],result['e_Plx'][0],float(catalog.entries[name][FASTSTARS.LUM_DIST][0]['value']),float(catalog.entries[name][FASTSTARS.LUM_DIST][0]['e_value']))
                        catalog.entries[name].add_quantity(FASTSTARS.LUM_DIST, str(distance), e_value=str(distance_error), u_value='kpc', source=source, derived=True)
                    else:
                        catalog.log.warning(
                        '"{}" has distance from Astraatmadja & Bailer-Jones (2016) prior.'.format(name)) 
                        distance, distance_error = parallax_to_distance(result['Plx'][0],result['e_Plx'][0])
                        catalog.entries[name].add_quantity(FASTSTARS.LUM_DIST, str(distance), e_value=str(distance_error), u_value='kpc', source=source, derived=True)
                        
    catalog.log.warning(
                '"{}" have Gaia photometry and "{}" have Gaia astrometry.'.format(cntgphot,cntgast)) 
    catalog.journal_entries()
    
    # Reactivate warnings
    warnings.filterwarnings("default")
    return
