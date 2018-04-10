"""Query against the Gaia database.
"""
import re

import numpy as np
from astropy.coordinates import SkyCoord as coord
import astropy.units as un
from astroquery.gaia import Gaia
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
    
silentgaiaquery = silent(Gaia.query_object)


def do_gaia(catalog):
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
            c=coord(Mradec,unit=(un.hourangle, un.deg),frame='icrs')

            result = silentgaiaquery(c, radius=un.Quantity(2.0, un.arcsecond))
            if len(result['dist'])>0:
                cntgphot += 1
                source = catalog.entries[name].add_source(bibcode='2016A&A...595A...2G')
                catalog.entries[name].add_photometry(
                        time=55000,
                        u_time='MJD',
                        telescope='Gaia',
                        band='G',
                        magnitude=str(result['phot_g_mean_mag'][0]),
                        source=source)
                if result['parallax'][0] != '--':
                    cntgast += 1
                    catalog.log.warning(
                        '"{}" has Gaia astrometry.'.format(name)) 
                    catalog.entries[name].add_quantity(FASTSTARS.PROPER_MOTION_RA, str(result['pmra'][0]), source, e_value=str(result['pmra_error'][0]), u_value='mas/yr')
                    catalog.entries[name].add_quantity(FASTSTARS.PROPER_MOTION_DEC, str(result['pmdec'][0]), source, e_value=str(result['pmdec_error'][0]), u_value='mas/yr')
                    catalog.entries[name].add_quantity(FASTSTARS.PARALLAX, str(result['parallax'][0]), source, e_value=str(result['parallax_error'][0]), u_value='mas')
                   
                    # Convert parallax to distance
                    if (FASTSTARS.LUM_DIST in catalog.entries[name]):
                        catalog.log.warning(
                        '"{}" has distance from photmetric distance prior.'.format(name)) 
                        distance, distance_error = parallax_to_distance(result['parallax'][0],result['parallax_error'][0],float(catalog.entries[name][FASTSTARS.LUM_DIST][0]['value']),float(catalog.entries[name][FASTSTARS.LUM_DIST][0]['e_value']))
                        catalog.entries[name].add_quantity(FASTSTARS.LUM_DIST, str(distance), e_value=str(distance_error), u_value='kpc', source=source, derived=True)
                    else:
                        catalog.log.warning(
                        '"{}" has distance from Astraatmadja & Bailer-Jones (2016) prior.'.format(name)) 
                        distance, distance_error = parallax_to_distance(result['parallax'][0],result['parallax_error'][0])
                        catalog.entries[name].add_quantity(FASTSTARS.LUM_DIST, str(distance), e_value=str(distance_error), u_value='kpc', source=source, derived=True)
                        
    catalog.log.warning(
                '"{}" have Gaia photometry and "{}" have Gaia astrometry.'.format(cntgphot,cntgast)) 
    catalog.journal_entries()
    
    # Reactivate warnings
    warnings.filterwarnings("default")
    return
