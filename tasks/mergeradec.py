"""Merge stars based on close ra/dec prior to cleanup."""
import re
import statistics
import warnings
from decimal import Decimal
from math import log10, pi, sqrt
import numpy as np

from astrocats.catalog.quantity import QUANTITY
from astrocats.catalog.utils import (get_sig_digits, is_number, pbar,
                                     pretty_num, tprint, uniq_cdl)
from astropy import units as un
from astropy.coordinates import SkyCoord as coord
from astropy.cosmology import Planck15 as cosmo
from astropy.cosmology import z_at_value
from scipy.spatial import cKDTree

from ..constants import CLIGHT, KM
from ..faststars import FASTSTARS


def do_mergeradec(catalog):
    """Merge stars based on close ra/dec prior to cleanup."""
    task_str = catalog.get_current_task_str()

    # Set preferred names, calculate some columns based on imported data,
    # sanitize some fields
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

        # Set the preferred name, switching to that name if name changed.
        #name = catalog.entries[name].set_preferred_name()

        #aliases = catalog.entries[name].get_aliases()

        if (FASTSTARS.RA not in catalog.entries[name] or
                FASTSTARS.DEC not in catalog.entries[name]):
            continue
        else:
            Mnames.append(name)
            Mradec.append(str(catalog.entries[name][FASTSTARS.RA][0]['value'])+str(catalog.entries[name][FASTSTARS.DEC][0]['value']))
    
    c=coord(Mradec,unit=(un.hourangle, un.deg))
    
    ### Construct tree to look for duplicates
    # Convert to pixel space
    pixx = np.cos(c.dec.rad)*np.cos(c.ra.rad)
    pixy = np.cos(c.dec.rad)*np.sin(c.ra.rad)
    pixz = np.sin(c.dec.rad)
    pix = np.vstack([pixx,pixy,pixz]).T

    # Construct tree
    pixtree = cKDTree(pix)

    # Query pairs
    pixarcsec = np.tan(5.*np.pi/180./3600.)
    pairs=pixtree.query_pairs(pixarcsec)
    arrpairs = np.array(list(pairs))
    
    # Attempt to merge pairs
    for i,j in arrpairs:
        try:
            catalog.copy_entry_to_entry(
                catalog.entries[Mnames[i]], catalog.entries[Mnames[j]])
            del catalog.entries[Mnames[i]]
            print("`{}` and `{}` merged".
                       format(Mnames[i],Mnames[j]))
        except KeyError:
            print("`{}` and `{}` pair already broken".
                       format(Mnames[i],Mnames[j]))
    
    catalog.save_caches()

    return
