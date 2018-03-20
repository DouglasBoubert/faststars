# -*- coding: utf-8 -*-
"""ASCII datafiles.

Often produced from LaTeX tables in the original papers,
but sometimes provided as supplementary datafiles on the journal webpages.
"""
import csv
import os
import re
from datetime import datetime
from decimal import Decimal
from glob import glob

from astrocats.catalog.photometry import PHOTOMETRY, set_pd_mag_from_counts
from astrocats.catalog.utils import (is_number, jd_to_mjd, make_date_string,
                                     pbar, pbar_strings)
from astropy import units as u
from astropy.coordinates import SkyCoord as coord
from astropy.io.ascii import read
from astropy.time import Time as astrotime

from ..supernova import SUPERNOVA


def do_ascii(catalog):
    """Process ASCII files extracted from datatables of published works."""
    task_str = catalog.get_current_task_str()

    # 2017ApJ...836...60L
    datafile = os.path.join(catalog.get_current_task_repo(), 'ASCII',
                            '2017ApJ...836...60L-tab1.cds')
    data = read(datafile, format='cds')
    for row in pbar(data, task_str):
        oname = row['ID']
        name, source = catalog.new_entry(oname, bibcode='2017ApJ...836...60L')
        photodict = {
            PHOTOMETRY.TIME: str(row['MJD']),
            PHOTOMETRY.BAND: row['Filter'],
            PHOTOMETRY.U_TIME: 'MJD',
            PHOTOMETRY.MAGNITUDE: row['mag'],
            PHOTOMETRY.E_MAGNITUDE: row['e_mag'],
            PHOTOMETRY.SOURCE: source
        }
        catalog.entries[name].add_photometry(**photodict)
    catalog.journal_entries()

    return
