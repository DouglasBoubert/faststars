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

from ..faststars import FASTSTARS
from ..utils import rgc_to_dhel


def do_ascii(catalog):
    """Process ASCII files extracted from datatables of published works."""
    task_str = catalog.get_current_task_str()

    # 2007ApJ...660..311B
    datafile = os.path.join(catalog.get_current_task_repo(), 'ASCII',
                            'ApJ_660_311_table1.csv')
    data = read(datafile, format='csv')
    for row in pbar(data, task_str):
        oname = 'SDSS'+str(row['SDSS'])[1:]
        name, source = catalog.new_entry(oname, bibcode='2007ApJ...660..311B')
        gallon = float(str(row['Glon']))
        gallat = float(str(row['Glat']))
        ra, dec = coord(
            l=gallon * u.degree, b=gallat * u.degree,
            frame='galactic').icrs.to_string(
                'hmsdms', sep=':').split()
        catalog.entries[name].add_quantity(
            FASTSTARS.RA, ra, source=source)
        catalog.entries[name].add_quantity(
            FASTSTARS.DEC, dec, source=source)
        catalog.entries[name].add_quantity(
            FASTSTARS.VELOCITY, str(row['Vhelio']), source=source)
        galrad_MS = float(str(row['Ra']))
        galrad_BHB = float(str(row['Rb']))
        dhel_MS = rgc_to_dhel(galrad_MS,gallon,gallat)
        dhel_BHB = rgc_to_dhel(galrad_BHB,gallon,gallat)
        catalog.entries[name].add_quantity(
            FASTSTARS.LUM_DIST, str(dhel_MS), lower_limit=str(dhel_BHB), source=source)
        catalog.entries[name].add_quantity(
            FASTSTARS.CLAIMED_TYPE, "pBHVS", source=source)
    catalog.journal_entries()

    return
