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
        oname = 'SDSS'+str(row['SDSS'])
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
            FASTSTARS.LUM_DIST, str(dhel_MS), lower_limit=str(dhel_BHB), u_value='kpc', source=source)
        catalog.entries[name].add_quantity(
            FASTSTARS.CLAIMED_TYPE, "pBHVS", source=source)
    catalog.journal_entries()
    
    # 2009ApJ...690.1369B
    datafile = os.path.join(catalog.get_current_task_repo(), 'ASCII',
                            'apj292642t1_ascii.csv')
    data = read(datafile, format='csv')
    for row in pbar(data, task_str):
        oname = row['Catalog']
        name, source = catalog.new_entry(oname, bibcode='2009ApJ...690.1369B')
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
            FASTSTARS.VELOCITY, str(row['Vhel']), e_value='17.0', source=source) # The 17.0 is the upper bound of the systematic error.
        galrad = float(str(row['RGC']))
        dhel = rgc_to_dhel(galrad,gallon,gallat)
        catalog.entries[name].add_quantity(
            FASTSTARS.LUM_DIST, str(dhel), u_value='kpc', source=source)
        catalog.entries[name].add_quantity(
            FASTSTARS.SPECTRAL_TYPE, row['Type'], source=source)
        if str(row['ID'])[:3]=='HVS':
            catalog.entries[name].add_quantity(
                FASTSTARS.CLAIMED_TYPE, "HVS", source=source)
            catalog.entries[name].add_quantity(
                FASTSTARS.ALIAS, row['ID'], source=source)
        else:
            catalog.entries[name].add_quantity(
                FASTSTARS.CLAIMED_TYPE, row['ID'], source=source)
    catalog.journal_entries()
    
    # 2012ApJ...744L..24L
    datafile = os.path.join(catalog.get_current_task_repo(), 'ASCII',
                            'apjl415156t1t2.csv')
    data = read(datafile, format='csv')
    for row in pbar(data, task_str):
        oname = str(row['Catalog'])
        name, source = catalog.new_entry(oname, bibcode='2012ApJ...744L..24L')
        catalog.entries[name].add_quantity(
                FASTSTARS.ALIAS, row['ID'], source=source)
        radec = oname.strip('SDSSJ')
        print(radec)
        radec = radec[0:2]+' '+radec[2:4]+' '+radec[4:9]+' '+radec[9:12]+' '+radec[12:14]+' '+radec[14:18]
        print(radec)
        ra, dec = coord(radec, 
                unit=(u.hourangle, u.deg)).to_string(
                'hmsdms', sep=':').split()
        catalog.entries[name].add_quantity(
            FASTSTARS.RA, ra, source=source)
        catalog.entries[name].add_quantity(
            FASTSTARS.DEC, dec, source=source)
        catalog.entries[name].add_quantity(
            FASTSTARS.VELOCITY, str(row['Vhel']), e_value=str(row['e_Vhel']), source=source)
        catalog.entries[name].add_quantity(
            FASTSTARS.LUM_DIST, str(row['Dhel']), e_value=str(row['e_Dhel']), u_value='kpc', source=source)
        catalog.entries[name].add_quantity(
            FASTSTARS.PROPER_MOTION_RA, str(row['pmra']), e_value=str(row['e_pmra']), u_value='mas/yr', source=source)
        catalog.entries[name].add_quantity(
            FASTSTARS.PROPER_MOTION_DEC, str(row['pmdec']), e_value=str(row['e_pmdec']), u_value='mas/yr', source=source)
        catalog.entries[name].add_quantity(
            FASTSTARS.CLAIMED_TYPE, "pHVS", source=source)
        catalog.entries[name].add_quantity(
            FASTSTARS.SPECTRAL_TYPE, "F", source=source)
    catalog.journal_entries()

    return
