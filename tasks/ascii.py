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
        if str(row['e_Vhel'])!='NA':
            catalog.entries[name].add_quantity(
                FASTSTARS.VELOCITY, str(row['Vhel']), e_value=str(row['e_Vhel']), source=source)
        else:
            catalog.entries[name].add_quantity(
                FASTSTARS.VELOCITY, str(row['Vhel']), source=source)
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
        radec = radec[0:2]+' '+radec[2:4]+' '+radec[4:9]+' '+radec[9:12]+' '+radec[12:14]+' '+radec[14:18]
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
    
    # 2012ApJ...751...55B
    datafile = os.path.join(catalog.get_current_task_repo(), 'ASCII',
                            'apj427101t1_ascii.csv')
    data = read(datafile, format='csv')
    for row in pbar(data, task_str):
        oname = str(row['Catalog']).replace(' ','')
        name, source = catalog.new_entry(oname, bibcode='2012ApJ...751...55B')
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
        if str(row['e_Vhel'])!='NA':
            catalog.entries[name].add_quantity(
                FASTSTARS.VELOCITY, str(row['Vhel']), e_value=str(row['e_Vhel']), source=source)
            catalog.entries[name].add_quantity(
                FASTSTARS.SPECTRAL_TYPE, row['Type'], source=source) # Only the new HVSs have definite spectral types
        else:
            catalog.entries[name].add_quantity(
                FASTSTARS.VELOCITY, str(row['Vhel']), source=source)
        galrad = float(str(row['RGC']))
        dhel = rgc_to_dhel(galrad,gallon,gallat)
        catalog.entries[name].add_quantity(
            FASTSTARS.LUM_DIST, str(dhel), u_value='kpc', source=source)
        if str(row['ID'])[:3]=='HVS':
            catalog.entries[name].add_quantity(
                FASTSTARS.CLAIMED_TYPE, "HVS", source=source)
            catalog.entries[name].add_quantity(
                FASTSTARS.ALIAS, row['ID'], source=source)
        else:
            catalog.entries[name].add_quantity(
                FASTSTARS.CLAIMED_TYPE, row['ID'], source=source)
    catalog.journal_entries()
    
    # 2014ApJ...780....7P
    datafile = os.path.join(catalog.get_current_task_repo(), 'ASCII',
                            'apj485719t1_ascii.csv')
    data = read(datafile, format='csv')
    for row in pbar(data, task_str):
        oname = 'SDSS'+str(row['Catalog']).replace(' ','')
        name, source = catalog.new_entry(oname, bibcode='2014ApJ...780....7P')
        catalog.entries[name].add_quantity(
                FASTSTARS.ALIAS, 'Pal'+str(row['Pal']), source=source)
        radec = oname.strip('SDSSJ')
        radec = radec[0:2]+' '+radec[2:4]+' '+radec[4:9]+' '+radec[9:12]+' '+radec[12:14]+' '+radec[14:18]
        ra, dec = coord(radec, 
                unit=(u.hourangle, u.deg)).to_string(
                'hmsdms', sep=':').split()
        catalog.entries[name].add_quantity(
            FASTSTARS.RA, ra, source=source)
        catalog.entries[name].add_quantity(
            FASTSTARS.DEC, dec, source=source)
        catalog.entries[name].add_quantity(
            FASTSTARS.VELOCITY, str(row['Vhel']), source=source)
        catalog.entries[name].add_quantity(
            FASTSTARS.LUM_DIST, str(row['Dhel']), u_value='kpc', source=source)
        catalog.entries[name].add_quantity(
            FASTSTARS.CLAIMED_TYPE, "pHVS", source=source)
        catalog.entries[name].add_quantity(
            FASTSTARS.SPECTRAL_TYPE, "G/K", source=source)
    catalog.journal_entries()
    
    # 2014ApJ...787...89B
    datafile = os.path.join(catalog.get_current_task_repo(), 'ASCII',
                            'apj494602t1_ascii.csv')
    data = read(datafile, format='csv')
    for row in pbar(data, task_str):
        oname = str(row['Catalog']).replace(' ','')
        name, source = catalog.new_entry(oname, bibcode='2014ApJ...787...89B')
        radec = oname.strip('SDSSJ')
        radec = radec[0:2]+' '+radec[2:4]+' '+radec[4:9]+' '+radec[9:12]+' '+radec[12:14]+' '+radec[14:19]
        ra, dec = coord(radec, 
                unit=(u.hourangle, u.deg)).to_string(
                'hmsdms', sep=':').split()
        catalog.entries[name].add_quantity(
            FASTSTARS.RA, ra, source=source)
        catalog.entries[name].add_quantity(
            FASTSTARS.DEC, dec, source=source)
        catalog.entries[name].add_quantity(
                FASTSTARS.VELOCITY, str(row['Vhel']), e_value=str(row['e_Vhel']), source=source)
        galrad = float(str(row['RGC']))
        errgalrad = float(str(row['e_RGC']))
        dhel = rgc_to_dhel(galrad,gallon,gallat)
        dhel_lo = rgc_to_dhel(galrad-errgalrad,gallon,gallat)
        dhel_hi = rgc_to_dhel(galrad+errgalrad,gallon,gallat)
        catalog.entries[name].add_quantity(
            FASTSTARS.LUM_DIST, str(dhel), e_lower_value=dhel-dhel_lo, e_upper_value=dhel_hi-dhel, u_value='kpc', source=source)
        if str(row['ID'])!='pBHVS':
            catalog.entries[name].add_quantity(
                FASTSTARS.CLAIMED_TYPE, "HVS", source=source)
            catalog.entries[name].add_quantity(
                FASTSTARS.ALIAS, 'HVS'+str(row['ID']), source=source)
        else:
            catalog.entries[name].add_quantity(
                FASTSTARS.CLAIMED_TYPE, 'pBHVS', source=source)
    catalog.journal_entries()
    
    # 2014EAS....67..255Z
    datafile = os.path.join(catalog.get_current_task_repo(), 'ASCII',
                            '1501.07824.csv')
    data = read(datafile, format='csv')
    for row in pbar(data, task_str):
        oname = str(row['ID'])
        name, source = catalog.new_entry(oname, bibcode='2014EAS....67..255Z')
        catalog.entries[name].add_quantity(
                FASTSTARS.VELOCITY, str(row['Vhel']), e_value=str(row['e_Vhel']), source=source)
        catalog.entries[name].add_quantity(
            FASTSTARS.PROPER_MOTION_RA, str(row['pmra']), e_value=str(row['e_pmra']), u_value='mas/yr', source=source)
        catalog.entries[name].add_quantity(
            FASTSTARS.PROPER_MOTION_DEC, str(row['pmdec']), e_value=str(row['e_pmdec']), u_value='mas/yr', source=source)
        catalog.entries[name].add_quantity(
            FASTSTARS.LUM_DIST, str(row['Dhel']), u_value='kpc', source=source) # This distance may have some metallicity dependent uncertainty?
        catalog.entries[name].add_quantity(
                FASTSTARS.CLAIMED_TYPE, 'pHVS', source=source)
        catalog.entries[name].add_quantity(
                FASTSTARS.SPECTRAL_TYPE, 'K/M', source=source)
    catalog.journal_entries()
    
    # 2014ApJ...789L...2Z
    datafile = os.path.join(catalog.get_current_task_repo(), 'ASCII',
                            'apjl496832t1t2t3_ascii.csv')
    data = read(datafile, format='csv')
    for row in pbar(data, task_str):
        oname = str(row['ID'])
        name, source = catalog.new_entry(oname, bibcode='2014ApJ...789L...2Z')
        catalog.entries[name].add_quantity(
            FASTSTARS.ALIAS, str(row['Catalog']), source=source)
        radec = str(row['Catalog']).strip('J')
        radec = radec[0:2]+' '+radec[2:4]+' '+radec[4:9]+' '+radec[9:12]+' '+radec[12:14]+' '+radec[14:18]
        ra, dec = coord(radec, 
                unit=(u.hourangle, u.deg)).to_string(
                'hmsdms', sep=':').split()
        catalog.entries[name].add_quantity(
            FASTSTARS.RA, ra, source=source)
        catalog.entries[name].add_quantity(
            FASTSTARS.DEC, dec, source=source)
        catalog.entries[name].add_quantity(
                FASTSTARS.VELOCITY, str(row['Vrb']), source=source)
        catalog.entries[name].add_quantity(
            FASTSTARS.PROPER_MOTION_RA, str(row['pmra']), e_value=str(row['e_pmra']), u_value='mas/yr', source=source)
        catalog.entries[name].add_quantity(
            FASTSTARS.PROPER_MOTION_DEC, str(row['pmdec']), e_value=str(row['e_pmdec']), u_value='mas/yr', source=source)
        catalog.entries[name].add_quantity(
            FASTSTARS.LUM_DIST, str(row['Dhel']), u_value='kpc', source=source) # This distance may have some metallicity dependent uncertainty?
        catalog.entries[name].add_quantity(
                FASTSTARS.CLAIMED_TYPE, 'pHVS', source=source)
        catalog.entries[name].add_quantity(
                FASTSTARS.SPECTRAL_TYPE, str(row['Type']), source=source)
    catalog.journal_entries()
    
    # 2014ApJ...794..146T
    datafile = os.path.join(catalog.get_current_task_repo(), 'ASCII',
                            'ApJ794146.csv')
    data = read(datafile, format='csv')
    for row in pbar(data, task_str):
        oname = str(row['Catalog']).replace(' ','')
        name, source = catalog.new_entry(oname, bibcode='2014ApJ...794..146T')
        radec = oname.strip('SDSSJ')
        radec = radec[0:2]+' '+radec[2:4]+' '+radec[4:9]+' '+radec[9:12]+' '+radec[12:14]+' '+radec[14:18]
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
            FASTSTARS.PROPER_MOTION_RA, str(row['pmra']), e_value=str(row['e_pmra']), u_value='mas/yr', source=source)
        catalog.entries[name].add_quantity(
            FASTSTARS.PROPER_MOTION_DEC, str(row['pmdec']), e_value=str(row['e_pmdec']), u_value='mas/yr', source=source)
        catalog.entries[name].add_quantity(
            FASTSTARS.LUM_DIST, str(row['Dhel']), e_value=str(row['e_Dhel']), u_value='kpc', source=source)
        catalog.entries[name].add_quantity(
                FASTSTARS.CLAIMED_TYPE, 'pHVS', source=source)
        catalog.entries[name].add_quantity(
                FASTSTARS.SPECTRAL_TYPE, str(row['Type']), source=source)
    catalog.journal_entries()
    
    # 2014ApJ...794..145S
    datafile = os.path.join(catalog.get_current_task_repo(), 'ASCII',
                            'apj501503t8_ascii.csv')
    data = read(datafile, format='csv')
    for row in pbar(data, task_str):
        oname = str(row['Catalog'])
        oname = oname[:4]+'J'+oname[4:]
        name, source = catalog.new_entry(oname, bibcode='2014ApJ...794..145S')
        radec = oname.strip('SDSSJ')
        radec = radec[0:2]+' '+radec[2:4]+' '+radec[4:8]+' '+radec[8:11]+' '+radec[11:13]+' '+radec[13:17]
        ra, dec = coord(radec, 
                unit=(u.hourangle, u.deg)).to_string(
                'hmsdms', sep=':').split()
        catalog.entries[name].add_quantity(
            FASTSTARS.RA, ra, source=source)
        catalog.entries[name].add_quantity(
            FASTSTARS.DEC, dec, source=source)
        catalog.entries[name].add_quantity(
                FASTSTARS.VELOCITY, str(row['Vhel']), source=source)
        catalog.entries[name].add_quantity(
            FASTSTARS.PROPER_MOTION_RA, str(row['pmra']), e_value=str(row['e_pmra']), u_value='mas/yr', source=source)
        catalog.entries[name].add_quantity(
            FASTSTARS.PROPER_MOTION_DEC, str(row['pmdec']), e_value=str(row['e_pmdec']), u_value='mas/yr', source=source)
        catalog.entries[name].add_quantity(
            FASTSTARS.LUM_DIST, str(row['Dhel']), e_value=str(row['e_Dhel']), u_value='kpc', source=source)
        catalog.entries[name].add_quantity(
                FASTSTARS.CLAIMED_TYPE, 'pHVS', source=source)
        catalog.entries[name].add_quantity(
                FASTSTARS.SPECTRAL_TYPE, str(row['Type']), source=source)
    catalog.journal_entries()
    
    # 2015MNRAS.447.2046H
    datafile = os.path.join(catalog.get_current_task_repo(), 'ASCII',
                            'hawkins2015.csv')
    data = read(datafile, format='csv')
    for row in pbar(data, task_str):
        oname = str(row['Catalog']).strip(' ')
        oname = 'RAVE'+oname
        name, source = catalog.new_entry(oname, bibcode='2015MNRAS.447.2046H')
        radec = oname.strip('RAVEJ')
        radec = radec[0:2]+' '+radec[2:4]+' '+radec[4:8]+' '+radec[8:11]+' '+radec[11:13]+' '+radec[13:15]
        ra, dec = coord(radec, 
                unit=(u.hourangle, u.deg)).to_string(
                'hmsdms', sep=':').split()
        catalog.entries[name].add_quantity(
            FASTSTARS.RA, ra, source=source)
        catalog.entries[name].add_quantity(
            FASTSTARS.DEC, dec, source=source)
        catalog.entries[name].add_quantity(
                FASTSTARS.VELOCITY, str(row['Vhel']).strip(' '), e_value=str(row['e_Vhel']).strip(' '), source=source)
        catalog.entries[name].add_quantity(
            FASTSTARS.PROPER_MOTION_RA, str(row['pmra']).strip(' '), e_value=str(row['e_pmra']).strip(' '), u_value='mas/yr', source=source)
        catalog.entries[name].add_quantity(
            FASTSTARS.PROPER_MOTION_DEC, str(row['pmdec']).strip(' '), e_value=str(row['e_pmdec']).strip(' '), u_value='mas/yr', source=source)
        catalog.entries[name].add_quantity(
            FASTSTARS.LUM_DIST, str(row['Dhel']).strip(' '), e_value=str(row['e_Dhel']).strip(' '), u_value='pc', source=source)
        catalog.entries[name].add_quantity(
                FASTSTARS.CLAIMED_TYPE, str(row['Claim']).strip(' '), source=source)
    catalog.journal_entries()
    
    # 2015A&A...576L..14Z
    datafile = os.path.join(catalog.get_current_task_repo(), 'ASCII',
                            'zeigerer.csv')
    data = read(datafile, format='csv')
    for row in pbar(data, task_str):
        oname = 'Pal'+str(row['Pal'])
        name, source = catalog.new_entry(oname, bibcode='2015A&A...576L..14Z')
        catalog.entries[name].add_quantity(
            FASTSTARS.PROPER_MOTION_RA, str(row['pmra']), e_value=str(row['e_pmra']), u_value='mas/yr', source=source)
        catalog.entries[name].add_quantity(
            FASTSTARS.PROPER_MOTION_DEC, str(row['pmdec']), e_value=str(row['e_pmdec']), u_value='mas/yr', source=source)
        catalog.entries[name].add_quantity(
                FASTSTARS.CLAIMED_TYPE, 'nHVS', source=source)
    catalog.journal_entries()
    
    # 2015ApJ...804...49B
    datafile = os.path.join(catalog.get_current_task_repo(), 'ASCII',
                            'apj510826t1_ascii.csv')
    data = read(datafile, format='csv')
    for row in pbar(data, task_str):
        oname = str(row['Catalog']).strip(' ')
        if oname[1]==':':
            # Add the leading 0
            oname = '0'+oname
        oname = 'SDSSJ'+oname
        name, source = catalog.new_entry(oname.replace(':',''), bibcode='2015ApJ...804...49B')
        #print(str(row['ID']))
        catalog.entries[name].add_quantity(
            FASTSTARS.ALIAS, str(row['ID']), source=source)
        catalog.entries[name].add_quantity(
            FASTSTARS.PROPER_MOTION_RA, str(row['pmra']).strip(' '), e_value=str(row['e_pmra']).strip(' '), u_value='mas/yr', source=source)
        catalog.entries[name].add_quantity(
            FASTSTARS.PROPER_MOTION_DEC, str(row['pmdec']).strip(' '), e_value=str(row['e_pmdec']).strip(' '), u_value='mas/yr', source=source)
        if str(row['newspec']) == 'y':
            radec = oname.strip('SDSSJ')
            ra, dec = radec[:10], radec[10:]
            catalog.entries[name].add_quantity(
                FASTSTARS.RA, ra, source=source)
            catalog.entries[name].add_quantity(
                FASTSTARS.DEC, dec, source=source)
            catalog.entries[name].add_quantity(
                FASTSTARS.VELOCITY, str(row['Vhel']), e_value=str(row['e_Vhel']), source=source)
            catalog.entries[name].add_quantity(
                FASTSTARS.LUM_DIST, str(row['Dhel']), e_value=str(row['e_Dhel']), u_value='kpc', source=source)
            catalog.entries[name].add_quantity(
                    FASTSTARS.CLAIMED_TYPE, str(row['Claim']), source=source)
            catalog.entries[name].add_quantity(
                    FASTSTARS.SPECTRAL_TYPE, str(row['Type']), source=source)
    catalog.journal_entries()
    
    # 2015RAA....15.1364L
    datafile = os.path.join(catalog.get_current_task_repo(), 'ASCII',
                            'li2015.csv')
    data = read(datafile, format='csv')
    for row in pbar(data, task_str):
        oname = str(row['Catalog'])
        oname = 'LAMOST'+oname
        name, source = catalog.new_entry(oname, bibcode='2015RAA....15.1364L')
        radec = oname.strip('LAMOSTJ')
        radec = radec[0:2]+' '+radec[2:4]+' '+radec[4:9]+' '+radec[9:12]+' '+radec[12:14]+' '+radec[14:]
        ra, dec = coord(radec, 
                unit=(u.hourangle, u.deg)).to_string(
                'hmsdms', sep=':').split()
        catalog.entries[name].add_quantity(
            FASTSTARS.ALIAS, 'Li'+str(row['ID']), source=source)
        catalog.entries[name].add_quantity(
            FASTSTARS.RA, ra, source=source)
        catalog.entries[name].add_quantity(
            FASTSTARS.DEC, dec, source=source)
        catalog.entries[name].add_quantity(
                FASTSTARS.VELOCITY, str(row['Vhel']), e_value=str(row['e_Vhel']), source=source)
        catalog.entries[name].add_quantity(
            FASTSTARS.PROPER_MOTION_RA, str(row['pmra']), e_value=str(row['e_pmra']), u_value='mas/yr', source=source)
        catalog.entries[name].add_quantity(
            FASTSTARS.PROPER_MOTION_DEC, str(row['pmdec']), e_value=str(row['e_pmdec']), u_value='mas/yr', source=source)
        catalog.entries[name].add_quantity(
            FASTSTARS.LUM_DIST, str(row['Dhel']), e_value=str(row['e_Dhel']), u_value='kpc', source=source)
        catalog.entries[name].add_quantity(
                FASTSTARS.CLAIMED_TYPE, 'pHVS', source=source)
        catalog.entries[name].add_quantity(
                FASTSTARS.SPECTRAL_TYPE, 'F/G/K', source=source)
    catalog.journal_entries()
    
    # 2015AJ....150...77V
    datafile = os.path.join(catalog.get_current_task_repo(), 'ASCII',
                            'vickers2015.csv')
    data = read(datafile, format='csv')
    for row in pbar(data, task_str):
        oname = str(row['Catalog'])
        oname = 'SDSS'+oname
        name, source = catalog.new_entry(oname, bibcode='2015AJ....150...77V')
        radec = oname.strip('SDSSJ')
        radec = radec[0:2]+' '+radec[2:4]+' '+radec[4:9]+' '+radec[9:12]+' '+radec[12:14]+' '+radec[14:]
        ra, dec = coord(radec, 
                unit=(u.hourangle, u.deg)).to_string(
                'hmsdms', sep=':').split()
        catalog.entries[name].add_quantity(
            FASTSTARS.RA, ra, source=source)
        catalog.entries[name].add_quantity(
            FASTSTARS.DEC, dec, source=source)
        catalog.entries[name].add_quantity(
                FASTSTARS.CLAIMED_TYPE, str(row['Claim']), source=source)
        catalog.entries[name].add_quantity(
                FASTSTARS.SPECTRAL_TYPE, 'F/G/K/M', source=source)
    catalog.journal_entries()
    
    return
