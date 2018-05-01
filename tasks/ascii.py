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

    # Shen unpublished
    #datafile = os.path.join(catalog.get_current_task_repo(), 'ASCII',
    #                        'shen-unpublished.csv')
    #with open(datafile) as f:
    #    data = []
    #    for row in f.readlines():
    #        cols = row.split(',')
    #        data.append({'name': 'GaiaDR2-' + cols[0], 'ra': cols[5], 'e_ra': cols[6], 'dec': cols[7], 'e_dec': cols[8]})###

    #for row in pbar(data, task_str):
    #    oname = str(row['name'])
    #    name, source = catalog.new_entry(oname, srcname='Shen et al. 2018')
    #    ra, dec = coord(float(row['ra']), float(row['dec']), frame='icrs', unit='deg').to_string('hmsdms', sep=':', precision=8, precision=8).split()
    #    catalog.entries[name].add_quantity(FASTSTARS.RA, ra, source, e_value=row['e_ra'], u_e_value='mas')
    #    catalog.entries[name].add_quantity(FASTSTARS.DEC, dec, source, e_value=row['e_dec'], u_e_value='mas')
    #    catalog.entries[name].add_quantity(FASTSTARS.DISCOVERER, 'Gaia', source)

    # 2007ApJ...660..311B
    datafile = os.path.join(catalog.get_current_task_repo(), 'ASCII',
                            'ApJ_660_311_table1.csv')
    data = read(datafile, format='csv')
    for row in pbar(data, task_str):
        oname = 'SDSS'+str(row['SDSS'])
        name, source = catalog.new_entry(oname, bibcode='2007ApJ...660..311B')
        if (FASTSTARS.DISCOVERER not in catalog.entries[name]):
            catalog.entries[name].add_quantity(FASTSTARS.DISCOVERER,'Warren R. Brown, Margaret J. Geller, Scott J. Kenyon, Michael J. Kurtz, Benjamin C. Bromley', source)
            catalog.entries[name].add_quantity(FASTSTARS.DISCOVER_DATE,str(2007), source)
        gallon = float(str(row['Glon']))
        gallat = float(str(row['Glat']))
        radec = oname.strip('SDSSJ')
        radec = radec[0:2]+' '+radec[2:4]+' '+radec[4:9]+' '+radec[9:12]+' '+radec[12:14]+' '+radec[14:18]
        ra, dec = coord(radec, 
                unit=(u.hourangle, u.deg)).to_string(
                'hmsdms', sep=':', precision=8).split()
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
            FASTSTARS.LUM_DIST, str(dhel_MS), u_value='kpc', source=source,derived=True)
        catalog.entries[name].add_quantity(
            FASTSTARS.LUM_DIST, str(dhel_BHB), u_value='kpc', source=source,derived=True)
        catalog.entries[name].add_quantity(
            FASTSTARS.SPECTRAL_TYPE, str(row['Sp']), source=source)
        catalog.entries[name].add_quantity(
            FASTSTARS.STELLAR_CLASS, "d", source=source)
        catalog.entries[name].add_quantity(
            FASTSTARS.STELLAR_CLASS, "bhb", source=source)
        catalog.entries[name].add_quantity(
            FASTSTARS.STELLAR_CLASS, "g", source=source)
    catalog.journal_entries()
    
    # 2009ApJ...690.1639B
    datafile = os.path.join(catalog.get_current_task_repo(), 'ASCII',
                            'apj292642t1_ascii.csv')
    data = read(datafile, format='csv')
    for row in pbar(data, task_str):
        oname = row['Catalog']
        name, source = catalog.new_entry(oname, bibcode='2009ApJ...690.1639B')
        if (FASTSTARS.DISCOVERER not in catalog.entries[name]):
            catalog.entries[name].add_quantity(FASTSTARS.DISCOVERER,'Warren R. Brown, Margaret J. Geller, Scott J. Kenyon', source)
            catalog.entries[name].add_quantity(FASTSTARS.DISCOVER_DATE,str(2009), source)
        gallon = float(str(row['Glon']))
        gallat = float(str(row['Glat']))
        if (oname!='US708') & (oname!='HE0437-5439'):
            radec = oname.strip('SDSSJ')
            radec = radec[0:2]+' '+radec[2:4]+' '+radec[4:9]+' '+radec[9:12]+' '+radec[12:14]+' '+radec[14:18]
            ra, dec = coord(radec, 
                unit=(u.hourangle, u.deg)).to_string(
                'hmsdms', sep=':', precision=8).split()
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
            FASTSTARS.LUM_DIST, str(dhel), u_value='kpc', source=source,derived=True)
        sptype = str(row['Type']).split('/')
        for SPTYPE in sptype:
            if SPTYPE == "sdO":
                catalog.entries[name].add_quantity(
                    FASTSTARS.SPECTRAL_TYPE, "O", source=source)
                catalog.entries[name].add_quantity(
                    FASTSTARS.STELLAR_CLASS, "sd", source=source)
            elif SPTYPE == "BHB":
                catalog.entries[name].add_quantity(
                    FASTSTARS.STELLAR_CLASS, "bhb", source=source)
            else:
                catalog.entries[name].add_quantity(
                    FASTSTARS.SPECTRAL_TYPE, SPTYPE, source=source)
                catalog.entries[name].add_quantity(
                    FASTSTARS.STELLAR_CLASS, "d", source=source)
        if str(row['ID'])[:3]=='HVS':
            catalog.entries[name].add_quantity(
                FASTSTARS.ALIAS, row['ID'], source=source)
    catalog.journal_entries()
    
    # 2012ApJ...744L..24L
    datafile = os.path.join(catalog.get_current_task_repo(), 'ASCII',
                            'apjl415156t1t2.csv')
    data = read(datafile, format='csv')
    for row in pbar(data, task_str):
        oname = str(row['Catalog'])
        name, source = catalog.new_entry(oname, bibcode='2012ApJ...744L..24L')
        if (FASTSTARS.DISCOVERER not in catalog.entries[name]):
            catalog.entries[name].add_quantity(FASTSTARS.DISCOVERER,'Yinbi Li, Ali Luo, Gang Zhao, Youjun Lu, Juanjuan Ren, Fang Zuo', source)
            catalog.entries[name].add_quantity(FASTSTARS.DISCOVER_DATE,str(2012), source)
        catalog.entries[name].add_quantity(
                FASTSTARS.ALIAS, row['ID'], source=source)
        radec = oname.strip('SDSSJ')
        radec = radec[0:2]+' '+radec[2:4]+' '+radec[4:9]+' '+radec[9:12]+' '+radec[12:14]+' '+radec[14:18]
        ra, dec = coord(radec, 
                unit=(u.hourangle, u.deg)).to_string(
                'hmsdms', sep=':', precision=8).split()
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
            FASTSTARS.SPECTRAL_TYPE, "F", source=source)
        catalog.entries[name].add_quantity(
            FASTSTARS.STELLAR_CLASS, "d", source=source)
    catalog.journal_entries()
    
    # 2012ApJ...751...55B
    datafile = os.path.join(catalog.get_current_task_repo(), 'ASCII',
                            'apj427101t1_ascii.csv')
    data = read(datafile, format='csv')
    for row in pbar(data, task_str):
        oname = str(row['Catalog']).replace(' ','')
        name, source = catalog.new_entry(oname, bibcode='2012ApJ...751...55B')
        if (FASTSTARS.DISCOVERER not in catalog.entries[name]):
            catalog.entries[name].add_quantity(FASTSTARS.DISCOVERER,'Warren R. Brown, Margaret J. Geller, Scott J. Kenyon', source)
            catalog.entries[name].add_quantity(FASTSTARS.DISCOVER_DATE,str(2012), source)
        if name == 'SDSSJ112359.47+751807.73':
            sourcewarren = catalog.entries[name].add_source(name='Warren R. Brown (private comm.)',secondary=True)
            catalog.entries[name].add_quantity(FASTSTARS.SPECTRAL_TYPE, 'B9', source=sourcewarren)
        gallon = float(str(row['Glon']))
        gallat = float(str(row['Glat']))
        if (oname!='US708') & (oname!='HE0437-5439'):
            radec = oname.strip('SDSSJ')
            radec = radec[0:2]+' '+radec[2:4]+' '+radec[4:9]+' '+radec[9:12]+' '+radec[12:14]+' '+radec[14:]
            ra, dec = coord(radec, 
                unit=(u.hourangle, u.deg)).to_string(
                'hmsdms', sep=':', precision=8).split()
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
            FASTSTARS.LUM_DIST, str(dhel), u_value='kpc', source=source, derived=True)
        sptype = str(row['Type']).split('/')
        for SPTYPE in sptype:
            if SPTYPE != "NA":
                catalog.entries[name].add_quantity(
                    FASTSTARS.SPECTRAL_TYPE, SPTYPE, source=source)
                catalog.entries[name].add_quantity(
                    FASTSTARS.STELLAR_CLASS, "d", source=source)
        if str(row['ID'])[:3]=='HVS':
            catalog.entries[name].add_quantity(
                FASTSTARS.ALIAS, row['ID'], source=source)
    catalog.journal_entries()
    
    # 2014ApJ...780....7P
    datafile = os.path.join(catalog.get_current_task_repo(), 'ASCII',
                            'apj485719t1_ascii.csv')
    data = read(datafile, format='csv')
    for row in pbar(data, task_str):
        oname = 'SDSS'+str(row['Catalog']).replace(' ','')
        name, source = catalog.new_entry(oname, bibcode='2014ApJ...780....7P')
        if (FASTSTARS.DISCOVERER not in catalog.entries[name]):
            catalog.entries[name].add_quantity(FASTSTARS.DISCOVERER,'Lauren E. Palladino, Katharine J. Schlesinger, Kelly Holley-Bockelmann, Carlos Allende Prieto, Timothy C. Beers, Young Sun Lee, Donald P. Schneider', source)
            catalog.entries[name].add_quantity(FASTSTARS.DISCOVER_DATE,str(2014), source)
        catalog.entries[name].add_quantity(
                FASTSTARS.ALIAS, 'Pal'+str(row['Pal']), source=source)
        radec = oname.strip('SDSSJ')
        radec = radec[0:2]+' '+radec[2:4]+' '+radec[4:9]+' '+radec[9:12]+' '+radec[12:14]+' '+radec[14:18]
        ra, dec = coord(radec, 
                unit=(u.hourangle, u.deg)).to_string(
                'hmsdms', sep=':', precision=8).split()
        catalog.entries[name].add_quantity(
            FASTSTARS.RA, ra, source=source)
        catalog.entries[name].add_quantity(
            FASTSTARS.DEC, dec, source=source)
        catalog.entries[name].add_quantity(
            FASTSTARS.VELOCITY, str(row['Vhel']), source=source)
        catalog.entries[name].add_quantity(
            FASTSTARS.LUM_DIST, str(row['Dhel']), u_value='kpc', source=source)
        catalog.entries[name].add_quantity(
            FASTSTARS.SPECTRAL_TYPE, "G", source=source)
        catalog.entries[name].add_quantity(
            FASTSTARS.SPECTRAL_TYPE, "K", source=source)
        catalog.entries[name].add_quantity(
            FASTSTARS.STELLAR_CLASS, "d", source=source)
    catalog.journal_entries()
    
    # 2014ApJ...787...89B
    datafile = os.path.join(catalog.get_current_task_repo(), 'ASCII',
                            'apj494602t1_ascii.csv')
    data = read(datafile, format='csv')
    for row in pbar(data, task_str):
        oname = str(row['Catalog']).replace(' ','')
        name, source = catalog.new_entry(oname, bibcode='2014ApJ...787...89B')
        if (FASTSTARS.DISCOVERER not in catalog.entries[name]):
            catalog.entries[name].add_quantity(FASTSTARS.DISCOVERER,'Warren R. Brown, Margaret J. Geller, Scott J. Kenyon', source)
            catalog.entries[name].add_quantity(FASTSTARS.DISCOVER_DATE,str(2014), source)
        radec = oname.strip('SDSSJ')
        radec = radec[0:2]+' '+radec[2:4]+' '+radec[4:9]+' '+radec[9:12]+' '+radec[12:14]+' '+radec[14:19]
        ra, dec = coord(radec, 
                unit=(u.hourangle, u.deg)).to_string(
                'hmsdms', sep=':', precision=8).split()
        catalog.entries[name].add_quantity(
            FASTSTARS.RA, ra, source=source)
        catalog.entries[name].add_quantity(
            FASTSTARS.DEC, dec, source=source)
        catalog.entries[name].add_quantity(
                FASTSTARS.VELOCITY, str(row['Vhel']), e_value=str(row['e_Vhel']), source=source)
        if FASTSTARS.SPECTRAL_TYPE not in catalog.entries[name]:
            sourceboubert = catalog.entries[name].add_source(name='Boubert et al. (2018, in prep)',secondary=True)
            catalog.entries[name].add_quantity(FASTSTARS.SPECTRAL_TYPE, 'A', source=sourceboubert)
            catalog.entries[name].add_quantity(FASTSTARS.SPECTRAL_TYPE, 'B', source=sourceboubert)
        galrad = float(str(row['RGC']))
        errgalrad = float(str(row['e_RGC']))
        dhel = rgc_to_dhel(galrad,gallon,gallat)
        dhel_lo = rgc_to_dhel(galrad-errgalrad,gallon,gallat)
        dhel_hi = rgc_to_dhel(galrad+errgalrad,gallon,gallat)
        catalog.entries[name].add_quantity(
            FASTSTARS.LUM_DIST, str(dhel), e_lower_value=dhel-dhel_lo, e_upper_value=dhel_hi-dhel, u_value='kpc', source=source, derived=True)
        if str(row['ID'])!='pBHVS':
            catalog.entries[name].add_quantity(
                FASTSTARS.ALIAS, 'HVS'+str(row['ID']), source=source)
        if str(row['ID'])=='22':
            catalog.entries[name].add_quantity(
                FASTSTARS.SPECTRAL_TYPE, "B", source=source)
            catalog.entries[name].add_quantity(
                FASTSTARS.STELLAR_CLASS, "d", source=source)
            catalog.entries[name].add_quantity(
                FASTSTARS.STELLAR_CLASS, "bhb", source=source)
            catalog.entries[name].add_quantity(
                FASTSTARS.STELLAR_CLASS, "g", source=source)
        elif str(row['ID'])=='23':
            catalog.entries[name].add_quantity(
                FASTSTARS.SPECTRAL_TYPE, "B", source=source)
            catalog.entries[name].add_quantity(
                FASTSTARS.STELLAR_CLASS, "d", source=source)
            catalog.entries[name].add_quantity(
                FASTSTARS.STELLAR_CLASS, "bhb", source=source)
            catalog.entries[name].add_quantity(
                FASTSTARS.STELLAR_CLASS, "g", source=source)
        elif str(row['ID'])=='24':
            catalog.entries[name].add_quantity(
                FASTSTARS.SPECTRAL_TYPE, "B", source=source)
            catalog.entries[name].add_quantity(
                FASTSTARS.STELLAR_CLASS, "d", source=source)
    catalog.journal_entries()
    
    # 2014EAS....67..255Z
    datafile = os.path.join(catalog.get_current_task_repo(), 'ASCII',
                            '1501.07824.csv')
    data = read(datafile, format='csv')
    for row in pbar(data, task_str):
        oname = str(row['ID'])
        name, source = catalog.new_entry(oname, bibcode='2014EAS....67..255Z')
        if (FASTSTARS.DISCOVERER not in catalog.entries[name]):
            catalog.entries[name].add_quantity(FASTSTARS.DISCOVERER, 'Y. Q. Zhang, M. C. Smith, J. L. Carlin', source)
            catalog.entries[name].add_quantity(FASTSTARS.DISCOVER_DATE,str(2014), source)
        catalog.entries[name].add_quantity(
                FASTSTARS.VELOCITY, str(row['Vhel']), e_value=str(row['e_Vhel']), source=source)
        catalog.entries[name].add_quantity(
            FASTSTARS.PROPER_MOTION_RA, str(row['pmra']), e_value=str(row['e_pmra']), u_value='mas/yr', source=source)
        catalog.entries[name].add_quantity(
            FASTSTARS.PROPER_MOTION_DEC, str(row['pmdec']), e_value=str(row['e_pmdec']), u_value='mas/yr', source=source)
        catalog.entries[name].add_quantity(
            FASTSTARS.LUM_DIST, str(row['Dhel']), u_value='kpc', source=source) # This distance may have some metallicity dependent uncertainty?
        catalog.entries[name].add_quantity(
                FASTSTARS.SPECTRAL_TYPE, 'K', source=source)
        catalog.entries[name].add_quantity(
                FASTSTARS.SPECTRAL_TYPE, 'M', source=source)
        catalog.entries[name].add_quantity(
                FASTSTARS.STELLAR_CLASS, 'd', source=source)
    catalog.journal_entries()
    
    # 2014ApJ...789L...2Z
    datafile = os.path.join(catalog.get_current_task_repo(), 'ASCII',
                            'apjl496832t1t2t3_ascii.csv')
    data = read(datafile, format='csv')
    for row in pbar(data, task_str):
        oname = str(row['ID'])
        name, source = catalog.new_entry(oname, bibcode='2014ApJ...789L...2Z')
        if (FASTSTARS.DISCOVERER not in catalog.entries[name]):
            catalog.entries[name].add_quantity(FASTSTARS.DISCOVERER,'Jing Zhong, Li Chen, Chao Liu, Richard de Grijs, Jinliang Hou, Shiyin Shen, Zhengyi Shao, Jing Li, Ali Luo, Jianrong Shi, Haotong Zhang, Ming Yang, Licai Deng, Ge Jin, Yong Zhang, Yonghui Hou, Zhenchao Zhang', source)
            catalog.entries[name].add_quantity(FASTSTARS.DISCOVER_DATE,str(2014), source)
        catalog.entries[name].add_quantity(
            FASTSTARS.ALIAS, str(row['Catalog']), source=source)
        radec = str(row['Catalog']).strip('J')
        radec = radec[0:2]+' '+radec[2:4]+' '+radec[4:9]+' '+radec[9:12]+' '+radec[12:14]+' '+radec[14:18]
        ra, dec = coord(radec, 
                unit=(u.hourangle, u.deg)).to_string(
                'hmsdms', sep=':', precision=8).split()
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
                FASTSTARS.SPECTRAL_TYPE, str(row['Type']), source=source)
        catalog.entries[name].add_quantity(
                FASTSTARS.STELLAR_CLASS, 'd', source=source)
    catalog.journal_entries()
    
    # 2014ApJ...794..146T
    datafile = os.path.join(catalog.get_current_task_repo(), 'ASCII',
                            'ApJ794146.csv')
    data = read(datafile, format='csv')
    for row in pbar(data, task_str):
        oname = str(row['Catalog']).replace(' ','')
        name, source = catalog.new_entry(oname, bibcode='2014ApJ...794..146T')
        if (FASTSTARS.DISCOVERER not in catalog.entries[name]):
            catalog.entries[name].add_quantity(FASTSTARS.DISCOVERER,'Christopher A. Theissen, Andrew A. West', source)
            catalog.entries[name].add_quantity(FASTSTARS.DISCOVER_DATE,str(2014), source)
        radec = oname.strip('SDSSJ')
        radec = radec[0:2]+' '+radec[2:4]+' '+radec[4:9]+' '+radec[9:12]+' '+radec[12:14]+' '+radec[14:18]
        ra, dec = coord(radec, 
                unit=(u.hourangle, u.deg)).to_string(
                'hmsdms', sep=':', precision=8).split()
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
                FASTSTARS.SPECTRAL_TYPE, str(row['Type']), source=source)
        catalog.entries[name].add_quantity(
                FASTSTARS.STELLAR_CLASS, 'd', source=source)
    catalog.journal_entries()
    
    # 2014ApJ...794..145S
    datafile = os.path.join(catalog.get_current_task_repo(), 'ASCII',
                            'apj501503t8_ascii.csv')
    data = read(datafile, format='csv')
    for row in pbar(data, task_str):
        oname = str(row['Catalog'])
        oname = oname[:4]+'J'+oname[4:]
        name, source = catalog.new_entry(oname, bibcode='2014ApJ...794..145S')
        if (FASTSTARS.DISCOVERER not in catalog.entries[name]):
            catalog.entries[name].add_quantity(FASTSTARS.DISCOVERER,'Antonia S. Savcheva, Andrew A. West, John J. Bochanski', source)
            catalog.entries[name].add_quantity(FASTSTARS.DISCOVER_DATE,str(2014), source)
        radec = oname.strip('SDSSJ')
        radec = radec[0:2]+' '+radec[2:4]+' '+radec[4:8]+' '+radec[8:11]+' '+radec[11:13]+' '+radec[13:17]
        ra, dec = coord(radec, 
                unit=(u.hourangle, u.deg)).to_string(
                'hmsdms', sep=':', precision=8).split()
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
        sptype = str(row['Type'])
        catalog.entries[name].add_quantity(
                FASTSTARS.SPECTRAL_TYPE, sptype[2:], source=source)
        catalog.entries[name].add_quantity(
                FASTSTARS.STELLAR_CLASS, "sd", source=source)
    catalog.journal_entries()
    
    # 2015MNRAS.447.2046H
    datafile = os.path.join(catalog.get_current_task_repo(), 'ASCII',
                            'hawkins2015.csv')
    data = read(datafile, format='csv')
    for row in pbar(data, task_str):
        oname = str(row['Catalog']).strip(' ')
        oname = 'RAVE'+oname
        name, source = catalog.new_entry(oname, bibcode='2015MNRAS.447.2046H')
        if (FASTSTARS.DISCOVERER not in catalog.entries[name]):
            catalog.entries[name].add_quantity(FASTSTARS.DISCOVERER,'K. Hawkins, G. Kordopatis, G. Gilmore, T. Masseron, R. F. G. Wyse, G. Ruchti, O. Bienaymé, J. Bland-Hawthorn, C. Boeche, K. Freeman, B. K. Gibson, E. K. Grebel, A. Helmi, A. Kunder, U. Munari, J. F. Navarro, Q. A. Parker, W. A. Reid, R. D. Scholz, G. Seabroke, A. Siebert, M. Steinmetz, F. Watson, T. Zwitter', source)
            catalog.entries[name].add_quantity(FASTSTARS.DISCOVER_DATE,str(2015), source)
        radec = oname.strip('RAVEJ')
        radec = radec[0:2]+' '+radec[2:4]+' '+radec[4:8]+' '+radec[8:11]+' '+radec[11:13]+' '+radec[13:15]
        ra, dec = coord(radec, 
                unit=(u.hourangle, u.deg)).to_string(
                'hmsdms', sep=':', precision=8).split()
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
            FASTSTARS.LUM_DIST, str(float(row['Dhel'])/1e3).strip(' '), e_value=str(float(row['e_Dhel'])/1e3).strip(' '), u_value='kpc', source=source)
        catalog.entries[name].add_quantity(
                FASTSTARS.STELLAR_CLASS, "g", source=source)
    catalog.journal_entries()
    
    # 2015A&A...576L..14Z
    datafile = os.path.join(catalog.get_current_task_repo(), 'ASCII',
                            'zeigerer.csv')
    data = read(datafile, format='csv')
    for row in pbar(data, task_str):
        oname = 'Pal'+str(row['Pal'])
        name, source = catalog.new_entry(oname, bibcode='2015A&A...576L..14Z')
        if (FASTSTARS.DISCOVERER not in catalog.entries[name]):
            catalog.entries[name].add_quantity(FASTSTARS.DISCOVERER,'E. Ziegerer, M. Volkert, U. Heber, A. Irrgang, B. T. Gänsicke, S. Geier', source)
            catalog.entries[name].add_quantity(FASTSTARS.DISCOVER_DATE,str(2015), source)
        catalog.entries[name].add_quantity(
            FASTSTARS.PROPER_MOTION_RA, str(row['pmra']), e_value=str(row['e_pmra']), u_value='mas/yr', source=source)
        catalog.entries[name].add_quantity(
            FASTSTARS.PROPER_MOTION_DEC, str(row['pmdec']), e_value=str(row['e_pmdec']), u_value='mas/yr', source=source)
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
        if (FASTSTARS.DISCOVERER not in catalog.entries[name]):
            catalog.entries[name].add_quantity(FASTSTARS.DISCOVERER,'Warren R. Brown, Jay Anderson, Oleg Y. Gnedin, Howard E. Bond, Margaret J. Geller, Scott J. Kenyon', source)
            catalog.entries[name].add_quantity(FASTSTARS.DISCOVER_DATE,str(2015), source)
        catalog.entries[name].add_quantity(
            FASTSTARS.ALIAS, str(row['ID']), source=source)
        catalog.entries[name].add_quantity(
            FASTSTARS.PROPER_MOTION_RA, str(row['pmra']).strip(' '), e_value=str(row['e_pmra']).strip(' '), u_value='mas/yr', source=source)
        catalog.entries[name].add_quantity(
            FASTSTARS.PROPER_MOTION_DEC, str(row['pmdec']).strip(' '), e_value=str(row['e_pmdec']).strip(' '), u_value='mas/yr', source=source)
        catalog.entries[name].add_quantity(
            FASTSTARS.LUM_DIST, str(row['Dhel']).strip(' '), e_value=str(row['e_Dhel']).strip(' '), u_value='kpc', source=source)
        if str(row['newspec']) == 'y':
            radec = oname.strip('SDSSJ')
            ra, dec = radec[:11], radec[11:]
            catalog.entries[name].add_quantity(
                FASTSTARS.RA, ra, source=source)
            catalog.entries[name].add_quantity(
                FASTSTARS.DEC, dec, source=source)
            catalog.entries[name].add_quantity(
                FASTSTARS.VELOCITY, str(row['Vhel']), e_value=str(row['e_Vhel']), source=source)
            #catalog.entries[name].add_quantity(
            #    FASTSTARS.LUM_DIST, str(row['Dhel']), e_value=str(row['e_Dhel']), u_value='kpc', source=source)
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
        if (FASTSTARS.DISCOVERER not in catalog.entries[name]):
            catalog.entries[name].add_quantity(FASTSTARS.DISCOVERER,'Yin-Bi Li, A.-Li Luo, Gang Zhao, You-Jun Lu, Peng Wei, Bing Du, Xiang Li, Yong-Heng Zhao, Zhan-Wen Han, Bo Wang, Yue Wu, Yong Zhang, Yong-Hui Hou, Yue-Fei Wang, Ming Yang', source)
            catalog.entries[name].add_quantity(FASTSTARS.DISCOVER_DATE,str(2015), source)
        radec = oname.strip('LAMOSTJ')
        radec = radec[0:2]+' '+radec[2:4]+' '+radec[4:9]+' '+radec[9:12]+' '+radec[12:14]+' '+radec[14:]
        ra, dec = coord(radec, 
                unit=(u.hourangle, u.deg)).to_string(
                'hmsdms', sep=':', precision=8).split()
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
                FASTSTARS.SPECTRAL_TYPE, 'F', source=source)
        catalog.entries[name].add_quantity(
                FASTSTARS.SPECTRAL_TYPE, 'G', source=source)
        catalog.entries[name].add_quantity(
                FASTSTARS.SPECTRAL_TYPE, 'K', source=source)
        catalog.entries[name].add_quantity(
                FASTSTARS.STELLAR_CLASS, 'd', source=source)
    catalog.journal_entries()
    
    # 2015AJ....150...77V
    datafile = os.path.join(catalog.get_current_task_repo(), 'ASCII',
                            'vickers2015.csv')
    data = read(datafile, format='csv')
    for row in pbar(data, task_str):
        oname = str(row['Catalog'])
        oname = 'SDSS'+oname
        name, source = catalog.new_entry(oname, bibcode='2015AJ....150...77V')
        if (FASTSTARS.DISCOVERER not in catalog.entries[name]):
            catalog.entries[name].add_quantity(FASTSTARS.DISCOVERER,'John J.Vickers, Martin C. Smith, Eva K. Grebel', source)
            catalog.entries[name].add_quantity(FASTSTARS.DISCOVER_DATE,str(2015), source)
        radec = oname.strip('SDSSJ')
        radec = radec[0:2]+' '+radec[2:4]+' '+radec[4:9]+' '+radec[9:12]+' '+radec[12:14]+' '+radec[14:]
        ra, dec = coord(radec, 
                unit=(u.hourangle, u.deg)).to_string(
                'hmsdms', sep=':', precision=8).split()
        catalog.entries[name].add_quantity(
            FASTSTARS.RA, ra, source=source)
        catalog.entries[name].add_quantity(
            FASTSTARS.DEC, dec, source=source)
        catalog.entries[name].add_quantity(
                FASTSTARS.SPECTRAL_TYPE, 'F', source=source)
        catalog.entries[name].add_quantity(
                FASTSTARS.SPECTRAL_TYPE, 'G', source=source)
        catalog.entries[name].add_quantity(
                FASTSTARS.SPECTRAL_TYPE, 'K', source=source)
        catalog.entries[name].add_quantity(
                FASTSTARS.SPECTRAL_TYPE, 'M', source=source)
        catalog.entries[name].add_quantity(
                FASTSTARS.STELLAR_CLASS, 'd', source=source)
    catalog.journal_entries()
    
    # 2015ApJ...813...26F
    datafile = os.path.join(catalog.get_current_task_repo(), 'ASCII',
                            'favia2015.csv')
    data = read(datafile, format='csv')
    for row in pbar(data, task_str):
        oname = str(row['Catalog'])
        oname = 'SDSS'+oname
        name, source = catalog.new_entry(oname, bibcode='2015ApJ...813...26F')
        if (FASTSTARS.DISCOVERER not in catalog.entries[name]):
            catalog.entries[name].add_quantity(FASTSTARS.DISCOVERER,'Andrej Favia, Andrew A. West, Christopher A. Theissen', source)
            catalog.entries[name].add_quantity(FASTSTARS.DISCOVER_DATE,str(2015), source)
        radec = oname.strip('SDSSJ')
        radec = radec[0:2]+' '+radec[2:4]+' '+radec[4:9]+' '+radec[9:12]+' '+radec[12:14]+' '+radec[14:]
        ra, dec = coord(radec, 
                unit=(u.hourangle, u.deg)).to_string(
                'hmsdms', sep=':', precision=8).split()
        catalog.entries[name].add_quantity(
            FASTSTARS.ALIAS, 'RdM'+str(row['ID']), source=source)
        catalog.entries[name].add_quantity(
            FASTSTARS.RA, ra, source=source)
        catalog.entries[name].add_quantity(
            FASTSTARS.DEC, dec, source=source)
        catalog.entries[name].add_quantity(
                FASTSTARS.VELOCITY, str(row['Vhel']), e_value=str(row['e_Vhel']), source=source)
        catalog.entries[name].add_quantity(
            FASTSTARS.PROPER_MOTION_RA, str(row['pmra1']), e_value=str(row['e_pmra1']), u_value='mas/yr', source=source)
        catalog.entries[name].add_quantity(
            FASTSTARS.PROPER_MOTION_DEC, str(row['pmdec1']), e_value=str(row['e_pmdec1']), u_value='mas/yr', source=source)
        catalog.entries[name].add_quantity(
            FASTSTARS.LUM_DIST, str(float(row['Dhel'])/1e3), e_value=str(float(row['e_Dhel'])/1e3), u_value='kpc', source=source)
        catalog.entries[name].add_quantity(
                FASTSTARS.SPECTRAL_TYPE, 'M'+str(row['Type']), source=source)
        catalog.entries[name].add_quantity(
                FASTSTARS.STELLAR_CLASS, 'd', source=source)
    catalog.journal_entries()
    
    # 2017ApJ...847L...9H
    datafile = os.path.join(catalog.get_current_task_repo(), 'ASCII',
                            'huang2017.csv')
    data = read(datafile, format='csv')
    for row in pbar(data, task_str):
        oname = str(row['Catalog'])
        name, source = catalog.new_entry(oname, bibcode='2017ApJ...847L...9H')
        if (FASTSTARS.DISCOVERER not in catalog.entries[name]):
            catalog.entries[name].add_quantity(FASTSTARS.DISCOVERER,'Y. Huang, X.-W. Liu, H.-W. Zhang, B.-Q. Chen, M.-S. Xiang, C. Wang, H.-B. Yuan, Z.-J. Tian, Y.-B. Li, B. Wang', source)
            catalog.entries[name].add_quantity(FASTSTARS.DISCOVER_DATE,str(2017), source)
        radec = oname.strip('J')
        radec = radec[0:2]+' '+radec[2:4]+' '+radec[4:9]+' '+radec[9:12]+' '+radec[12:14]+' '+radec[14:]
        ra, dec = coord(radec, 
                unit=(u.hourangle, u.deg)).to_string(
                'hmsdms', sep=':', precision=8).split()
        catalog.entries[name].add_quantity(
            FASTSTARS.ALIAS, str(row['Alias']), source=source)
        catalog.entries[name].add_quantity(
            FASTSTARS.RA, ra, source=source)
        catalog.entries[name].add_quantity(
            FASTSTARS.DEC, dec, source=source)
        catalog.entries[name].add_quantity(
                FASTSTARS.VELOCITY, str(row['Vhel']), e_value=str(row['e_Vhel']), source=source)
        catalog.entries[name].add_quantity(
            FASTSTARS.LUM_DIST, str(row['Dhel']), e_value=str(row['e_Dhel']), u_value='kpc', source=source)
        sptype = str(row['Type'])
        catalog.entries[name].add_quantity(
                FASTSTARS.SPECTRAL_TYPE, sptype[:2], source=source)
        lumclass = sptype[2:].split('/')
        for LC in lumclass:
            if LC == "IV":
                catalog.entries[name].add_quantity(
                    FASTSTARS.STELLAR_CLASS, 'sg', source=source)
            elif LC == "V":
                catalog.entries[name].add_quantity(
                    FASTSTARS.STELLAR_CLASS, 'd', source=source)
                    
    # 2017MNRAS.470.1388M
    datafile = os.path.join(catalog.get_current_task_repo(), 'ASCII',
                            'marchetti2017.csv')
    data = read(datafile, format='csv')
    for row in pbar(data, task_str):
        oname = 'TYC '+str(row['Tycho 2 ID']).strip(' ')
        name, source = catalog.new_entry(oname, bibcode='2017MNRAS.470.1388M')
        if (FASTSTARS.DISCOVERER not in catalog.entries[name]):
            catalog.entries[name].add_quantity(FASTSTARS.DISCOVERER,'T. Marchetti, E. M. Rossi, G. Kordopatis, A. G. A. Brown, A. Rimoldi, E. Starkenburg, K. Youakim, R. Ashley', source)
            catalog.entries[name].add_quantity(FASTSTARS.DISCOVER_DATE,str(2017), source)
        sourcegaia = catalog.entries[name].add_source(bibcode='2016A&A...595A...2G')
        radec = str(row['RADEC'])
        ra, dec = coord(radec, 
                unit=(u.hourangle, u.deg)).to_string(
                'hmsdms', sep=':', precision=8).split()
        catalog.entries[name].add_quantity(
            FASTSTARS.RA, ra, source=sourcegaia)
        catalog.entries[name].add_quantity(
            FASTSTARS.DEC, dec, source=sourcegaia)
        catalog.entries[name].add_quantity(
                FASTSTARS.VELOCITY, str(row['HRV']), e_value=str(row['e_HRV']), source=source)
        #catalog.entries[name].add_quantity(
        #    FASTSTARS.LUM_DIST, str(float(str(row['d']))/1e3), e_lower_value=str(float(str(row['e_low_d']))/1e3), e_upper_value=str(float(str(row['e_upp_d']))/1e3), u_value='kpc', source=source)
        
        if str(row['dspec'])!='--':
            catalog.entries[name].add_quantity(
                FASTSTARS.LUM_DIST, str(float(str(row['dspec']))/1e3), e_value=str(float(str(row['e_dspec']))/1e3), u_value='kpc', source=source)
        sptype = str(row['SpectralType'])
        if sptype != '--':
            catalog.entries[name].add_quantity(
                FASTSTARS.SPECTRAL_TYPE, sptype, source=source)
        stellarclass = str(row['StellarClass']).split('/')
        if str(row['StellarClass'])!='--':
            for SC in stellarclass:
                catalog.entries[name].add_quantity(
                    FASTSTARS.STELLAR_CLASS, SC, source=source)
                    
    # 2018arXiv180410607M
    datafile = os.path.join(catalog.get_current_task_repo(), 'ASCII',
                            'marchetti2018.txt')
    data = read(datafile)
    for row in pbar(data, task_str):
        oname = 'GaiaDR2 '+str(row['source_id']).strip(' ')
        name, source = catalog.new_entry(oname, bibcode='2018arXiv180410607M')
        if (FASTSTARS.DISCOVERER not in catalog.entries[name]):
            catalog.entries[name].add_quantity(FASTSTARS.DISCOVERER,'T. Marchetti, E. M. Rossi, A. G. A. Brown', source)
            catalog.entries[name].add_quantity(FASTSTARS.DISCOVER_DATE,str(2018), source)
        #sourcegaia = catalog.entries[name].add_source(bibcode='2016A&A...595A...2G')
        ra, dec = coord(ra=float(row['ra'])*u.deg,dec=float(row['dec'])*u.deg,frame='icrs').to_string(
                'hmsdms', sep=':', precision=8).split()
        catalog.entries[name].add_quantity(
            FASTSTARS.RA, ra, source=source)
        catalog.entries[name].add_quantity(
            FASTSTARS.DEC, dec, source=source)
        #catalog.entries[name].add_quantity(
        #        FASTSTARS.VELOCITY, str(row['HRV']), e_value=str(row['e_HRV']), source=source)
        #catalog.entries[name].add_quantity(
        #    FASTSTARS.LUM_DIST, str(float(str(row['d']))/1e3), e_lower_value=str(float(str(row['e_low_d']))/1e3), e_upper_value=str(float(str(row['e_upp_d']))/1e3), u_value='kpc', source=source)
        
        #if str(row['dspec'])!='--':
        #    catalog.entries[name].add_quantity(
        #        FASTSTARS.LUM_DIST, str(float(str(row['dspec']))/1e3), e_value=str(float(str(row['e_dspec']))/1e3), u_value='kpc', source=source)
        #sptype = str(row['SpectralType'])
        #if sptype != '--':
        #    catalog.entries[name].add_quantity(
        #        FASTSTARS.SPECTRAL_TYPE, sptype, source=source)
        #stellarclass = str(row['StellarClass']).split('/')
        #if str(row['StellarClass'])!='--':
        #    for SC in stellarclass:
        #        catalog.entries[name].add_quantity(
        #            FASTSTARS.STELLAR_CLASS, SC, source=source)
    catalog.journal_entries()
    return
