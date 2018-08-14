# -*- coding: utf-8 -*-
"""General data import tasks."""
import json
import os
from copy import deepcopy
from decimal import Decimal
from glob import glob

from astrocats.catalog.spectrum import SPECTRUM
from astrocats.catalog.utils import jd_to_mjd, pbar_strings
from astropy.io import fits
from astropy.time import Time as astrotime

from ..faststars import FASTSTARS, FastStars


def do_internal(catalog):
    """Load events from files in the 'internal' repository, and save them."""
    task_str = catalog.get_current_task_str()
    path_pattern = os.path.join(catalog.get_current_task_repo(), '*.json')
    files = glob(path_pattern)
    catalog.log.debug("found {} files matching '{}'".format(
        len(files), path_pattern))
    for datafile in pbar_strings(files, task_str):
        new_entry = FastStars.init_from_file(
            catalog, path=datafile, clean=True, merge=True)

        name = new_entry[FASTSTARS.NAME]
        old_name = None

        for alias in new_entry.get_aliases():
            if catalog.entry_exists(alias):
                old_name = catalog.get_preferred_name(alias)
                if catalog.entries[old_name]._stub:
                    catalog.add_entry(old_name)
                break

        if old_name:
            old_entry = deepcopy(catalog.entries[old_name])
            catalog.copy_entry_to_entry(new_entry, old_entry)
            catalog.entries[old_name] = old_entry
        else:
            catalog.entries[name] = new_entry

        catalog.journal_entries()

    return


def do_external_fits_spectra(catalog):
    fpath = catalog.get_current_task_repo()
    with open(os.path.join(fpath, 'fits', 'meta.json'), 'r') as f:
        metadict = json.loads(f.read())

    fureps = {'erg/cm2/s/A': 'erg/s/cm^2/Angstrom'}
    task_str = catalog.get_current_task_str()
    path_pattern = os.path.join(
        catalog.get_current_task_repo(), 'fits', '*.fits')
    files = glob(path_pattern)
    for datafile in files:
        filename = datafile.split('/')[-1]
        if filename == 'meta.json':
            continue
        hdulist = fits.open(datafile)
        for oi, obj in enumerate(hdulist[0].header):
            if any(x in ['.', '/'] for x in obj):
                del (hdulist[0].header[oi])
        hdulist[0].verify('silentfix')
        hdrkeys = list(hdulist[0].header.keys())
        # print(hdrkeys)
        name = ''
        if filename in metadict:
            if 'name' in metadict[filename]:
                name = metadict[filename]['name']
        if not name:
            name = hdulist[0].header['OBJECT']
        if 'bibcode' in metadict[filename]:
            name, source = catalog.new_entry(
                name, bibcode=metadict[filename]['bibcode'])
        elif 'donator' in metadict[filename]:
            name, source = catalog.new_entry(
                name, srcname=metadict[filename]['donator'])
        else:
            if 'OBSERVER' in hdrkeys:
                name, source = catalog.new_entry(
                    name, srcname=hdulist[0].header['OBSERVER'])
            else:
                name = catalog.add_entry(name)
                source = catalog.entries[name].add_self_source()
        # for key in hdulist[0].header.keys():
        #     print(key, hdulist[0].header[key])
        if hdulist[0].header['SIMPLE']:
            if 'JD' in hdrkeys:
                mjd = str(jd_to_mjd(Decimal(str(hdulist[0].header['JD']))))
            elif 'MJD' in hdrkeys:
                mjd = str(hdulist[0].header['MJD'])
            elif 'DATE-OBS' in hdrkeys:
                if 'T' in hdulist[0].header['DATE-OBS']:
                    dateobs = hdulist[0].header['DATE-OBS'].strip()
                elif 'UTC-OBS' in hdrkeys:
                    dateobs = hdulist[0].header['DATE-OBS'].strip(
                    ) + 'T' + hdulist[0].header['UTC-OBS'].strip()
                mjd = str(astrotime(dateobs, format='isot').mjd)
            else:
                raise ValueError("Couldn't find JD/MJD for spectrum.")
            w0 = hdulist[0].header['CRVAL1']
            if hdulist[0].header['NAXIS'] == 1:
                wd = hdulist[0].header['CDELT1']
                fluxes = [str(x) for x in list(hdulist[0].data)]
                errors = False
            elif hdulist[0].header['NAXIS'] == 2:
                wd = hdulist[0].header['CD1_1']
                fluxes = [str(x) for x in list(hdulist[0].data)[0]]
                errors = False
            elif hdulist[0].header['NAXIS'] == 3:
                wd = hdulist[0].header['CD1_1']
                fluxes = [str(x) for x in list(hdulist[0].data)[0][0]]
                errors = [str(x) for x in list(hdulist[0].data)[3][0]]
            else:
                print('Warning: Skipping FITS spectrum `{}`.'.format(filename))
                continue
            waves = [str(w0 + wd * x) for x in range(0, len(fluxes))]
        else:
            raise ValueError('Non-simple FITS import not yet supported.')
        if 'BUNIT' in hdrkeys:
            fluxunit = hdulist[0].header['BUNIT']
            if fluxunit in fureps:
                fluxunit = fureps[fluxunit]
        else:
            if max([float(x) for x in fluxes]) < 1.0e-5:
                fluxunit = 'erg/s/cm^2/Angstrom'
            else:
                fluxunit = 'Uncalibrated'
        specdict = {
            SPECTRUM.U_WAVELENGTHS: 'Angstrom',
            SPECTRUM.WAVELENGTHS: waves,
            SPECTRUM.TIME: mjd,
            SPECTRUM.U_TIME: 'MJD',
            SPECTRUM.FLUXES: fluxes,
            SPECTRUM.U_FLUXES: fluxunit,
            SPECTRUM.FILENAME: filename,
            SPECTRUM.SOURCE: source
        }
        if 'TELESCOP' in hdrkeys:
            specdict[SPECTRUM.TELESCOPE] = hdulist[0].header['TELESCOP']
        if 'INSTRUME' in hdrkeys:
            specdict[SPECTRUM.INSTRUMENT] = hdulist[0].header['INSTRUME']
        if 'AIRMASS' in hdrkeys:
            specdict[SPECTRUM.AIRMASS] = hdulist[0].header['AIRMASS']
        if errors:
            specdict[SPECTRUM.ERRORS] = errors
            specdict[SPECTRUM.U_ERRORS] = fluxunit
        if 'SITENAME' in hdrkeys:
            specdict[SPECTRUM.OBSERVATORY] = hdulist[0].header['SITENAME']
        elif 'OBSERVAT' in hdrkeys:
            specdict[SPECTRUM.OBSERVATORY] = hdulist[0].header['OBSERVAT']
        if 'OBSERVER' in hdrkeys:
            specdict[SPECTRUM.OBSERVER] = hdulist[0].header['OBSERVER']
        catalog.entries[name].add_spectrum(**specdict)
        hdulist.close()
        catalog.journal_entries()
    return
