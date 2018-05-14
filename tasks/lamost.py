"""Import tasks for LAMOST."""
import os
from decimal import Decimal

import astropy.coordinates as coord
import astropy.units as un
import astropy.constants as con
import requests
from astrocats.catalog.spectrum import SPECTRUM
from astrocats.catalog.photometry import PHOTOMETRY
from astrocats.catalog.utils import jd_to_mjd, pbar, is_number, pretty_num
from astropy.io import fits
from astropy.time import Time as astrotime
from astroquery.vizier import Vizier

from ..faststars import FASTSTARS


def do_lamost(catalog):
    """Import spectra from LAMOST."""
    task_str = catalog.get_current_task_str()

    # Set preferred names, calculate some columns based on imported data,
    # sanitize some fields
    keys = list(catalog.entries.keys())

    viz = Vizier(columns=["**"])

    fureps = {'erg/cm2/s/A': 'erg/s/cm^2/Angstrom'}

    c_kms = con.c.cgs.value / 1.0e5

    for oname in pbar(keys, task_str):
        # Some events may be merged in cleanup process, skip them if
        # non-existent.

        if (FASTSTARS.RA not in catalog.entries[oname] or
                FASTSTARS.DEC not in catalog.entries[oname]):
            continue
        else:
            result = viz.query_region(
                coord.SkyCoord(
                    ra=catalog.entries[oname][FASTSTARS.RA][0]['value'],
                    dec=catalog.entries[oname][FASTSTARS.DEC][0]['value'],
                    unit=(un.hourangle, un.deg), frame='icrs'),
                width="2s",
                catalog="V/149/dr2")

            if not result.keys():
                continue
            tab = result['V/149/dr2']

            star = None
            for row in tab:
                if row['Class'].lower() == 'star':
                    star = row
                    break
            if not star:
                continue

            try:
                name, source = catalog.new_entry(
                    oname, bibcode='2016yCat.5149....0L', srcname='LAMOST',
                    url='http://dr3.lamost.org/')
            except Exception:
                catalog.log.warning(
                    '"{}" was not found, suggests merge occurred in cleanup '
                    'process.'.format(oname))
                continue

            if row['SubClass'] is not 'Non':
                catalog.entries[name].add_quantity(
                    FASTSTARS.SPECTRAL_TYPE, row['SubClass'], source=source)

            if is_number(row['z']):
                catalog.entries[name].add_quantity(
                    FASTSTARS.REDSHIFT, str(row['z']), e_value=str(row['e_z']),
                    source=source)
                catalog.entries[name].add_quantity(
                    FASTSTARS.VELOCITY,
                    pretty_num(float(row['z']) * c_kms, sig=5),
                    e_value=pretty_num(float(row['e_z'] * c_kms), sig=5),
                    source=source)

            mag_types = list(row['magType'].replace('psf_', ''))

            nmt = []
            nmi = 0
            for mt in mag_types:
                if is_number(mt):
                    nmt[nmi - 1] += mt
                else:
                    nmt += mt
                    nmi += 1
            mag_types = [x.upper() if x in [
                'b', 'v', 'j', 'h'] else x for x in nmt]

            for mi, mt in enumerate(mag_types):
                snrf = 'snr' + mt.lower()
                if snrf in row.columns and float(row[snrf]) < 3:
                    continue
                photodict = {
                    PHOTOMETRY.TIME: str(row['MJD']),
                    PHOTOMETRY.U_TIME: 'MJD',
                    PHOTOMETRY.BAND: mt,
                    PHOTOMETRY.TELESCOPE: 'LAMOST',
                    PHOTOMETRY.MAGNITUDE: str(row['mag' + str(mi + 1)]),
                    PHOTOMETRY.SOURCE: source
                }
                if snrf in row.columns:
                    photodict[PHOTOMETRY.E_MAGNITUDE] = str(
                        Decimal('2.5') *
                        (Decimal('1') + Decimal('1') /
                         Decimal(str(row[snrf]))).log10())[:5]
                catalog.entries[name].add_photometry(**photodict)

            vname = row['PlanId']

            ffile = ('spec-' +
                     row['LMJD'] + '-' + vname + '_sp' + row['spId'] + '-' +
                     row['FiberId'] + '.fits.gz')

            furl = 'http://dr3.lamost.org/sas/fits/' + vname + '/' + ffile

            datafile = os.path.join(
                catalog.get_current_task_repo(), 'LAMOST', ffile)

            if not os.path.exists(datafile):
                fr = requests.get(furl)

                open(datafile, 'wb').write(fr.content)

            hdulist = fits.open(datafile)
            for oi, obj in enumerate(hdulist[0].header):
                if any(x in ['.', '/'] for x in obj):
                    del (hdulist[0].header[oi])
            hdulist[0].verify('silentfix')
            hdrkeys = list(hdulist[0].header.keys())
            # print(hdrkeys)
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
                if hdulist[0].header['NAXIS'] == 2:
                    waves = [str(x) for x in list(hdulist[0].data)[2]]
                    fluxes = [str(x) for x in list(hdulist[0].data)[1]]
                else:
                    print('Warning: Skipping FITS spectrum `{}`.'.format(
                        datafile))
                    continue
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
                SPECTRUM.FILENAME: ffile,
                SPECTRUM.SOURCE: source
            }
            if 'TELESCOP' in hdrkeys:
                specdict[SPECTRUM.TELESCOPE] = hdulist[0].header['TELESCOP']
            if 'INSTRUME' in hdrkeys:
                specdict[SPECTRUM.INSTRUMENT] = hdulist[0].header['INSTRUME']
            if 'SITENAME' in hdrkeys:
                specdict[SPECTRUM.OBSERVATORY] = hdulist[0].header['SITENAME']
            elif 'OBSERVAT' in hdrkeys:
                specdict[SPECTRUM.OBSERVATORY] = hdulist[0].header['OBSERVAT']
            if 'OBSERVER' in hdrkeys:
                specdict[SPECTRUM.OBSERVER] = hdulist[0].header['OBSERVER']
            if 'AIRMASS' in hdrkeys:
                specdict[SPECTRUM.AIRMASS] = hdulist[0].header['AIRMASS']
            catalog.entries[name].add_spectrum(**specdict)
            hdulist.close()
            catalog.journal_entries()

    return
