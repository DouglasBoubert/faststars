"""Cleanup catalog before final write to disk."""
import re
import statistics
import warnings
from decimal import Decimal
from math import log10, pi, sqrt

from astrocats.catalog.quantity import QUANTITY
from astrocats.catalog.utils import (get_sig_digits, is_number, pbar,
                                     pretty_num, tprint, uniq_cdl)
from astropy import units as un
from astropy.coordinates import SkyCoord as coord
from astropy.cosmology import Planck15 as cosmo
from astropy.cosmology import z_at_value

from ..constants import CLIGHT, KM
from ..faststars import FASTSTARS


def do_cleanup(catalog):
    """Cleanup catalog after importing all data."""
    task_str = catalog.get_current_task_str()

    # Set preferred names, calculate some columns based on imported data,
    # sanitize some fields
    keys = list(catalog.entries.keys())

    cleanupcnt = 0
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
        name = catalog.entries[name].set_preferred_name()

        aliases = catalog.entries[name].get_aliases()

        if (FASTSTARS.RA not in catalog.entries[name] or
                FASTSTARS.DEC not in catalog.entries[name]):
            prefixes = [
                'SDSS'
            ]
            for alias in aliases:
                for prefix in prefixes:
                    if (alias.startswith(prefix) and
                            is_number(alias.replace(prefix, '')[:6])):
                        noprefix = alias.split(':')[-1].replace(
                            prefix, '').replace('.', '')
                        decsign = '+' if '+' in noprefix else '-'
                        noprefix = noprefix.replace('+', '|').replace('-', '|')
                        nops = noprefix.split('|')
                        if len(nops) < 2:
                            continue
                        rastr = nops[0]
                        decstr = nops[1]
                        ra = ':'.join([rastr[:2], rastr[2:4], rastr[4:6]]) + \
                            ('.' + rastr[6:] if len(rastr) > 6 else '')
                        dec = (
                            decsign + ':'.join(
                                [decstr[:2], decstr[2:4], decstr[4:6]]) +
                            ('.' + decstr[6:] if len(decstr) > 6 else ''))
                        if catalog.args.verbose:
                            tprint('Added ra/dec from name: ' + ra + ' ' + dec)
                        source = catalog.entries[name].add_self_source()
                        catalog.entries[name].add_quantity(
                            FASTSTARS.RA, ra, source, derived=True)
                        catalog.entries[name].add_quantity(
                            FASTSTARS.DEC, dec, source, derived=True)
                        break
                if FASTSTARS.RA in catalog.entries[name]:
                    break

        catalog.entries[name].sanitize()
        catalog.journal_entries(bury=True, final=True, gz=True)
        cleanupcnt = cleanupcnt + 1
        if catalog.args.travis and cleanupcnt % 1000 == 0:
            break

    catalog.save_caches()

    return
