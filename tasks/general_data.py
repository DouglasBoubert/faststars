# -*- coding: utf-8 -*-
"""General data import tasks."""
import json
import os
from collections import OrderedDict
from copy import deepcopy
from decimal import Decimal
from glob import glob

from astrocats.catalog.photometry import PHOTOMETRY
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
