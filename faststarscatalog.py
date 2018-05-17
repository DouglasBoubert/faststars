"""FastStars specific catalog class."""
import codecs
import json
import os
from collections import OrderedDict
from datetime import datetime
from subprocess import check_output

from astrocats.catalog.catalog import Catalog
from astrocats.catalog.quantity import QUANTITY
from astrocats.catalog.utils import read_json_arr, read_json_dict

from .faststars import FASTSTARS, FastStars
from .utils import name_clean


class FastStarsCatalog(Catalog):
    """Catalog class for `FastStars` objects."""

    class PATHS(Catalog.PATHS):
        """Paths to catalog inputs/outputs."""

        PATH_BASE = os.path.abspath(os.path.dirname(__file__))

        def __init__(self, catalog):
            """Initialize paths."""
            super(FastStarsCatalog.PATHS, self).__init__(catalog)
            # auxiliary datafiles
            self.TYPE_SYNONYMS = os.path.join(
                self.PATH_INPUT, 'type-synonyms.json')
            self.SOURCE_SYNONYMS = os.path.join(
                self.PATH_INPUT, 'source-synonyms.json')
            self.URL_REDIRECTS = os.path.join(
                self.PATH_INPUT, 'url-redirects.json')
            self.BIBERRORS = os.path.join(self.PATH_INPUT, 'biberrors.json')
            # cached datafiles
            self.BIBAUTHORS = os.path.join(
                self.PATH_OUTPUT, 'cache', 'bibauthors.json')
            self.EXTINCT = os.path.join(
                self.PATH_OUTPUT, 'cache', 'extinctions.json')

        def get_repo_years(self):
            """Return an array of years based upon output repositories."""
            repo_folders = self.get_repo_output_folders(bones=False)
            repo_years = [int(repo_folders[x][-4:])
                          for x in range(len(repo_folders))]
            repo_years[0] -= 1
            return repo_years

    class SCHEMA(object):
        """Define the HASH/URL associated with the present schema."""

        HASH = (check_output(['git', '-C', 'astrocats/faststars',
                              'log', '-n', '1', '--format="%h"',
                              '--', 'SCHEMA.md'])
                .decode('ascii').strip().strip('"').strip())
        URL = ('https://github.com/astrocatalogs/faststars/blob/' + HASH +
               '/SCHEMA.md')

    def __init__(self, args, log):
        """Initialize catalog."""
        # Initialize super `astrocats.catalog.catalog.Catalog` object
        super(FastStarsCatalog, self).__init__(args, log)
        self.proto = FastStars
        self._load_aux_data()
        return

    def should_bury(self, name):
        """Determine whether a fast star should be "buried".

        For fast stars, objects that have enough data such that they can be
        definitively determined to be "bound" are buried.
        """
        bury_entry = False

        if (FASTSTARS.BOUND_PROBABILITY in self.entries[name] and
            not self.entries[name][
                FASTSTARS.BOUND_PROBABILITY][0].get(
                    QUANTITY.UPPER_LIMIT, False) and float(
                    self.entries[name][FASTSTARS.BOUND_PROBABILITY][0][
                        QUANTITY.VALUE]) >= 0.9998):
            bury_entry = True

        return (bury_entry, True)

    def _load_aux_data(self):
        """Load auxiliary dictionaries for use in this catalog."""
        # Create/Load auxiliary dictionaries
        self.nedd_dict = OrderedDict()
        self.bibauthor_dict = read_json_dict(self.PATHS.BIBAUTHORS)
        self.biberror_dict = read_json_dict(self.PATHS.BIBERRORS)
        self.extinctions_dict = read_json_dict(self.PATHS.EXTINCT)
        self.source_syns = read_json_dict(self.PATHS.SOURCE_SYNONYMS)
        self.url_redirs = read_json_dict(self.PATHS.URL_REDIRECTS)
        self.type_syns = read_json_dict(self.PATHS.TYPE_SYNONYMS)
        # Create/Load auxiliary arrays
        #self.nonsnetypes = read_json_arr(self.PATHS.NON_SNE_TYPES)
        return

    def save_caches(self):
        """Save caches to JSON files."""
        jsonstring = json.dumps(self.bibauthor_dict, indent='\t',
                                separators=(',', ':'), ensure_ascii=False)
        with codecs.open(self.PATHS.BIBAUTHORS, 'w', encoding='utf8') as f:
            f.write(jsonstring)
        jsonstring = json.dumps(self.extinctions_dict, indent='\t',
                                separators=(',', ':'), ensure_ascii=False)
        with codecs.open(self.PATHS.EXTINCT, 'w', encoding='utf8') as f:
            f.write(jsonstring)

    def clean_entry_name(self, name):
        """Clean entry's name."""
        return name_clean(name)
