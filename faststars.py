"""FastStars class."""
import warnings
from collections import OrderedDict
from decimal import Decimal

import numpy as np
from astrocats.catalog.entry import ENTRY, Entry
from astrocats.catalog.key import KEY_TYPES, Key
from astrocats.catalog.photometry import PHOTOMETRY
from astrocats.catalog.quantity import QUANTITY
from astrocats.catalog.source import SOURCE
from astrocats.catalog.utils import (bib_priority, get_sig_digits,
                                     get_source_year, is_integer, is_number,
                                     jd_to_mjd, listify, make_date_string,
                                     pretty_num, uniq_cdl)
from astropy.time import Time as astrotime
from six import string_types

from .constants import MAX_VISUAL_BANDS
from .utils import frame_priority, host_clean, radec_clean


class FASTSTARS(ENTRY):
    """FastStars `Key` child class."""
    DISCOVERY_DATE = Key('discoverdate', KEY_TYPES.STRING)
    PROPER_MOTION_RA = Key('propermotionra', KEY_TYPES.NUMERIC,
        replace_better=True)
    PROPER_MOTION_DEC = Key('propermotiondec', KEY_TYPES.NUMERIC,
        replace_better=True)
    PARALLAX = Key('parallax', KEY_TYPES.NUMERIC)
    CLAIMED_TYPE = Key('claimedtype',
                       KEY_TYPES.STRING,
                       kind_preference=['HVS'],
                       replace_better=True)
    SPECTRAL_TYPE = Key('spectraltype',
                        KEY_TYPES.STRING,
                        replace_better=True)
    ERRORS = Key('errors')


class FastStars(Entry):
    """FastStars `Entry` child class.

    NOTE: OrderedDict data is just the `name` values from the JSON file.
          I.e. it does not include the highest nesting level
          { name: DATA }, it *just* includes DATA

    FIX: check that no stored values are empty/invalid (delete key in that
         case?)
    """

    _KEYS = FASTSTARS

    def __init__(self, catalog, name=None, stub=False):
        """Initialize `FastStars`."""
        super(FastStars, self).__init__(catalog, name, stub=stub)
        return

    def _append_additional_tags(self, name, sources, quantity):
        """Append additional bits of data to an existing quantity when a newly
        added quantity is found to be a duplicate
        """
        svalue = quantity.get(QUANTITY.VALUE, '')
        serror = quantity.get(QUANTITY.E_VALUE, '')
        sprob = quantity.get(QUANTITY.PROB, '')
        skind = quantity.get(QUANTITY.KIND, '')

        for ii, ct in enumerate(self[name]):
            if ct[QUANTITY.VALUE] == svalue and sources:
                if ct.get(QUANTITY.KIND, '') != skind:
                    return
                for source in sources.split(','):
                    if (source not in
                            self[name][ii][QUANTITY.SOURCE].split(',')):
                        self[name][ii][QUANTITY.SOURCE] += ',' + source
                        if serror and QUANTITY.E_VALUE not in self[name][ii]:
                            self[name][ii][QUANTITY.E_VALUE] = serror
                        if sprob and QUANTITY.PROB not in self[name][ii]:
                            self[name][ii][QUANTITY.PROB] = sprob
                return

    def _clean_quantity(self, quantity):
        """Clean quantity value before it is added to entry."""
        value = quantity.get(QUANTITY.VALUE, '').strip()
        error = quantity.get(QUANTITY.E_VALUE, '').strip()
        unit = quantity.get(QUANTITY.U_VALUE, '').strip()
        kinds = [x.strip() for x in listify(quantity.get(QUANTITY.KIND, []))]
        key = quantity._key

        if not value:
            return False

        if error and (not is_number(error) or float(error) < 0):
            raise ValueError(self[self._KEYS.NAME] + "'s quanta " + key +
                             ' error value must be a number and positive.')

        # Set default units
        if not unit and key == self._KEYS.VELOCITY:
            unit = 'km/s'
        if not unit and key == self._KEYS.RA:
            unit = 'hours'
        if not unit and key == self._KEYS.DEC:
            unit = 'degrees'
        if not unit and key in [self._KEYS.LUM_DIST, self._KEYS.COMOVING_DIST]:
            unit = 'Mpc'


        if is_number(value):
            value = '%g' % Decimal(value)
        if error:
            error = '%g' % Decimal(error)

        if value:
            quantity[QUANTITY.VALUE] = value
        if error:
            quantity[QUANTITY.E_VALUE] = error
        if unit:
            quantity[QUANTITY.U_VALUE] = unit
        if kinds:
            quantity[QUANTITY.KIND] = kinds if len(kinds) > 1 else kinds[0]
        elif QUANTITY.KIND in quantity:
            del (quantity[QUANTITY.KIND])

        return True

    def add_quantity(self,
                     quantities,
                     value,
                     source,
                     forcereplacebetter=False,
                     **kwargs):
        """Add `Quantity` to `FastStars`."""
        success = super(FastStars, self).add_quantity(
            quantities, value, source, **kwargs)

        if not success:
            return

        for quantity in listify(quantities):
            my_quantity_list = self.get(quantity, [])

            if ((forcereplacebetter or quantity.replace_better) and
                    len(my_quantity_list) > 1):

                # The quantity that was just added should be last in the list
                added_quantity = my_quantity_list.pop()

                newquantities = []
                isworse = True

                if type(quantity) != Key:
                    isworse = False
                elif quantity.type == KEY_TYPES.NUMERIC:
                    newsig = get_sig_digits(added_quantity[QUANTITY.VALUE])
                    for ct in my_quantity_list:
                        addct = False
                        checke = False
                        if (len(quantity.kind_preference) > 0 and not set(
                                listify(ct.get(QUANTITY.KIND, [])))
                                .isdisjoint(quantity.kind_preference) and
                                not set(
                                    listify(
                                        added_quantity.get(QUANTITY.KIND,
                                                            [])))
                                .isdisjoint(quantity.kind_preference)):
                            aqi = min([
                                quantity.kind_preference.index(x)
                                for x in listify(added_quantity[
                                    QUANTITY.KIND])
                            ])
                            qqi = min([
                                quantity.kind_preference.index(x)
                                for x in listify(ct[QUANTITY.KIND])
                            ])
                            if aqi > qqi:
                                addct = True
                            if aqi == qqi:
                                checke = True
                            if aqi <= qqi:
                                isworse = False
                        else:
                            checke = True
                        if checke and QUANTITY.E_VALUE in ct:
                            if QUANTITY.E_VALUE in added_quantity:
                                if (float(added_quantity[QUANTITY.E_VALUE])
                                        >= float(ct[QUANTITY.E_VALUE])):
                                    addct = True
                                if (float(added_quantity[QUANTITY.E_VALUE])
                                        <= float(ct[QUANTITY.E_VALUE])):
                                    isworse = False
                        else:
                            if (checke and
                                    QUANTITY.E_VALUE in added_quantity):
                                isworse = False
                            else:
                                oldsig = get_sig_digits(ct[QUANTITY.VALUE])
                                if oldsig >= newsig:
                                    addct = True
                                if newsig >= oldsig:
                                    isworse = False
                        if addct:
                            newquantities.append(ct)
                elif quantity.type == KEY_TYPES.STRING:
                    for ct in my_quantity_list:
                        addct = False
                        if (len(quantity.kind_preference) > 0 and not set(
                                listify(ct.get(QUANTITY.KIND, [])))
                                .isdisjoint(quantity.kind_preference) and
                                not set(
                                    listify(
                                        added_quantity.get(QUANTITY.KIND,
                                                            [])))
                                .isdisjoint(quantity.kind_preference)):
                            aqi = min([
                                quantity.kind_preference.index(x)
                                for x in listify(added_quantity[
                                    QUANTITY.KIND])
                            ])
                            qqi = min([
                                quantity.kind_preference.index(x)
                                for x in listify(ct[QUANTITY.KIND])
                            ])
                            if aqi >= qqi:
                                addct = True
                            if aqi <= qqi:
                                isworse = False
                        else:
                            addct = True
                            isworse = False
                        if addct:
                            newquantities.append(ct)

                if isworse:
                    self._log.info("Removing quantity '{}' with value '{}' "
                                   "and kind '{}' determined to be worse than "
                                   "existing alternative values.".format(
                                       quantity, added_quantity[
                                           QUANTITY.VALUE],
                                       added_quantity.get(QUANTITY.KIND, '')))
                else:
                    newquantities.append(added_quantity)
                if len(newquantities) > 0:
                    self[quantity] = newquantities

        return True

    def add_source(self, **kwargs):
        # Sanitize some fields before adding source
        # Replace reference names and URLs using dictionaries.
        if (kwargs.get(SOURCE.BIBCODE, []) and
                len(kwargs[SOURCE.BIBCODE]) != 19):
            raise ValueError("Bibcode '{}' must be exactly 19 characters "
                             "long".format(kwargs[SOURCE.BIBCODE]))

        # if SOURCE.NAME not in kwargs:
        #     kwargs[SOURCE.NAME] = kwargs[SOURCE.BIBCODE]

        if SOURCE.NAME in kwargs:

            for rep in self.catalog.source_syns:
                if kwargs[SOURCE.NAME] in self.catalog.source_syns[rep]:
                    kwargs[SOURCE.NAME] = rep
                    break

        if SOURCE.URL in kwargs:
            for rep in self.catalog.url_redirs:
                if kwargs[SOURCE.URL] in self.catalog.url_redirs[rep]:
                    kwargs[SOURCE.URL] = rep
                    break

        return super(FastStars, self).add_source(**kwargs)

    def priority_prefixes(self):
        """Prefixes to given priority to when merging duplicate entries.
        """
        return ('AT', 'SN')

    def add_self_source(self):
        return self.add_source(
            bibcode=self.catalog.OSC_BIBCODE,
            name=self.catalog.OSC_NAME,
            url=self.catalog.OSC_URL,
            secondary=True)

    def _get_save_path(self, bury=False):
        """Return the path that this Entry should be saved to.

        Determines output repository based on the name (i.e. the year) of the
        hypervelocity.  If `bury` is true, then this entry is saved to the
        'boneyard'.
        """
        self._log.debug("_get_save_path(): {}".format(self.name()))
        filename = self.get_filename(self[self._KEYS.NAME])

        # Put non-hypervelocity in the boneyard
        if bury:
            outdir = self.catalog.PATHS.get_repo_boneyard()

        # Get normal repository save directory
        else:
            repo_folders = self.catalog.PATHS.get_repo_output_folders(bones=False)
            outdir = repo_folders[0]

        return outdir, filename

    def sanitize(self):
        super(FastStars, self).sanitize()

        # Calculate some columns based on imported data, sanitize some fields
        name = self[self._KEYS.NAME]
        aliases = self.get_aliases()


        if self._KEYS.SOURCES in self:
            for source in self[self._KEYS.SOURCES]:
                if SOURCE.BIBCODE in source:
                    import urllib
                    from html import unescape
                    # First sanitize the bibcode
                    if len(source[SOURCE.BIBCODE]) != 19:
                        source[SOURCE.BIBCODE] = urllib.parse.unquote(
                            unescape(source[SOURCE.BIBCODE])).replace('A.A.',
                                                                      'A&A')
                    if source[SOURCE.BIBCODE] in self.catalog.biberror_dict:
                        source[SOURCE.BIBCODE] = \
                            self.catalog.biberror_dict[source[SOURCE.BIBCODE]]

                    if (source[SOURCE.BIBCODE] not in
                            self.catalog.bibauthor_dict):
                        bibcode = source[SOURCE.BIBCODE]
                        adsquery = (self.catalog.ADS_BIB_URL +
                                    urllib.parse.quote(bibcode) +
                                    '&data_type=Custom&format=%253m%20%25(y)')
                        bibcodeauthor = ''
                        try:
                            response = urllib.request.urlopen(adsquery)
                            html = response.read().decode('utf-8')
                            hsplit = html.split("\n")
                            if len(hsplit) > 5:
                                bibcodeauthor = hsplit[5]
                        except:
                            pass

                        if not bibcodeauthor:
                            warnings.warn(
                                "Bibcode didn't return authors, not converting"
                                "this bibcode.")

                        self.catalog.bibauthor_dict[bibcode] = unescape(
                            bibcodeauthor).strip()

            for source in self[self._KEYS.SOURCES]:
                if (SOURCE.BIBCODE in source and
                        source[SOURCE.BIBCODE] in self.catalog.bibauthor_dict
                        and
                        self.catalog.bibauthor_dict[source[SOURCE.BIBCODE]]):
                    source[SOURCE.REFERENCE] = self.catalog.bibauthor_dict[
                        source[SOURCE.BIBCODE]]
                if (SOURCE.NAME not in source and SOURCE.BIBCODE in source and
                        source[SOURCE.BIBCODE]):
                    source[SOURCE.NAME] = source[SOURCE.BIBCODE]

        #if self._KEYS.REDSHIFT in self:
        #    self[self._KEYS.REDSHIFT] = list(
        #        sorted(
        #            self[self._KEYS.REDSHIFT],
        #            key=lambda q: frame_priority(q, self._KEYS.REDSHIFT)))

        if self._KEYS.VELOCITY in self:
            self[self._KEYS.VELOCITY] = list(
                sorted(
                    self[self._KEYS.VELOCITY],
                    key=lambda q: frame_priority(q, self._KEYS.VELOCITY)))

        if self._KEYS.CLAIMED_TYPE in self:
            self[self._KEYS.CLAIMED_TYPE] = self.ct_list_prioritized()

        # Renumber and reorder sources
        if self._KEYS.SOURCES in self:
            # Sort sources reverse-chronologically
            self[self._KEYS.SOURCES] = sorted(
                self[self._KEYS.SOURCES], key=lambda x: bib_priority(x))

            # Assign new aliases to match new order
            source_reps = OrderedDict(
                [[x[SOURCE.ALIAS], str(i + 1)]
                 for i, x in enumerate(self[self._KEYS.SOURCES])])
            for i, source in enumerate(self[self._KEYS.SOURCES]):
                self[self._KEYS.SOURCES][i][SOURCE.ALIAS] = source_reps[source[
                    SOURCE.ALIAS]]

            # Change sources to match new aliases
            for key in self.keys():
                if self._KEYS.get_key_by_name(key).no_source:
                    continue
                for item in self[key]:
                    aliases = [
                        str(y)
                        for y in sorted(
                            int(source_reps[x])
                            for x in item[item._KEYS.SOURCE].split(','))
                    ]
                    item[item._KEYS.SOURCE] = ','.join(aliases)

    def clean_internal(self, data):
        """Clean input data from the 'FastStars/input/internal' repository.

        FIX: instead of making changes in place to `dirty_event`, should a new
             event be created, values filled, then returned??
        FIX: currently will fail if no bibcode and no url
        """
        self._log.debug("clean_internal(): {}".format(self.name()))

        def_source_dict = {}
        # Find source that will be used as default
        sources = data.get(self._KEYS.SOURCES, [])
        if sources:
            def_source_dict = sources[0]
            allow_alias = False
            if SOURCE.ALIAS in def_source_dict:
                del (def_source_dict[SOURCE.ALIAS])
        else:
            # If there are no existing sources, add OSC as one
            self.add_self_source()
            sources = self.get(self._KEYS.SOURCES, [])
            def_source_dict = sources[0]
            allow_alias = True


        dist_key = 'distinctfrom'
        if dist_key in data:
            distincts = data.pop(dist_key)
            if ((isinstance(distincts, list) and
                 isinstance(distincts[0], string_types))):
                source = self.add_self_source()
                for df in distincts:
                    self.add_quantity(self._KEYS.DISTINCT_FROM, df, source)
            else:
                data[dist_key] = list(distincts)

        # Go through all remaining keys in 'dirty' event, and make sure
        # everything is a quantity with a source (OSC if no other)
        for key in data.keys():
            # The following line should be used to replace the above once keys
            # returns the superclass keys too
            if self._KEYS.get_key_by_name(key).no_source:
                pass
            elif key == self._KEYS.PHOTOMETRY:
                for p, photo in enumerate(data[self._KEYS.PHOTOMETRY]):
                    if photo.get(PHOTOMETRY.U_TIME) == 'JD':
                        data[self._KEYS.PHOTOMETRY][p][
                            PHOTOMETRY.U_TIME] = 'MJD'
                        data[self._KEYS.PHOTOMETRY][p][PHOTOMETRY.TIME] = str(
                            jd_to_mjd(Decimal(photo['time'])))
                    if QUANTITY.SOURCE not in photo:
                        if not def_source_dict:
                            raise ValueError("No sources found, can't add "
                                             "photometry.")
                        source = self.add_source(
                            allow_alias=allow_alias, **def_source_dict)
                        data[self._KEYS.PHOTOMETRY][p][
                            QUANTITY.SOURCE] = source
            else:
                for qi, quantity in enumerate(data[key]):
                    if QUANTITY.SOURCE not in quantity:
                        if not def_source_dict:
                            raise ValueError("No sources found, can't add "
                                             "quantity.")
                        source = self.add_source(
                            allow_alias=allow_alias, **def_source_dict)
                        data[key][qi][QUANTITY.SOURCE] = source

        return data


    def set_preferred_name(self):
        """Set preferred name of faststar.

        Highest preference goes to names of the form 'SN####AA'.
        Otherwise base the name on whichever survey is the 'discoverer'.

        FIX: This function needs to be heavily modified
        """
        name = self[self._KEYS.NAME]
        newname = ''
        aliases = self.get_aliases()
        # if there are no other options to choose from, skip
        if len(aliases) <= 1:
            return name
        # If the name is already in the form 'SN####AA' then keep using
        # that
        #if (name.startswith('SN') and
        #    ((is_number(name[2:6]) and not is_number(name[6:])) or
        #     (is_number(name[2:5]) and not is_number(name[5:])))):
        #    return name
        # Otherwise, use the shortest name.
        if not newname:
            newname = min(aliases, key=len)
        if newname and name != newname:
            file_entry = None
            # Make sure new name doesn't already exist
            if newname in self.catalog.entries:
                if self.catalog.entries[newname]._stub:
                    file_entry = self.init_from_file(
                        self.catalog, name=newname)
                else:
                    file_entry = self.catalog.entries[newname]

            if file_entry:
                self._log.info("`{}` already exists, copying `{}` to it".
                               format(newname, name))
                # Douglas had to add this try-except because some entries had already been deleted.
                try:
                    self.catalog.copy_entry_to_entry(
                        self.catalog.entries[name], file_entry)
                    del self.catalog.entries[name]
                except KeyError:
                    self._log.info("`{}` has already been coped to `{}`".
                               format(name, newname))
                self.catalog.entries[newname] = file_entry
            else:
                self._log.info("Changing entry from name '{}' to preferred"
                               " name '{}'".format(name, newname))
                self.catalog.entries[newname] = self.catalog.entries[name]
                self.catalog.entries[newname][self._KEYS.NAME] = newname
                del self.catalog.entries[name]
            return newname

        return name

    def ct_list_prioritized(self):
        ct_list = list(
            sorted(
                self[self._KEYS.CLAIMED_TYPE],
                key=lambda key: self._ct_priority(key)))
        return ct_list

    def _ct_priority(self, attr):
        aliases = attr['source'].split(',')
        max_source_year = -10000
        vaguetypes = ['CC', 'I']
        if attr[QUANTITY.VALUE] in vaguetypes:
            return -max_source_year
        for alias in aliases:
            if alias == 'D':
                continue
            source = self.get_source_by_alias(alias)
            if SOURCE.BIBCODE in source:
                source_year = get_source_year(source)
                if source_year > max_source_year:
                    max_source_year = source_year
        return -max_source_year
