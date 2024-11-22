# general class for exoplanet populations
from ..imports import *

# from ..telescopes import *
# from ..models import *
from .column_descriptions import *
from .pineda_skew import make_skew_samples_from_lowerupper

# these are keywords that can be set for a population
default_plotkw = dict(
    color="black",
    alpha=1,
    zorder=0,
    marker="o",
    linewidth=1,
    respond_to_color=True,
    respond_to_size=True,
    exact=False,
    label_planets=False,
    filled=True,
    outlined=False,
)

# what keywords can we set for the population plotkw?
allowed_plotkw = list(default_plotkw.keys())
allowed_plotkw += [
    "s",
    "c",
    "cmap",
    "norm",
    "vmin",
    "vmax" "outlined",
    "filled",
    "markeredgewidth",
    "markeredgecolor",
]


def _clean_column(raw_column):
    """
    Convert a Column into an array disconnected from its original Table.

    Parameters
    ----------
    raw_column : astropy.table.column.Column
        The column to clean.

    Returns
    -------
    cleaned_column : array, Quantity, Time, SkyCoord, ...
        The data as an array-like object, without
        connection to its original table.
    """
    # strip the connection to the table column
    if isinstance(raw_column, Column):
        cleaned_column = np.array(raw_column)
    else:
        cleaned_column = type(raw_column)(raw_column)
        # This is equivalent to
        #
        # if isinstance(raw_column, u.Quantity):
        #   cleaned_column = u.Quantity(raw_column)
        #
        # I'm not 100% sure why this is necessary, but
        # in some tests it seems a QTable will return a
        # Quantity or SkyCoord or Time object, but still
        # with some sort of connection to being a Column
        # from a Table (and therefore potentially some
        # sneaky link that we want to erase). Re-initializing
        # an object seems to erase that connection.
    return cleaned_column


class Population(Talker):
    """
    Populations of astronomical objects might contain
        Exoplanets (planets with host stars),
        Solar System objects (major/minor planets),
        Stars (single or multiples),
    or more!
    """

    def __init__(self, standard, label=None, **plotkw):
        """
        Initialize a Population of exoplanets from a standardized table.

        This creates a population from a standardized data table,
        effectively just storing a table in the right format inside
        the object. To simplify indexing via system name without having
        to stress about capitalization and/or spaces, the population
        will define a `tidyname` and `tidyhostname` columns, and use
        those as table indices for searching via name.

        Parameters
        ----------
        standard : astropy.table.QTable, str
            If a Table, the standardized data table.
            If a string, the filename to the standardized data table.
        label : str
            A string by which this table can be identified. This will be
            used for plots and for saving data files, so please try to
            pick something informative and unique!
        **plotkw : dict
            All other keyword arguments will be taken as plotting suggestions.
        """

        if isinstance(standard, Table) or isinstance(standard, Row):
            # a standardized table with a minimum set of columns we can expect
            self.standard = QTable(copy.deepcopy(standard))
            # (the deepcopy seems to be needed to reset indexing,
            #  so that astropy.table doesn't try to track slices across
            #  too many subsets, because it was getting lost in weird ways)

            # store a label for this population
            self.label = label

            # keywords to use for plotting
            self._plotkw = plotkw

        elif isinstance(standard, str):
            filename = standard
            self.standard = QTable(ascii.read(filename))
            self.label = self.standard.meta["label"]
            self._plotkw = self.standard.meta["plotkw"]

        # define some cleaned names and hostnames, for indexing
        try:
            self.standard["tidyname"]
        except KeyError:
            self.standard["tidyname"] = [
                clean(x).lower() for x in self.standard["name"]
            ]
        try:
            self.standard["tidyhostname"]
        except KeyError:
            self.standard["tidyhostname"] = [
                clean(x).lower() for x in self.standard["hostname"]
            ]

        # make sure the table is searchable via names
        self._make_sure_index_exists("tidyname")
        self._make_sure_index_exists("tidyhostname")

        # test that indexing still works
        name = self.tidyname[0]
        self.standard.loc[name]

        # define internal lists of column names
        self._populate_colname_summary()

    def _populate_colname_summary(self):
        """
        Populate a dictionary summarizing the types of columns in `.standard`.
        """
        uncertainty_suffixes = [
            "_uncertainty_lower",
            "_uncertainty_upper",
            "_uncertainty",
        ]
        limit_suffixes = ["_lower_limit", "_upper_limit"]
        reference_suffixes = ["_reference"]
        all_suffixes = uncertainty_suffixes + limit_suffixes + reference_suffixes

        def ends_with(s, suffixes):
            return np.any([s[-len(suffix) :] == suffix for suffix in suffixes])

        def remove_suffixes(s):
            for suffix in all_suffixes:
                s = s.replace(suffix, "")
            return s

        self._colnames = {}
        self._colnames["everything"] = np.unique(
            [remove_suffixes(x) for x in self.standard.colnames]
        )
        self._colnames["with uncertainties"] = np.unique(
            [
                remove_suffixes(x)
                for x in self.standard.colnames
                if ends_with(x, uncertainty_suffixes)
            ]
        )
        self._colnames["without uncertainties"] = np.unique(
            [
                x
                for x in self._colnames["everything"]
                if x not in self._colnames["with uncertainties"]
            ]
        )
        # self._colnames["with limits"] = np.unique(
        #    [
        #        remove_suffixes(x)
        #        for x in self.standard.colnames
        #        if ends_with(x, limit_suffixes)
        #    ]
        # )
        self._colnames["with references"] = np.unique(
            [
                remove_suffixes(x)
                for x in self.standard.colnames
                if ends_with(x, reference_suffixes)
            ]
        )

    def print_column_summary(self):
        """
        Print a summary of columns that come directly from the `.standard` table.
        """
        print(f"The following columns are present in the internal `.standard` table:\n")
        for k in self._colnames:
            print(f"{k} =\n{self._colnames[k]}\n")

    def _list_table_indices(self):
        """
        Return a list of keys being used as table indices.

        Core populations should likely have just ['tidyname', 'tidyhostname']
        but internal `.individual_references` populations could have lots,
        one for each quantity with an reference that might be chosen.

        Returns
        -------
        index_keys : list
            The list of keys that are being used as an index.
        """
        return [x.columns[0].name for x in self.standard.indices]

    def _make_sure_index_exists(self, k):
        """
        Add a new key as an index, but don't add it twice.

        Parameters
        ----------
        k : str
            The new key to add.
        """
        if k not in self._list_table_indices():
            self.standard.add_index(k)

    @property
    def _fileprefix(self):
        """
        Define a fileprefix for this population, to be used
        for setting the filename of the standardized population.

        Return
        ------
        """
        return clean(self.label)

    @property
    def _standardized_data_path(self):
        """
        Define the filepath for the standardized table.
        """
        return os.path.join(directories["data"], f"standardized-{self._fileprefix}.txt")

    def save(self, filename=None, overwrite=True):
        """
        Save this population out to a file.

        This saves the standardized data table from this population
        out to a file, along with metadata needed to be loaded
        back into as a drop-in replacement for a live population.

        Parameters
        ----------
        filename : str
            The filepath to which the population should be saved.
        overwrite : bool
            Whether to automatically overwrite existing populations.

        Examples
        --------
        >>> from exoatlas import *
        >>> one = Exoplanets()['GJ1132b']
        >>> one.save('pop.txt')
        >>> the_same_one = Population('pop.txt')
        """

        # be a little fussy about overwriting automatic filenames
        if filename == None:
            filename = f"exoatlas-population-{self._fileprefix}.ecsv"
            overwrite = False
            if os.path.exists(filename):
                warnings.warn(
                    f"{filename} will not be overwritten unless you explicitly provide a filename."
                )

                # save the table as an ascii table for humans to read
        to_save = copy.deepcopy(self.standard)
        to_save.meta["label"] = self.label
        to_save.meta["plotkw"] = self._plotkw

        to_save.write(filename, format="ascii.ecsv", overwrite=overwrite)
        print(
            f"""
        Saved {self} to {filename}.
        It can be reloaded with `x = Population('{filename}')`
        """
        )

    def sort(self, x, reverse=False):
        """
        Sort this population by some key or attribute.

        This sorts the population in place, meaning that the
        Population object from which it is called will be modified.
        Nothing will be returned.

        Parameters
        ----------
        x : str
            The key by which to sort the table.
        reverse : bool
            Whether to reverse the sort order.
                `reverse == False` means low to high
                `reverse == True` means high to low
        """

        # get the values by which to sort the population
        to_sort = getattr(self, x)

        # define the sorting indices
        i = np.argsort(to_sort)
        if reverse:
            i = i[::-1]

        # reorder the standardized table
        self.standard = self.standard[i]

    def __add__(self, other):
        """
        Create a new population by adding two together:

            `bigger = this + other`

        Parameters
        ----------
        other : Population
            The population to be tacked onto this one.

        Returns
        -------
        bigger : Population
            A new population, consisting of all the planets
            in `this` population and some extra ones added
            from `other`.

        """

        # skip any warnings that pop up
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

            #  create a new table, joining both together
            table = join(
                self.standard.filled(), other.standard.filled(), join_type="outer"
            )
            # TO-DO, I'm not 100% sure why the tables need to be `filled()` here; shouldn't they already be?

            # create an informative label
            label = f"{self.label} + {other.label}"

        # create and return the new population
        return type(self)(standard=table, label=label)

    def _remove_by_key(self, other, key="tidyname"):
        """
        Create a new population by removing some rows from here:

            `smaller = this - other`

        Parameters
        ----------
        other : Population
            The population of planets to be removed from
            `this` population to create a new `smaller` one.

        Returns
        -------
        smaller : Population
            A subset of `this` population, where some rows
            have been removed.
        """

        # skip any warnings that pop up
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

            #  create a new table, joining both together
            table = setdiff(self.standard, other.standard, keys=key)

            # create an informative label
            label = f"{self.label} - {other.label}"

        # create and return the new population
        return type(self)(standard=table, label=label)

    def __sub__(self, other):
        """
        Create a new population by removing some rows from here:

            `smaller = this - other`

        Parameters
        ----------
        other : Population
            The population of planets to be removed from
            `this` population to create a new `smaller` one.

        Returns
        -------
        smaller : Population
            A subset of `this` population, where some rows
            have been removed.
        """
        return self._remove_by_key(other)

    def __getitem__(self, key):
        """
        Create a subpopulation of planets by indexing, slicing, or masking.

        Parameters
        ----------
        key : int, list, array, slice
            An index, slice, or mask to select a subset of the population.

        Returns
        -------
        subset : Population
            A subset of the population, as set by the index, slice, or mask.
        """

        try:
            # if the key is an index/slice/mask, return it
            if self.label is None:
                label = None
            else:
                label = f"Subset of {self.label}"
            subset = type(self)(
                standard=self.standard[key], label=label, **self._plotkw
            )

            # if the key is a column, raise an error
            if type(key) in self.standard.colnames:
                raise IndexError(
                    f"""
                You seem to be trying to access a column from this
                population via `pop[{key}]`. For clarity, all `[]`
                indexing is reserved for selecting subsets of the
                population.

                To access your particular column, please try either
                `pop.{key}` or `pop.standard[{key}]` to return a
                1D array of the entries in that column.
                """
                )
        except KeyError:
            # use a string or a list of strings make a subset by planet name
            # FIXME - maybe we should make this say more when it's making a sneaky choice for us?
            try:
                subset = self.create_subset_by_name(key)
            except KeyError:
                subset = self.create_subset_by_hostname(key)

        return subset

    def create_subset_by_name(self, key):
        """
        Extract a subset of this population,
        based on one or more planet names.

        Parameters
        ----------
        key : strings, list of strings
            The name of a planet ("GJ1132b")
            or a list of planet names.

            (All names will be stripped of
            special characters and converted
            to lower case before indexing.)

        Returns
        -------
        subset : Population
            A new population containing
            some subset of the original.
        """

        # use a (list of) string(s) to index population by name
        if isinstance(key, str):
            # is it just one name?
            tidy = clean(key).lower()
        elif isinstance(key[0], str):
            # is it a list of names?
            tidy = [clean(k).lower() for k in key]

        # pull out rows by planet name
        subset = self.standard.loc["tidyname", tidy]

        # create a useful label for the population
        if isinstance(key, str):
            label = key
        elif isinstance(key[0], str):
            label = "+".join(key)

        # create that new sub-population
        return type(self)(standard=subset, label=label, **self._plotkw)

    def create_subset_by_hostname(self, key):
        """
        Extract a subset of this population,
        based on one or more planet hostnames.

        Parameters
        ----------
        key : strings, list of strings
            The hostname of a planet ("GJ1132")
            or a list of planet hostnames.

            (All names will be stripped of
            special characters and converted
            to lower case before indexing.)

        Returns
        -------
        subset : Population
            A new population containing
            some subset of the original.
        """

        # use a string or a list of strings to index the population by name
        if isinstance(key, str):
            # is it just one name?
            tidy = clean(key).lower()
        elif isinstance(key[0], str):
            # is it a list of names?
            tidy = [clean(k).lower() for k in key]

        # pull out rows by planet name
        subset = self.standard.loc["tidyhostname", tidy]

        # create a useful label for the population
        if isinstance(key, str):
            label = key
        elif isinstance(key[0], str):
            label = "+".join(key)

        # create that new sub-population
        return type(self)(standard=subset, label=label, **self._plotkw)

    def create_subset_by_position(
        self,
        coordinates,
        radius=1 * u.arcmin,
        use_proper_motions=False,
        return_indices=False,
    ):
        """
        Extract a subset of this population,
        by performing a spatial cross-match by
        RA and Dec. This will return all planets
        from this population that fall within
        the specified radius of at least one of
        the specified coordinates.

        Parameters
        ----------
        coordinates : astropy.coordinates.SkyCoord
            The sky coordinate (RA, Dec) or list
            of coordinates we want to search for
            nearby objects.

        radius : astropy.units.Quantity
            The angular radius around each position
            that we should include in each search.

        use_proper_motions : bool
            Should we use available proper motions,
            embedded in the skycoords, to propagate
            positions to a shared epoch before
            cross-matching? Alas, this ability
            is *not yet implemented*. FIXME!

        return_indices : bool
            Should we also return the indices
            of the original coordinates that
            were matched to existing positions?

        Returns
        -------
        subset : Population
            A new population containing a subset
            of the original, including *all* planets
            that fall within the 2D sky search space.

        """

        if use_proper_motions:
            raise NotImplementedError("No cross-matching with proper motions yet :-(")

        # create astropy coordinates for this population
        population_coords = SkyCoord(ra=self.ra, dec=self.dec)

        # do a spatial cross match on the sky
        #  (idx gives the index into coordinates,
        #   each corresponding to an entry in population_coords)
        idx, d2d, d3d = population_coords.match_to_catalog_sky(coordinates)

        # identify which systems are actually close on the sky
        match = d2d < radius

        # create new populations that are linked by spatial position
        i_match = match.nonzero()[0]
        # matched_coordinates = coordinates[idx[i_match]]
        subset = self.standard[i_match]

        # define a meaningful label
        label = f"Spatial Cross-Match ({len(coordinates)} positions, {radius} radius)"

        # create that new sub-population
        new_population = type(self)(standard=subset, label=label, **self._plotkw)

        # choose what to return
        if return_indices:
            i_from_original_coordinates = idx[i_match]
            return new_population, i_from_original_coordinates
        else:
            return new_population

    def get_lowerupper_uncertainty_from_table(self, key):
        """
        Return two arrays of lower and upper uncertainties,
        direct from a table.

        Parameters
        ----------
        key : str
            The quantity for which we want lower + upper uncertaintes.

        Returns
        -------
        lower : np.array, u.Quantity
            The magnitude of the lower uncertainties (x_{-lower}^{+upper})
        upper : np.array, u.Quantity
            The magnitude of the upper uncertainties (x_{-lower}^{+upper})
        """

        # first try for asymmetric table uncertainties
        try:
            lower = _clean_column(self.standard[f"{key}_uncertainty_lower"])
            upper = _clean_column(self.standard[f"{key}_uncertainty_upper"])
            return np.abs(lower), np.abs(upper)
        except KeyError:
            pass

        # next try for symmetric table uncertainties
        try:
            sym = _clean_column(self.standard[f"{key}_uncertainty"])
            return np.abs(sym), np.abs(sym)
        except KeyError:
            pass

        # then give up and return zeroes
        unc = 0 * _clean_column(self.standard[key])
        return unc, unc

    def get_values_from_table(self, key, distribution=False):
        """
        Retrieve values directly from the standardized table.

        This wrapper extracts values from the internal `.standard`
        table without performing any calculations or filtering.
        Some quantities might have explicit methods that override
        direct retrieval from the table, but this provides direct
        access to the table values no matter what.

        Parameters
        ----------
        key : str
            The quantity to extract. This must exactly match a
            column in `.standard`; if not, `KeyError` will be raised.
        distribution : bool
            Should it be returned as an astropy distribution,
            for uncertainty propagation?

        Returns
        -------
        values : array, astropy.units.Quantity, astropy.uncertainty.Distribution
            The requested table column.
                If it has units, it will be an astropy Quantity.
                If it has no units, it will be an astropy array.
                If `distribution` is True, it will be an astropy Distribution.
        """

        # extract the column from self.standard for this key
        try:
            raw_column = self.standard[key]
        except KeyError:
            raise KeyError(
                f"The column '{key}' wasn't found in {self}'s internal table."
            )

        # return a distribution or an array
        if distribution:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                # set the mode of the distribution to be the table value
                mu = _clean_column(raw_column)
                # get the lower + upper uncertainties from the table
                sigma_lower, sigma_upper = self.get_lowerupper_uncertainty_from_table(
                    key
                )
                samples = make_skew_samples_from_lowerupper(
                    mu=mu, sigma_lower=sigma_lower, sigma_upper=sigma_upper
                )
                return Distribution(samples)
        else:
            return _clean_column(raw_column)

    def __getattr__(self, key):
        """
        If an attribute/method isn't explicitly defined for a population,
        look for it as a column of the standardized table or as a plotting keyword.

        For example, `population.stellar_radius` will try to
        access `population.standard['stellar_radius']`.

        Parameters
        ----------
        key : str
            The attribute we're trying to get.

        Returns
        -------
        value : Quantity, str
            The array of quantities from the standardized table,
            or the plotting keyword.

        """
        if key in ["label", "_plotkw", "standard"]:
            raise RuntimeError(
                f"Yikes! It looks like `.{key}` isn't defined, but it should be!"
            )
        try:
            # extract the column from the standardized table
            return self.standard[key]
        except KeyError:
            # try to get a plotkw from this pop, from the plotting defaults, from None
            try:
                assert key in allowed_plotkw
                return self._plotkw.get(key, default_plotkw[key])
            except (AssertionError, KeyError):
                raise AttributeError(
                    f"""
                Alas, there seems to be no way to find `.{key}`
                as a table column, attribute, method, or property of {self}.
                """
                )

    def __setattr__(self, key, value):
        """
        Define what happens when we try to set an attribute via `pop.attr = x`.

        If the keyword is a pre-defined "plotting" keyword in `allowed_plotkw`,
        then we should save it in a special `_plotkw` dictionary. Otherwise,
        the attribute should be set as normal.

        Parameters
        ----------
        key : str
            The attribute we're trying to set.
        value : anything
            The value we're trying to give that attribute.
        """

        if key in allowed_plotkw:
            # store plotting keywords in a separate plotting dictionary
            self._plotkw[key] = value
        else:
            # otherwise, store attributes as normal for objects
            self.__dict__[key] = value

    def __repr__(self):
        """
        How should this object appear as a repr/str?
        """
        return f"✨ {self.label} | {len(self)} elements ✨"

    def get(self, key):
        """
        Return an array of values for a column.

        Parameters
        ----------
        key : str
            The name of the quantity we're trying to retrieve.

        Returns
        -------
        value : Quantity
            The values, as an array with units.
        """
        return getattr(self, key)

    def get_uncertainty(self, key):
        """
        Return an array of symmetric uncertainties on a column.

        Parameters
        ----------
        key : str
            The quantity for which we want uncertaintes.
        uncertainty : Quantity
            The uncertainties, as an array with units.
        """

        # first try for an `uncertainty_{key}` column
        try:
            return self.get(f"{key}_uncertainty")  # __getattr__
        except (
            KeyError,
            AssertionError,
            AtlasError,
            AttributeError,
        ):  # is including AttributeError a kludge?
            # this can be removed after debugging
            self._speak(f'no symmetric uncertainties found for "{key}"')

        # then try for crudely averaging asymmetric uncertainties
        try:
            lower = self.get(f"{key}_uncertainty_lower")  # __getattr__
            upper = self.get(f"{key}_uncertainty_upper")  # __getattr__
            avg = 0.5 * (np.abs(lower) + np.abs(upper))
            return avg
        except (KeyError, AssertionError, AtlasError, AttributeError):
            # this can be removed after debugging
            self._speak(f'no asymmetric uncertainties found for "{key}"')

        # then give up and return nans
        return np.nan * self.standard[key]

    def get_uncertainty_lowerupper(self, key):
        """
        Return two arrays of lower and upper uncertainties on a column.

        Parameters
        ----------
        key : str
            The quantity for which we want lower + upper uncertaintes.

        Returns
        -------
        lower : np.array
            The magnitude of the lower uncertainties (x_{-lower}^{+upper})
        upper : np.array
            The magnitude of the upper uncertainties (x_{-lower}^{+upper})
        """

        # first try for actual asymmetric uncertainties
        try:
            lower = self.__getattr__(f"{key}_uncertainty_lower")
            upper = self.__getattr__(f"{key}_uncertainty_upper")
            return np.abs(lower), np.abs(upper)
        except (KeyError, AssertionError, AttributeError):
            # this can be removed after debugging
            self._speak(f'no asymmetric uncertainties found for "{key}"')

        # first try for an `uncertainty_{key}` column
        try:
            sym = self.__getattr__(f"{key}_uncertainty")
            return np.abs(sym), np.abs(sym)
        except (KeyError, AssertionError, AttributeError):
            # this can be removed after debugging
            self._speak(f'no symmetric uncertainties found for "{key}"')

        # then give up and return nans
        unc = np.nan * self.__getattr__(key)
        return unc, unc

    def _validate_columns(self):
        """
        Make sure this standardized table has all the necessary columns.
        Summarize the amount of good data in each.
        """

        N = len(self.standard)
        for k in basic_columns:
            try:
                n = sum(self.standard[k].mask == False)
            except AttributeError:
                try:
                    n = sum(np.isfinite(self.standard[k]))
                except TypeError:
                    n = sum(self.standard[k] != "")
            self._speak(f"{k:>25} | {n:4}/{N} rows = {n/N:4.0%} are not empty")

    def _find_index(self, name):
        """
        Return index of a particular planet in the population.

        ??? = maybe this could/should be replaced with some table cleverness?
        """

        return np.array([clean(name) in clean(x) for x in self.name]).nonzero()[0]

    def update_values(self, planets, **kwargs):
        """
        Update values for one or more planets.

        This modifies the internal `.standard` table
        to update individual values. This is meant to
        be a tool that can be used to provide alternate,
        better, and/or unpublished planet parameters
        to a population.

        Changes made by this function will only last for
        the duration of the Python session in which they
        are run. If you want them to be more permanent,
        you will need to save out the population to a
        custom curated standardized file.

        Parameters
        ----------
        planets : str, list
            Names of one or more planets to update.
        kwargs : dict
            All other keyword arguments will be passed to update
            the values for the planet(s). Quantities should have
            appropriate units.

            If only one planet eis being updated, quantities should
            each be single values. If N planets are being updated,
            quantities should be (N,)-dimensional arrays in the
            same order as the planets.
        """

        # create `tidyname` keys to index by planet name
        if isinstance(planets, str):
            planets_to_index = [clean(planets).lower()]
        else:
            planets_to_index = [clean(p).lower() for p in planets]

        # extract just the subsection of the table relating to these planets
        i = self.standard.loc_indices[planets_to_index]

        # loop over keyword arguments
        for k, v in kwargs.items():

            # make sure it's a valid keyword value to assign
            # assert k in subset.colnames

            if k[-12:] == "_uncertainty":
                # nudge symmetric uncertainties into asymmetric form
                for suffix, sign in zip(["_lower", "_upper"], [-1, 1]):
                    old = self.standard[i][k + suffix] * 1
                    new = sign * np.abs(v)
                    self.standard[i][k + suffix] = new
                    print(f"{planets_to_index} | {k+suffix}: {old} > {new}")
            else:
                # update value in table
                old = self.standard[i][k] * 1
                new = v
                self.standard[i][k] = new
                print(f"{planets_to_index} | {k}: {old} > {new}")

            # warn if uncertainties should have been provided but weren't
            should_have_uncertainty = f"{k}_uncertainty_lower" in self.standard.colnames
            does_have_uncertainty = (f"{k}_uncertainty" in kwargs) or (
                (f"{k}_uncertainty_lower" in kwargs)
                and (f"{k}_uncertainty_upper" in kwargs)
            )
            if should_have_uncertainty and (does_have_uncertainty == False):
                warnings.warn(
                    f"'{k}' should probably have some uncertainties, which you didn't provide."
                )

    def __len__(self):
        """
        How many planets are in this population?
        """
        return len(self.standard)

    def create_table(
        self,
        desired_columns=[
            "name",
            "radius",
            "relative_insolation",
            "stellar_radius",
            "stellar_teff",
            "ra",
            "dec",
            "distance",
        ],
    ):
        """
        Create an astropy table based on this population,
        using a subset of columns, which may include ones
        that have been calculated as Population properties.

        Parameters
        ----------
        desired_columns : list
            The columns you want to include. Anything that
            can be accessed via Population.??? can be provided
            here as a string.

        Returns
        -------
        table : astropy.table.QTable
            A table, with those columns, in the same order
            as the Population itself.
        """
        # FIXME! need to add method support for arguments

        # create a dictionary with the desired columns
        d = {c: getattr(self, c) for c in desired_columns}

        # turn that into an astropy Table
        t = QTable(d)

        return t

    def create_planning_table(
        self,
        desired_columns=[
            "name",
            "ra",
            "dec",
            "period",
            "transit_midpoint",
            "transit_duration",
            "radius",
            "relative_insolation",
            "stellar_radius",
            "stellar_teff",
            "distance",
        ],
    ):
        """
        Create an astropy table to help plan transit observations.

        Parameters
        ----------
        desired_columns : list
            The columns you want to include. Anything that
            can be accessed via Population.??? can be provided
            here as a string.

        Returns
        -------
        table : astropy.table.QTable
            A table, with those columns, in the same order
            as the Population itself.
        """
        return self.create_table(desired_columns=desired_columns)
