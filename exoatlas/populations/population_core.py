# general class for exoplanet populations
from ..imports import *

# from ..telescopes import *
# from ..models import *
from .column_descriptions import *
from .pineda_skew import make_skew_samples_from_lowerupper, gaussian_central_1sigma

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
        name = self.standard["tidyname"][0]
        self.standard.loc[name]

        # define internal lists of column names
        self._populate_column_summaries()
        self._populate_column_methods()

        # define how many samples and iterations to use for uncertainty propagation
        self.targeted_fractional_uncertainty_precision = 0.05
        self._number_of_uncertainty_samples = 100 




    def _create_function_to_access_table_quantity(self, name):
        """
        For a given column name "x", add a `.x()` method to this Population.

        This internal wrapper creates a named method to access
        the data in a particular column of the standardized table.
        It's a little frivolous, but this mostly makes it easier to
        see variables via tab completion and for doing error propagation.
        This will be called at the initial creation of the Population,
        and whenever new columns are added via `.add_column`.

        Parameters
        ----------
        name : str
            The name of the column to add.
        """

        assert name in self.standard.colnames

        # create the method to extract a particular column
        def f(distribution=False, **kw):
            return self.get_values_from_table(name, distribution=distribution)

        # build a basic docstring for that method.
        try:
            unit = self.get_values_from_table(name).unit.to_string()
            unit_string = f", with units of '{unit}'"
        except AttributeError:
            unit_string = ""

        f.__doc__ = f"""
            The table quantity '{name}'{unit_string}.

            Parameters 
            ----------
            distribution : bool 
                Should we try to return an astropy Distribution, 
                for uncertainty propagation? If the table contains 
                quoted uncertainties, a distribution will be returned; 
                if not, a simple array will be.
            kw : dict 
                All other keyword arguments will be ignored.

            Returns 
            -------
            values : array-like, or Distribution
                The values for this table column. 
                If `distribution==True`, these values will be as an 
                astropy Distribution, built from the table uncertainties.
            """
        return f

    def _populate_column_methods(self):
        """
        Populate this object with one method for each table column.

        This wrapper converts table columns from the internal .standard
        table into callable methods, so that (for example), the data
        in `self.standard['radius']` can be retrieved as `self.radius()`.
        This is necessary for seamlesssly integrating core table quantities
        with derived quantities (some of which must be callable because
        they require keyword arguments), for making them appear as hints
        with tab-completion in ipython/jupyter, and for being able to
        attach a docstring to each.
        """

        # add all the basic table quantities as methods to the population
        for k in self._colnames["everything"]:
            # if no other function is defined, populate a table one
            if hasattr(self, k) == False:
                setattr(self, k, self._create_function_to_access_table_quantity(k))

    def _populate_column_summaries(self):
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

        def ends_with_any_of(s, suffixes):
            return np.any([s.endswith(suffix) for suffix in suffixes])

        def remove_suffixes(s):
            for suffix in all_suffixes:
                s = s.split(suffix)[0]
            return s

        self._colnames = {}
        self._colnames["everything"] = np.unique(
            [remove_suffixes(x) for x in self.standard.colnames]
        )
        self._colnames["with uncertainties"] = np.unique(
            [
                remove_suffixes(x)
                for x in self.standard.colnames
                if ends_with_any_of(x, uncertainty_suffixes)
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
        #        if ends_with_any_of(x, limit_suffixes)
        #    ]
        # )
        self._colnames["with references"] = np.unique(
            [
                remove_suffixes(x)
                for x in self.standard.colnames
                if ends_with_any_of(x, reference_suffixes)
            ]
        )

    def add_column(
        self,
        name="",
        data=None,
        uncertainty=None,
        uncertainty_lower=None,
        uncertainty_upper=None,
    ):
        """
        Add a column to the existing population.

        This wrapper adds a new column of data to the `.standard` table,
        populates columns for uncertainties if provided, and registers
        a new method for accessing that data.

        Parameters
        ----------
        name : str
            The name of this column; it must be able to be a valid Python variable name.
            If "x", we will add the column `.standard["x"]` and `.x()` as methods.
        data : astropy.units.Quantity
            An array of data to be added into this column. It must have the same
            size at the population, and its entries should be ordered the same
            as the population.
        uncertainty : astropy.units.Quantity
            (optional) Symmetric uncertainties associated with the column.
            If included
        uncertainty_lower : astropy.units.Quantity
            (optional, along with uncertainty_upper) The magnitude of the
            lower uncertainty, for asymmetric uncertainties.
        uncertainty_upper : astropy.units.Quantity
            (optional, along with uncertainty_lower) The magnitude of the
            upper uncertainty, for asymmetric uncertainties.
        """

        # warn if overwriting an existing column
        if name in self.standard.colnames:
            warnings.warn(
                f"'{name}' already exists in `.standard`; you're overwriting something!"
            )

        # put the data into the column
        self.standard[name] = data

        # add uncertainties, if provided
        if uncertainty is not None:
            self.standard[f"{name}_uncertainty"] = uncertainty
        elif (uncertainty_lower is not None) and (uncertainty_upper is not None):
            self.standard["{name}_uncertainty_lower"] = uncertainty_lower
            self.standard["{name}_uncertainty_upper"] = uncertainty_upper

        # register the method for the column
        if hasattr(self, name):
            warnings.warn(
                f"'.{name}()' already exists in `.standard` for this Population; please consider a different name!"
            )
        setattr(self, name, self._create_function_to_access_table_quantity(name))

    def add_calculation(self, name, function):
        """
        Add a new calculation to this population.

        This wrapper adds a new method to calculate a new quantity,
        defining a new function that can be called with access to the
        population's internal data. Uncertainty propagation can also
        be applied to any new calculations registered using
        `.add_calculation`.

        Parameters
        ----------
        name : str
            The name of this calculation; it must be able to be a valid Python variable name.
            If "x", we will add `.x()` as as a new method.
        function : function
            Python function that calculates a quantity. This function
            should look something like:

                def f(self, distribution=False, **kw):
                    '''
                    Brief Name for Quantity (unit)

                    A more detailed description of the quantity,
                    making the docstring for the function more
                    readable for everyone who might use it.
                    '''
                    x = self.yay(distribution=distribution)
                    y = self.wow(distribution=distribution)
                    return x*y

            where the `distribution` keywords enable uncertainty propagation.
        """
        if hasattr(self, name):
            warnings.warn(
                f"Eep! You're overwriting `.{name}()`. That might be OK, but you should be aware!"
            )

        setattr(self.__class__, name, function)

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
        population_coords = SkyCoord(ra=self.ra(), dec=self.dec())

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

    def get_uncertainty_lowerupper_from_table(self, key):
        """
        Return two arrays of lower and upper uncertainties, directly from table columns.

        If no uncertainties of any kind are found, a KeyError should be raised.

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

        # this is a bit of a kludge, but when doing error propagation it 
        # might be possible for this to get called with either 'stellar_luminosity' 
        # or 'stellar_luminosity_from_table'. To avoid triggering a KeyError, 
        # we should make sure to strip '_from_table' from the key:
        key = key.replace('_from_table', '')

        # first try for asymmetric table uncertainties
        try:
            lower = _clean_column(self.standard[f"{key}_uncertainty_lower"])
            upper = _clean_column(self.standard[f"{key}_uncertainty_upper"])
            return np.abs(lower), np.abs(upper)
        except KeyError:
            # next try for symmetric table uncertainties
            sym = _clean_column(self.standard[f"{key}_uncertainty"])
            return np.abs(sym), np.abs(sym)

    def get_uncertainty_from_table(self, key, **kw):
        """
        Return an array of symmetric uncertainties on a quantity,
        directly from table columns.

        This (very crudely) averages the lower and upper uncertainties
        to make a symmetric errorbar. For quantities with nearly symmetric
        errors this should be totally fine, but for ones with wildly
        asymmetric uncertainties, use `.get_uncertainty_lowerupper_from_table()`.
        This returns an array, not a function.

        Parameters
        ----------
        key : str
            The quantity for which we want uncertaintes.

        Returns
        -------
        uncertainty : Quantity
            The (symmetrically-averaged) uncertainties, as an array with units.
        """

        sigma_lower, sigma_upper = self.get_uncertainty_lowerupper_from_table(key, **kw)
        return 0.5 * (sigma_lower + sigma_upper)

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
                try:
                    # get the lower + upper uncertainties from the table
                    sigma_lower, sigma_upper = (
                        self.get_uncertainty_lowerupper_from_table(key)
                    )
                    samples = make_skew_samples_from_lowerupper(
                        mu=mu, sigma_lower=sigma_lower, sigma_upper=sigma_upper, N_samples=self._number_of_uncertainty_samples
                    )
                    return Distribution(samples)
                except KeyError:
                    # if uncertainties fail, give up and return simple array
                    return mu
        else:
            return _clean_column(raw_column)

    def __getattr__(self, key):
        """
        If an attribute/method isn't explicitly defined for a population,
        try to access it as `population.something` with this wrapper.

        Lots of data is stored inside an exoatlas Population; this hidden
        method will be called to try to access those data only if no other
        explicit definition defines that variable first. If we try to access
        an internal Population variable via `pop.x` or `getattr(pop, 'x')`,
        here's what happens, shown approximately in order of precedence:

            1) We look for an explicit definition that has been
            attached to the Population as a method/attribute. This may have been
            defined either with a `def x()` definition inside a module in
            `exoatlas/populations/calculations/` or an attribute/method that
            gets defined anywhere with `pop.x =`, potentially including
            by a user creating a new quantity function and assigning it.
            These definitions should show up with `pop.<tab>` in jupyter.

            2) We look for an implicit method definition that was created
            by `_populate_column_methods()`, where every data column in
            the `.standard` table gets its own function with basic docstring.
            These definitions should also show up with `pop.<tab>` in jupyter.

            3) We use this `__getattr__` function. It will be called only if
            no other definition anywhere overrides it. Practically, this is
            mostly used for retrieving obvious plotting keywords, which
            might be hidden inside the `._plotkw` internal dictionary, or
            to get uncertainties when requesting `.x_uncertainty()` or
            `.x_uncertainty_lowerupper()`.

        Parameters
        ----------
        key : str
            The attribute we're trying to get.

        Returns
        -------
        value : any
            The attribute we're trying to get.
        """

        # do a quick check that something essential isn't missing
        if key in ["label", "_plotkw", "standard"]:
            raise RuntimeError(
                f"""
                Yikes! It looks like `.{key}` isn't defined for this Population. 
                Ideally, this error should never been seen, but if it does, something's 
                gone dreadfully wrong. 
                """
            )

        # try to get a plotkw from this pop, from the plotting defaults, from None
        try:
            assert key in allowed_plotkw
            return self._plotkw.get(key, default_plotkw[key])
        except (AssertionError, KeyError):
            pass

        # check if we're asking for a quantity explicitly from the table
        if key.endswith("_from_table"):
            quantity_key = key.split("_from_table")[0]

            def f(distribution=False, **kw):
                return self.get_values_from_table(key=quantity_key, distribution=distribution)

            f.__docstring__ = f"""
            A function to return the column '.{quantity_key}'. 

            Parameters 
            ----------
            distribution : bool 
                Should this return a Distribution, for uncertainty propagation?

            Returns
            -------
            value : np.array, u.Quantity, astropy.uncertainty.Distribution
                The values, drawn from the table. 
            """
            return f

        # check if we're asking for an uncertainty
        if key.endswith("_uncertainty_lowerupper"):
            quantity_key = key.split("_uncertainty_lowerupper")[0]

            def f(**kw):
                return self.get_uncertainty_lowerupper(key=quantity_key, **kw)

            f.__docstring__ = f"""
            A function to return the two-sided uncertainty on '.{quantity_key}'. 

            Returns
            -------
            lower : np.array, u.Quantity
                The magnitude of the lower uncertainties (x_-lower^+upper)
            upper : np.array, u.Quantity
                The magnitude of the upper uncertainties (x_-lower^+upper)              
            """
            return f

        if key.endswith("_uncertainty"):
            quantity_key = key.split("_uncertainty")[0]

            def f(**kw):
                return self.get_uncertainty(key=quantity_key, **kw)

            f.__docstring__ = f"""
            A function to return the symmetric uncertainty on '.{quantity_key}'. 

            Returns
            -------
            sigma : np.array, u.Quantity
                The magnitude of the uncertainties, assumed to be symmetric.
            """
            return f

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

    def get(self, key, **kw):
        """
        Return an array property for a Population.

        This returns an array, not a function. For example, stellar
        radius might be retrieved either as `.stellar_radius()`
        or as `.get('stellar_radius')`. Quantities that can take
        keyword arguments can pass those keywords as in
        `.stellar_brightness(wavelength=5*u.micron)` or
        `.get('stellar_brightness', wavelength=5*u.micron)`.


        Parameters
        ----------
        key : str
            The name of the quantity we're trying to retrieve.

        kw : dict
            All additional keywords will be passed to the
            method that retrieves the quantities. One common
            (internal) use might be to pass `distribution=True`,
            to be used for uncertainty propagation.

        Returns
        -------
        value : Quantity
            The values, as an array with units.
        """
        return getattr(self, key)(**kw)

    def get_uncertainty_lowerupper(self, key, **kw):
        """
        Return two arrays of lower and upper uncertainties.

        For table column quantities with uncertainties, this should
        return uncertainties direct from `_uncertainty` or
        `_uncertainty_lower` + `_uncertainty_upper` columns in the table.
        This returns an array, not a function.

        For derived quantities, this uncertainty will be estimated from
        central 68% confidence interval of samples from the calculated
        distribution of the quantity.

        Parameters
        ----------
        key : str
            The quantity for which we want lower + upper uncertaintes.
        kw : dict
            All other keywords will be passed to `.get`.

        Returns
        -------
        lower : np.array
            The magnitude of the lower uncertainties (x_{-lower}^{+upper})
        upper : np.array
            The magnitude of the upper uncertainties (x_{-lower}^{+upper})
        """

        # first try for uncertainties direct from table
        try:
            sigma_lower, sigma_upper = self.get_uncertainty_lowerupper_from_table(key)
            return sigma_lower, sigma_upper
        except KeyError:
            mu = self.get(key, **kw)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=RuntimeWarning)
                f = self.targeted_fractional_uncertainty_precision
                total_number_of_samples = 1/f**2
                number_of_iterations = int(np.maximum(np.ceil(total_number_of_samples/self._number_of_uncertainty_samples), 1))

                sigma_lowers, sigma_uppers = [], []
                for i in range(number_of_iterations):
                    d = self.get(key, distribution=True, **kw)
                    # calculate uncertainties from percentiles if possible...
                    if isinstance(d, Distribution):
                        lower, upper = d.pdf_percentiles(
                            100
                            * np.array(
                                [
                                    0.5 - gaussian_central_1sigma / 2,
                                    0.5 + gaussian_central_1sigma / 2,
                                ]
                            )
                        )
                        sigma_lower, sigma_upper = mu - lower, upper - mu
                    # ...otherwise just force the uncertainties to zero
                    else:
                        sigma_lower, sigma_upper = 0 * mu, 0 * mu
                    sigma_lowers.append(sigma_lower)
                    sigma_uppers.append(sigma_upper)
                average_sigma_lower = np.mean(sigma_lowers, axis=0)
                average_sigma_upper = np.mean(sigma_uppers, axis=0)

                return average_sigma_lower, average_sigma_upper

    def get_uncertainty(self, key, **kw):
        """
        Return an array of symmetric uncertainties on a quantity.

        This (very crudely) averages the lower and upper uncertainties
        to make a symmetric errorbar. For quantities with nearly symmetric
        errors this should be totally fine, but for ones with wildly
        asymmetric uncertainties, use `.get_uncertainty_lowerupper()`.
        This returns an array, not a function.

        Parameters
        ----------
        key : str
            The quantity for which we want uncertaintes.

        Returns
        -------
        uncertainty : Quantity
            The (symmetrically-averaged) uncertainties, as an array with units.
        """

        sigma_lower, sigma_upper = self.get_uncertainty_lowerupper(key, **kw)
        return 0.5 * (sigma_lower + sigma_upper)

    def get_fractional_uncertainty(self, key, **kw):
        """
        Return an array of estimates of the fractional uncertainty on a quantity.

        This little wrapper calculates the ratio sigma_x/x, to make it
        easier to select subsets on the basis of something like
        "it's been measured to better than 20% precision".

        Parameters
        ----------
        key : str
            The quantity for which we want fractional uncertainty.
        **kw : dict
            Additional keywords will be passed to the

        Returns
        -------
        fractional_uncertainty : Quantity
            The (symmetrically-averaged) uncertainties, as an array with units.
        """
        x = self.get(key, **kw)
        sigma_x = self.get_uncertainty(key, **kw)
        fractional_uncertainty = sigma_x / x
        return fractional_uncertainty

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
                    n = sum(np.atleast_1d(self.standard[k] != ""))
            self._speak(f"{k:>25} | {n:4}/{N} rows = {n/N:4.0%} are not empty")

    def _find_index(self, name):
        """
        Return index of a particular planet in the population.

        ??? = maybe this could/should be replaced with some table cleverness?
        """

        return np.array([clean(name) in clean(x) for x in self.name()]).nonzero()[0]

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
        d = {c: getattr(self, c)() for c in desired_columns}

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

    def _choose_calculation(
        self,
        methods=[],
        distribution=False,
        how_to_choose="preference",
        visualize=False,
        **kw,
    ):
        """
        Choose values element-wise from multiple functions.

        This wrapper helps merge together different calculations
        of the same quantity, choosing sources based on preference,
        whether the calculation is finite, and/or uncertainties.

        Parameters
        ----------
        methods : list of strings
            The names of methods to consider for calculating
            the quantity, in order of preference. The order
            matters most if `how_to_choose=='preference'`
            (see below); if `how_to_choose=='precision'`
            then preference doesn't really matter (except
            in cases where two different options have
            identical fractional uncertainties).
        distribution : bool
            If False, return a simple array of values.
            If True, return an astropy.uncertainty.Distribution,
            which can be used for error propagation.
        how_to_choose : str
            A string describing how to pick values from
            among the different possible calculations.
            Options include:
                'preference' = pick the highest preference
                option that produces a finite value
                'precision' = pick the most precise value,
                base on the calculated mean uncertainty
        visualize : bool
            Should we visualize all the options and what gets chosen?
        **kw : dict
            All other keywords will be passed to the functions
            for calculating quantities. All functions must be able
            to accept the same set of keyword.
        """

        #
        errors_indicating_missing_data = (AttributeError, KeyError)

        # construct a list of arrays of values
        value_estimates = []
        for m in methods:
            try:
                this = self.get(m, distribution=distribution, **kw)
                value_estimates.append(this)
            except errors_indicating_missing_data:
                pass

        if how_to_choose == "preference":
            # loop through options, picking the first finite value
            for i, v in enumerate(value_estimates):
                if i == 0:
                    values = v
                isnt_good_yet = np.isfinite(values) == False
                values[isnt_good_yet] = v[isnt_good_yet]
        elif how_to_choose == "precision":
            # calculate (symmetric) fractional uncertainties for all options
            mu_estimates = []
            for m in methods:
                try:
                    this = self.get(m, distribution=False, **kw)
                    mu_estimates.append(this)
                except errors_indicating_missing_data:
                    pass
            uncertainty_estimates = []
            for m in methods:
                try:
                    this = self.get_uncertainty(m, **kw)
                    uncertainty_estimates.append(this)
                except errors_indicating_missing_data:
                    pass

            fractional_uncertainty_estimates = [
                u / v for v, u in zip(mu_estimates, uncertainty_estimates)
            ]
            # loop through options to select the one with smallest uncertainty
            for i, (v, f) in enumerate(
                zip(value_estimates, fractional_uncertainty_estimates)
            ):
                f[np.isnan(f)] = np.inf
                if i == 0:
                    values = v
                    current_fractional_uncertainty = f
                has_smaller_uncertainty = f < current_fractional_uncertainty
                current_fractional_uncertainty[has_smaller_uncertainty] = f[
                    has_smaller_uncertainty
                ]
                values[has_smaller_uncertainty] = v[has_smaller_uncertainty]
        else:
            raise ValueError(
                f"""
            "{how_to_choose}" is not a valid option choosing from among
            {methods}
            Only "preference" or "precision" are currently allowed. 
            """
            )

        if visualize:
            plt.figure(figsize=(8, 3))
            if distribution:
                for m, v in zip(methods, value_estimates):
                    plt.violinplot(
                        dataset=np.array(v.distribution.T), positions=np.arange(len(v))
                    )
                plt.plot(
                    values.pdf_median(),
                    color="black",
                    marker="o",
                    markerfacecolor="none",
                )
            else:
                for m, v in zip(methods, value_estimates):
                    plt.plot(v, alpha=0.5, label=m, marker=".")
                plt.plot(
                    values,
                    color="black",
                    marker="o",
                    markersize=10,
                    markerfacecolor="none",
                    linewidth=0,
                )
            plt.legend(frameon=False)

        return values

    from .calculations.planetary import (
        semimajoraxis_from_period,
        semimajoraxis_from_transit_scaled_semimajoraxis,
        semimajoraxis,
        scaled_semimajoraxis_from_semimajoraxis,
        scaled_semimajoraxis,
        eccentricity,
        argument_of_periastron,
        transit_impact_parameter_from_inclination,
        transit_impact_parameter,
        insolation,
        relative_insolation,
        log_relative_insolation,
        relative_cumulative_xuv_insolation,
        teq,
        planet_luminosity,
        transit_depth_from_radii,
        transit_depth,
        scaled_radius_from_radii,
        scaled_radius,
        transit_duration_from_orbit,
        transit_duration,
        mass_estimated_from_radius,
        radius_estimated_from_mass,
        kludge_mass,
        kludge_radius,
        kludge_stellar_age,
        surface_gravity,
        density,
        escape_velocity,
        orbital_velocity,
        impact_velocity,
        escape_parameter,
        scale_height,
    )

    from .calculations.stellar import (
        stellar_luminosity_from_radius_and_teff,
        stellar_luminosity,
        distance_modulus,
    )

    from .calculations.observability import (
        angular_separation,
        imaging_contrast,
        transmission_signal,
        emission_signal,
        reflection_signal,
        stellar_brightness,
        stellar_brightness_in_telescope_units,
        depth_uncertainty,
        _get_noise_and_unit,
        depth_snr,
        emission_snr,
        reflection_snr,
        transmission_snr,
    )
