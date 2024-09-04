# general class for exoplanet populations
from ..imports import *
from ..telescopes import *
from ..models import *

import string

basic_columns = ["name", "hostname", "ra", "dec", "distance"]

transit_columns = [
    "period",
    "semimajoraxis",
    "eccentricity",
    "omega",
    "inclination",
    "transit_midpoint",
    "transit_duration",
    "transit_depth",
    "stellar_teff",
    "stellar_mass",
    "stellar_radius",
    "radius",
    "mass",
    "transit_ar",
    "transit_b",
]

calculated_columns = [
    "a_over_rs",
    "b",
    "insolation",
    "relative_insolation",
    "log_relative_insolation",
    "teq",
    "planet_luminosity",
    "density",
    "surface_gravity",
    "distance_modulus",
    "escape_velocity",
    "escape_parameter",
    "angular_separation",
    "imaging_contrast",
    "stellar_luminosity",
]


table_columns = basic_columns + transit_columns
attribute_columns = table_columns + calculated_columns


method_columns = [
    "scale_height",
    "transmission_signal",
    "transmission_snr",
    "emission_signal",
    "emission_snr",
    "reflection_signal",
    "reflection_snr",
    "stellar_brightness",
    "stellar_brightness_in_telescope_units",
    "depth_uncertainty",
]

desired_columns = [
    "mass_uncertainty_upper",
    "mass_uncertainty_lower",
    "radius_uncertainty_upper",
    "radius_uncertainty_lower",
    "distance_uncertainty_upper",
    "distance_uncertainty_lower",
]

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
allowed_plotkw += ["s", "c", "cmap", "norm", "vmin", "vmax" "outlined", "filled"]


class Population(Talker):
    """
    Populations of astronomical objects might contain
        Exoplanets (planets with host stars),
        Solar System objects (major/minor planets),
        Stars (single or multiples),
    or more!
    """

    # kludge?
    _pithy = True

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
            self.plotkw = plotkw

        elif isinstance(standard, str):
            filename = standard
            self.standard = QTable(ascii.read(filename))
            self.label = self.standard.meta["label"]
            self.plotkw = self.standard.meta["plotkw"]

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

    def _list_table_indices(self):
        """
        Return a list of keys being used as table indices.
        """
        return [x.columns[0].name for x in self.standard.indices]

    def _make_sure_index_exists(self, k):
        """
        Add a key as an index, but don't add it twice.
        """
        if k not in self._list_table_indices():
            self.standard.add_index(k)

    @property
    def fileprefix(self):
        """
        Define a fileprefix for this population, to be used
        for setting the filename of the standardized population.
        """
        return clean(self.label)

    @property
    def standardized_data_path(self):
        """
        Define the filepath for the standardized table.
        """
        return os.path.join(directories["data"], f"standardized-{self.fileprefix}.txt")

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
            Whether

        Examples
        --------
        >>> from exoatlas import *
        >>> one = Exoplanets()['GJ1132b']
        >>> one.save('pop.txt')
        >>> the_same_one = Population('pop.txt')
        """

        # be a little fussy about overwriting automatic filenames
        if filename == None:
            filename = f"exoatlas-population-{self.label}.ecsv"
            overwrite = False
            if os.path.exists(filename):
                warnings.warn(
                    f"{filename} will not be overwritten unless you explicitly provide a filename."
                )

                # save the table as an ascii table for humans to read
        to_save = copy.deepcopy(self.standard)
        to_save.meta["label"] = self.label
        to_save.meta["plotkw"] = self.plotkw

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

        Parameters
        ----------
        x : str
            The key by which to sort the table.
        reverse : bool
            Whether to reverse the sort order.
                `reverse = False` means low to high
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

    def remove_by_key(self, other, key="tidyname"):
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
        return self.remove_by_key(other)

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
            subset = type(self)(standard=self.standard[key], label=label, **self.plotkw)

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
        return type(self)(standard=subset, label=label, **self.plotkw)

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
        return type(self)(standard=subset, label=label, **self.plotkw)

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
        RA and Dec. This will return all plannets
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
        new_population = type(self)(standard=subset, label=label, **self.plotkw)

        # choose what to return
        if return_indices:
            i_from_original_coordinates = idx[i_match]
            return new_population, i_from_original_coordinates
        else:
            return new_population

    def __getattr__(self, key):
        """
        If an attribute/method isn't defined for a population,
        look for it as a column of the standardized table.

        For example, `population.stellar_radius` will try to
        access `population.standard['stellar_radius']`.

        Parameters
        ----------
        key : str
            The attribute we're trying to get.
        """
        if key in ["label", "plotkw", "standard"]:
            raise RuntimeError(f"Yikes! {key}")
        try:
            # extract the column from the standardized table
            return self.standard[key]
        except KeyError:
            # try to get a plotkw from this pop, from the plotting defaults, from None
            try:
                assert key in allowed_plotkw
                return self.plotkw.get(key, default_plotkw[key])
            except (AssertionError, KeyError):
                raise AttributeError(
                    f"""
                Alas, there seems to be no way to find `.{key}`
                as an attribute or propetry of {self}.
                """
                )  # AtlasError

    def __setattr__(self, key, value):
        """
        Define what happens when we try to set an attribute via `pop.attr = x`.
        If the keyword is a pre-defined "plotting" keyword in `allowed_plotkw`,
        then we should save it in a special `plotkw` dictionary. Otherwise,
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
            self.plotkw[key] = value
        else:
            # otherwise, store attributes as normal for objects
            self.__dict__[key] = value

    def __repr__(self):
        """
        How should this object appear as a repr/str?
        """
        return f"<{self.label} | {self.n} elements>"

    def get(self, key):
        """
        Return an array of values for a column.
        """
        return getattr(self, key)

    def get_uncertainty(self, key):
        """
        Return an array of symmetric uncertainties on a column.

        Parameters
        ----------
        key : str
            The column for which we want errors.
        """

        # first try for an `uncertainty_{key}` column
        try:
            return self.__getattr__(f"{key}_uncertainty")
        except (
            KeyError,
            AssertionError,
            AtlasError,
            AttributeError,
        ):  # is including AttributeError a kludge?
            # this can be removed after debugging
            self.speak(f'no symmetric uncertainties found for "{key}"')

        # then try for crudely averaging asymmetric uncertainties
        try:
            lower = self.__getattr__(f"{key}_uncertainty_lower")
            upper = self.__getattr__(f"{key}_uncertainty_upper")
            avg = 0.5 * (np.abs(lower) + np.abs(upper))
            return avg
        except (KeyError, AssertionError, AtlasError, AttributeError):
            # this can be removed after debugging
            self.speak(f'no asymmetric uncertainties found for "{key}"')

        # then give up and return nans
        return np.nan * self.standard[key]

    def get_uncertainty_lowerupper(self, key):
        """
        Return two arrays of lower and upper uncertainties on a column.

        Parameters
        ----------
        key : str
            The column for which we want errors.

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
            self.speak(f'no asymmetric uncertainties found for "{key}"')

        # first try for an `uncertainty_{key}` column
        try:
            sym = self.__getattr__(f"{key}_uncertainty")
            return np.abs(sym), np.abs(sym)
        except (KeyError, AssertionError, AttributeError):
            # this can be removed after debugging
            self.speak(f'no symmetric uncertainties found for "{key}"')

        # then give up and return nans
        unc = np.nan * self.__getattr__(key)
        return unc, unc

    def single(self, name):
        """
        Create a subpopulation of a single planet.
        """

        # create a subset of the standardized table
        subset = self.standard.loc[name]

        # create a new object, from this subset
        return type(self)(standard=subset, label=name, **self.plotkw)

    def validate_columns(self):
        """
        Make sure this standardized table has all the necessary columns.
        Summarize the amount of good data in each.
        """

        N = len(self.standard)
        for k in table_columns:
            try:
                n = sum(self.standard[k].mask == False)
            except AttributeError:
                try:
                    n = sum(np.isfinite(self.standard[k]))
                except TypeError:
                    n = sum(self.standard[k] != "")
            self.speak(f"{k:>25} | {n:4}/{N} rows = {n/N:4.0%} are not empty")

    def find(self, name):
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

    @property
    def n(self):
        """
        How many planets are in this population?
        """
        return len(self.standard)

    def __len__(self):
        """
        How many planets are in this population?
        """
        return len(self.standard)

    @property
    def semimajor_axis(self):
        """
        Have a safe way to calculate the semimajor axis of planets,
        that fills in gaps as necessary. Basic strategy:

            First from table.
            Then from NVK3L.
            Then from a/R*.

        """

        # pull out the actual values from the table
        a = self.standard["semimajoraxis"].copy()

        # try to replace bad ones with NVK3L
        bad = np.isfinite(a) == False
        self.speak(f"{sum(bad)}/{self.n} semimajoraxes are missing")

        # calculate from the period and the stellar mass
        P = self.period[bad]
        M = self.stellar_mass[bad]
        G = con.G
        a[bad] = ((G * M * P**2 / 4 / np.pi**2) ** (1 / 3)).to("AU")

        # replace those that are still bad with the a/R*
        stillbad = np.isfinite(a) == False
        self.speak(f"{sum(stillbad)}/{self.n} are still missing after NVK3L")
        # (pull from table to avoid potential for recursion)
        try:
            a_over_rs = self.standard["transit_ar"][stillbad]
            rs = self.standard["stellar_radius"][stillbad]
            a[stillbad] = a_over_rs * rs
        except KeyError:
            pass

        return a

    @property
    def angular_separation(self):
        """
        Calculate the angular separation,
        simply as theta = a/D
        """

        a = self.semimajor_axis
        D = self.distance

        theta = np.arctan(a / D).to(u.arcsec)

        return theta

    @property
    def imaging_contrast(self):
        """
        What is the reflected light eclipse depth,
        for an albedo of 100%?

        But use a kludged radius
        """
        return 0.25 * (self.kludge_radius / self.semimajor_axis).decompose() ** 2

    @property
    def a_over_rs(self):
        """
        Have a safe way to calculate the scaled semimajor axis of planets,
        that fills in gaps as necessary. Basic strategy:

            First from table, mostly derived from transit.
            Then from the semimajor axis.
        """

        # pull out the values from the table
        a_over_rs = self.standard["transit_ar"].copy()

        # try to replace bad ones with NVK3L
        bad = np.isfinite(a_over_rs) == False
        self.speak(f"{sum(bad)}/{self.n} values for a/R* are missing")

        a = self.semimajor_axis[bad]
        R = self.stellar_radius[bad]
        a_over_rs[bad] = a / R

        stillbad = np.isfinite(a_over_rs) == False
        self.speak(f"{sum(stillbad)}/{self.n} are still missing after a and R*")

        return a_over_rs

    @property
    def stellar_luminosity(self):
        T = self.stellar_teff
        R = self.stellar_radius
        sigma = con.sigma_sb
        return (4 * np.pi * R**2 * sigma * T**4).to(u.Lsun)

    @property
    def e(self):
        """
        FIXME -- assumes are missing eccentricities are 0!
        """

        # pull out the actual values from the table
        e = self.standard["eccentricity"].copy()

        # try to replace bad ones with NVK3L
        bad = np.isfinite(e) == False
        self.speak(f"{sum(bad)}/{self.n} eccentricities are missing")
        self.speak(f"assuming they are all zero")
        e[bad] = 0

        return e

    @property
    def omega(self):
        """
        (FIXME! we need better longitudes of periastron)
        """

        # pull out the actual values from the table
        omega = self.standard["omega"].copy()

        # try to replace bad ones with NVK3L
        bad = np.isfinite(omega) == False
        self.speak(f"{sum(bad)}/{self.n} longitudes of periastron are missing")
        e_zero = self.e == 0
        self.speak(f"{sum(e_zero)} have eccentricities assumed to be 0")
        omega[e_zero] = 0 * u.deg

        return omega

    @property
    def b(self):
        """
        Transit impact parameter.
        (FIXME! split this into transit and occultation)
        """

        # pull out the actual values from the table
        b = self.standard["transit_b"].copy()

        # try to replace bad ones with NVK3L
        bad = np.isfinite(b) == False
        self.speak(f"{sum(bad)}/{self.n} impact parameters are missing")

        # calculate from the period and the stellar mass
        a_over_rs = self.a_over_rs[bad]
        i = self.standard["inclination"][bad]
        e = self.e[bad]
        omega = self.omega[bad]
        b[bad] = a_over_rs * np.cos(i) * ((1 - e**2) / (1 + e * np.sin(omega)))

        # report those that are still bad
        stillbad = np.isfinite(b) == False
        self.speak(f"{sum(stillbad)}/{self.n} are still missing after using i")

        return b

    # the 1360 W/m^2 that Earth receives from the Sun
    earth_insolation = (1 * u.Lsun / 4 / np.pi / u.AU**2).to(u.W / u.m**2)

    @property
    def insolation(self):
        """
        The insolation the planet receives, in W/m^2.
        """

        # calculate the average insolation the planet receives
        insolation = self.stellar_luminosity / 4 / np.pi / self.semimajor_axis**2
        return insolation.to(u.W / u.m**2)

    @property
    def relative_insolation(self):
        """
        The insolation the planet receives, relative to Earth.
        """
        return self.insolation / self.earth_insolation

    @property
    def log_relative_insolation(self):
        return np.log10(self.relative_insolation)

    @property
    def relative_cumulative_xuv_insolation(self):
        """
        A *very very very* approximate estimate for the cumulative XUV flux (J/m**2)
        felt by a planet over its lifetime. It comes from Zahnle + Catling (2017) Equation 27,
        where they say they did integrals over the Lammer et al. (2009) XUV flux relations.
        This effectively assumes that the early times dominate, so the time integral doesn't
        depend (?!?) on the age of the system. It's very rough!
        """
        xuv_proxy = (self.stellar_luminosity / u.Lsun) ** -0.6
        return self.relative_insolation * xuv_proxy

    @property
    def teq(self):
        """
        The equilibrium temperature of the planet.
        """
        f = self.insolation
        sigma = con.sigma_sb
        A = 1
        return ((f * A / 4 / sigma) ** (1 / 4)).to(u.K)

    @property
    def planet_luminosity(self):
        """
        The bolometric luminosity of the planet (assuming zero albedo).
        """
        return (self.teq**4 * con.sigma_sb * 4 * np.pi * self.radius**2).to(u.W)

    @property
    def transit_depth(self):
        """
        The depth of the transit
        (FIXME, clarify if this is 1.5-3.5 or what)
        """

        # pull out the actual values from the table
        d = self.standard["transit_depth"].copy()

        # try to replace bad ones with NVK3L
        bad = np.isfinite(d) == False
        self.speak(f"{sum(bad)}/{self.n} transit depths are missing")

        Rp = self.radius[bad]
        Rs = self.stellar_radius[bad]

        d[bad] = (Rp / Rs).decompose() ** 2

        # report those that are still bad
        stillbad = np.isfinite(d) == False
        self.speak(f"{sum(stillbad)}/{self.n} are still missing after Rp/Rs")

        return d

    @property
    def transit_duration(self):
        """
        The duration of the transit
        (FIXME, clarify if this is 1.5-3.5 or what)
        """

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

            # pull out the actual values from the table
            d = self.standard["transit_duration"].copy()

            # try to replace bad ones with NVK3L
            bad = np.isfinite(d) == False
            self.speak(f"{sum(bad)}/{self.n} transit durations are missing")

            P = self.period[bad]
            a_over_rs = self.a_over_rs[bad]
            b = self.b[bad]

            T0 = P / np.pi / a_over_rs
            T = T0 * np.sqrt(1 - b**2)

            e = self.e[bad]
            omega = self.omega[bad]
            factor = np.sqrt(1 - e**2) / (1 + e * np.sin(omega))

            d[bad] = (T * factor).to(u.day)

            # report those that are still bad
            stillbad = np.isfinite(d) == False
            self.speak(f"{sum(stillbad)}/{self.n} are still missing after P, a/R*, b")

            return d

    @property
    def kludge_mass(self):
        """
        Have a safe way to calculate the mass of planets,
        that fills in gaps as necessary. Basic strategy:

            First from table.
            Then from msini.
        """

        # pull out the actual values from the table
        M = self.standard["mass"].copy()

        # try to replace bad ones with NVK3L
        bad = np.isfinite(M) == False
        self.speak(f"{sum(bad)}/{self.n} masses are missing")

        # estimate from the msini
        try:
            M[bad] = self.msini[bad]
        except (KeyError, AssertionError, AtlasError, AttributeError):
            pass

        # replace those that are still bad with the a/R*
        stillbad = np.isfinite(M) == False
        self.speak(f"{sum(stillbad)}/{self.n} are still missing after msini")

        return M

    @property
    def kludge_radius(self):
        """
        Have a safe way to calculate the radii of planets,
        that fills in gaps as necessary. Basic strategy:

            First from table.
            Then from mass, via Chen & Kipping (2017).
        """

        # pull out the actual values from the table
        R = self.standard["radius"].copy()

        # try to replace bad ones with NVK3L
        bad = np.isfinite(R) == False
        self.speak(f"{sum(bad)}/{self.n} radii are missing")

        # estimate from Chen and Kipping
        try:
            M = self.kludge_mass
            R[bad] = estimate_radius(M[bad])
        except (KeyError, AssertionError, AtlasError, AttributeError):
            pass

        # replace those that are still bad with the a/R*
        stillbad = np.isfinite(R) == False
        self.speak(
            f"{sum(stillbad)}/{self.n} are still missing after Chen & Kipping (2017)"
        )

        return R

    @property
    def kludge_age(self):
        """
        Have a safe way to calculate the age of planets,
        that fills in gaps as necessary. Basic strategy:

            First from table.
            Then assume 5 Gyr.
        """

        # pull out the actual values from the table
        age = self.standard["stellar_age"].copy()

        # try to replace bad ones with NVK3L
        bad = np.isfinite(age) == False
        self.speak(f"{sum(bad)}/{self.n} ages are missing")

        # estimate from the msini
        try:
            age[bad] = 5 * u.Gyr
        except (KeyError, AssertionError, AtlasError, AttributeError):
            pass

        # replace those that are still bad with the a/R*
        stillbad = np.isfinite(age) == False
        self.speak(
            f"{sum(stillbad)}/{self.n} are still missing after blindly assuming 5Gyr for missing ages"
        )

        return age

    @property
    def surface_gravity(self):
        """
        (FIXME) -- make an assumption for planets without masses
        """

        G = con.G
        M = self.mass
        R = self.radius

        g = (G * M / R**2).to("m/s**2")
        return g

    @property
    def density(self):
        """
        The density of the planet.
        """
        mass = self.mass
        volume = 4 / 3 * np.pi * (self.radius) ** 3
        return (mass / volume).to("g/cm**3")

    @property
    def escape_velocity(self):
        """
        The escape velocity of the planet.
        """
        G = con.G
        M = self.mass
        R = self.radius
        return np.sqrt(2 * G * M / R).to("km/s")

    @property
    def orbital_velocity(self):
        return (2 * np.pi * self.semimajor_axis / self.period).to("km/s")

    @property
    def impact_velocity(self):
        return np.sqrt(self.orbital_velocity**2 + self.escape_velocity**2)

    @property
    def escape_parameter(self):
        """
        The Jeans atmospheric escape parameter for atomic hydrogen,
        at the equilibrium temperature of the planet.
        """
        k = con.k_B
        T = self.teq
        mu = 1
        m_p = con.m_p
        G = con.G
        M = self.mass
        R = self.radius

        e_thermal = k * T
        e_grav = G * M * m_p / R
        return (e_grav / e_thermal).decompose()

    @property
    def distance_modulus(self):
        """
        The distance modulus to the system, in magnitudes.
        """
        mu = 5 * np.log10(self.distance / (10 * u.pc))
        return mu

    def scale_height(self, mu=2.32):
        """
        The scale height of the atmosphere, at equilibrium temperature.
        """
        k = con.k_B
        T = self.teq
        m_p = con.m_p
        g = self.surface_gravity
        return (k * T / mu / m_p / g).to("km")

    def transmission_signal(self, mu=2.32, threshold=2):
        """
        What is the transit depth of 1 scale height of an
        atmosphere transiting in front of the star.

        Parameters
        ----------
        mu : float
            Mean molecular weight (default 2.2 for H/He)
        threshold : float
            By how many sigma must the planet mass be detected?

        """
        with np.errstate(invalid="ignore"):
            H = self.scale_height(mu)
            Rp = self.radius
            Rs = self.stellar_radius
            depth = (2 * H * Rp / Rs**2).decompose()

            dlnm = self.get_uncertainty("mass") / self.mass
            bad = dlnm > 1 / threshold
            depth[bad] = np.nan
            return depth

    def reflection_signal(self, albedo=0.1):
        """
        What is the reflected light eclipse depth,
        for an albedo of 100%?
        """
        return albedo * 0.25 * (self.radius / self.semimajor_axis).decompose() ** 2

    def emission_signal(self, wavelength=5 * u.micron):
        """
        What is the thermal emission eclipse depth,
        assuming Planck spectra for both star and planet?

        This calculation assumes a Bond albedo of 0
        and that heat is uniformly distributed over the planet.

        Parameters
        ----------
        wavelength : astropy.unit.Quantity
            The wavelength at which it should be calculated.
        """

        # create thermal emission sources for both star and planet
        import rainbowconnection as rc

        star = rc.Thermal(teff=self.stellar_teff, radius=self.stellar_radius)
        planet = rc.Thermal(teff=self.teq, radius=self.radius)

        # calculate the depth as the luminosity ratio
        depths = planet.spectrum(wavelength) / star.spectrum(wavelength)

        return depths

    def stellar_brightness(self, wavelength=5 * u.micron):
        """
        How many photons/s/m^2/micron do we receive from the star?

        This is calculated from the distance, radius, and
        stellar effective temperature of the stars.

        (It could be potentially be improved with PHOENIX
        model grids and/or cleverness with photometry.)

        Parameters
        ----------
        wavelength : astropy.unit.Quantity
            The wavelength at which it should be calculated.
        """

        # import some tools for easy cartoon spectra
        import rainbowconnection as rc

        # create source with right temperature, size, distance
        teff, radius = self.stellar_teff, self.stellar_radius
        star = rc.Thermal(teff=teff, radius=radius).at(self.distance)

        # calculate the energy flux
        flux_in_energy = star.spectrum(wavelength)

        # convert to photon flux
        photon_energy = con.h * con.c / wavelength / u.ph
        flux_in_photons = flux_in_energy / photon_energy

        # return the
        return flux_in_photons.to("ph s^-1 m^-2 micron^-1")

    def stellar_brightness_in_telescope_units(self, telescope_name="JWST", **kw):
        """
        The stellar brightness, converted to telescope units.

        Parameters
        ----------
        telescope_name : str
            The name of the telescope.

        wavelength : astropy.unit.Quantity
            The wavelength at which it should be calculated.

        R : float
            The spectral resolution at which the
            telescope will bin wavelengths.

        dt : astropy.units.quantity.Quantity
            The time over which the telescope exposes.
        """

        # what counts as 1 "telescope unit" (e.g. JWST at R=20 at 5 microns for 1 hour)
        telescope_unit = define_telescope_unit_by_name(telescope_name, **kw)

        # what's the photon flux (photons/m**2/s)
        flux_in_photons = self.stellar_brightness(telescope_unit.wavelength)

        # quote the brightness as (for example) gigaphotons/JWST at R=20 at 5 microns in 1 hour
        unit = lotsofphotons_unit / telescope_unit
        return flux_in_photons.to(unit)

    def depth_uncertainty(
        self, telescope_name="JWST", per_transit=False, dt=1 * u.hour, **kw
    ):
        """
        What is the transit/eclipse depth uncertainty
        with a particular telescope
        at a particular wavelength
        at a particular resolution?

        By default, this will be calculated for one transit.
        Optionally, it can be calculated for a given amount of time instead.

        Parameters
        ----------
        telescope_name : str
            The name of the telescope.

        per_transit : bool
            If True, calculate the depth uncertainty for one transit.
            If False, calculate the depth uncertainty for a certain amount
            of in-transit time. You likely want to specify `dt` as a
            keyword argument to set that amount of in-transit time.
            In either case, an out-of-transit baseline equal to the
            total in-transit time will be assumed. This means the actual
            time cost will be twice the transit duration or `dt` chosen,
            and the depth uncertainty will be a factor sqrt(2) larger
            than the pure photon noise binned to the relevant timescale.

        wavelength : astropy.unit.Quantity
            The wavelength at which it should be calculated.

        R : float
            The spectral resolution at which the
            telescope will bin wavelengths.

        dt : astropy.units.quantity.Quantity
            The time over which the telescope exposes. If `per_transit=True`,
            this will be ignored. Otherwise, it will set the total amount
            of in-transit time observed, assuming that an equal amount of
            time will *also* be observed out of transit.
        """

        # what counts as 1 "telescope unit" (e.g. JWST at R=20 at 5 microns for 1 hour)
        telescope_unit = define_telescope_unit_by_name(telescope_name, dt=dt, **kw)

        # what's the photon flux (photons/m**2/s)
        flux_in_photons = self.stellar_brightness(telescope_unit.wavelength)

        # what's the total collecting power?
        if per_transit:
            ratio_of_collecting_time = self.transit_duration / dt
        else:
            ratio_of_collecting_time = 1
        collecting_power = 1 * telescope_unit * ratio_of_collecting_time

        # what's the total number of photons collected during transit
        N = (flux_in_photons * collecting_power).to(u.ph).value

        # what's the flux uncertainty on the time scale of one transit?
        sigma = 1 / np.sqrt(N)

        # inflate by a factor of sqrt(2) for equal out-of-transit
        oot = np.sqrt(2)
        sigma_depth = sigma * oot

        return sigma_depth

    def _get_noise_and_unit(self, telescope_name="JWST", per_transit=False, **kw):
        """
        Tiny helper to get the noise and the telescope_unit
        for a telescope observation of a planet.
        """

        # figure out the noise
        noise = self.depth_uncertainty(
            telescope_name=telescope_name, per_transit=per_transit, **kw
        )

        # create a telescope unit (mostly to get a default wavelength)
        telescope_unit = define_telescope_unit_by_name(telescope_name, **kw)

        return noise, telescope_unit

    def depth_snr(self, telescope_name="JWST", **kw):
        """
        What's the approximate S/N for the detection of the planet's transit?
        """

        noise, telescope_unit = self._get_noise_and_unit(
            telescope_name=telescope_name, **kw
        )
        signal = self.transit_depth
        return signal / noise

    def emission_snr(self, telescope_name="JWST", **kw):
        """
        What's the approximate S/N for the detection of the
        thermal emission eclipse of a planet?
        """

        noise, telescope_unit = self._get_noise_and_unit(
            telescope_name=telescope_name, **kw
        )
        signal = self.emission_signal(wavelength=telescope_unit.wavelength)
        return signal / noise

    def reflection_snr(self, telescope_name="JWST", albedo=1, **kw):
        """
        What's the approximate S/N for the detection of the
        reflected light eclipse of a planet?
        """

        noise, telescope_unit = self._get_noise_and_unit(
            telescope_name=telescope_name, **kw
        )
        signal = self.reflection_signal(albedo=albedo)
        return signal / noise

    def transmission_snr(self, telescope_name="JWST", mu=2.32, threshold=2, **kw):
        """
        What's the approximate S/N for the detection of the
        reflected light eclipse of a planet?
        """

        noise, telescope_unit = self._get_noise_and_unit(
            telescope_name=telescope_name, **kw
        )
        signal = self.transmission_signal(mu=mu, threshold=threshold)
        return signal / noise

    def scatter(
        self, xname, yname, c=None, s=None, names=True, xlog=True, ylog=True, **kw
    ):
        """
        Quick tool to plot one parameter against another.
        """
        plt.ion()
        x, y = self.__getattr__(xname), self.__getattr__(yname)
        try:
            self.ax.cla()
        except:
            self.figure = plt.figure("Exoplanet Population")
            self.ax = plt.subplot()

        self.ax.set_xlabel(xname)
        self.ax.set_ylabel(yname)
        self.ax.scatter(x, y, c=c, s=s, **kw)
        if False:
            for i in range(len(x)):
                self.ax.text(x[i], y[i], self.table["NAME"][i])
        if xlog:
            plt.xscale("log")
        if ylog:
            plt.yscale("log")

        plt.draw()

    def thumbtack(self, maxr=1000, dr=100, labels=False):
        """Plot the planets as thumbtacks."""

        def scale(d):
            return np.array(d) ** 1.5

        r = scale(self.distance)
        x, y = r * np.cos(self.ra * np.pi / 180), r * np.sin(self.ra * np.pi / 180)
        plt.ion()
        plt.figure("thumbtacks")

        ax = plt.subplot()
        ax.cla()
        ax.set_aspect("equal")
        theta = np.linspace(0, 2 * np.pi, 1000)
        angle = -90 * np.pi / 180

        gridkw = dict(alpha=0.25, color="green")
        for originalradius in np.arange(dr, maxr * 2, dr):
            radii = scale(originalradius)

            ax.plot(radii * np.cos(theta), radii * np.sin(theta), linewidth=3, **gridkw)
            ax.text(
                radii * np.cos(angle),
                radii * np.sin(angle),
                "{0:.0f} pc".format(originalradius),
                rotation=90 + angle * 180 / np.pi,
                va="bottom",
                ha="center",
                size=13,
                weight="extra bold",
                **gridkw,
            )

        ax.plot(
            x, y, marker="o", alpha=0.5, color="gray", linewidth=0, markeredgewidth=0
        )
        close = (self.name == "WASP-94A b").nonzero()[
            0
        ]  # (self.distance < maxr).nonzero()[0]
        if labels:
            for c in close:
                plt.text(x[c], y[c], self.name[c])
        ax.set_xlim(-scale(maxr), scale(maxr))
        ax.set_ylim(-scale(maxr), scale(maxr))

    def compare(self, x="teq", y="radius", area="depth", color="stellar_radius"):
        xplot = self.__dict__[x]
        yplot = self.__dict__[y]
        sizeplot = self.__dict__[size]
        colorplot = self.__dict__[color]

        maxarea = 1000
        area = self.__dict__[area]
        sizeplot = np.sqrt(area / np.nanmax(area) * maxarea)

        plt.scatter(xplot, yplot, linewidth=0, marker="o", markersize=sizeplot)

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
