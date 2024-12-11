# exoplanet population of all "confirmed" exoplanets from exoplanet archive
from ...imports import *
from ..predefined import PredefinedPopulation
from .exoplanet_downloaders import *

__all__ = ["Exoplanets", "ExoplanetsPSCP", "ExoplanetsPS"]

suffixes = [
    "_reference",
    "",
    "_uncertainty_lower",
    "_uncertainty_upper",
    "_lower_limit",
    "_upper_limit",
]


def parse_reflink(x):
    """
    Convert a 'reflink/refname' from NASA Exoplanet Archive into
    a human-friendly string and a definitely-unique URL.

    Parameters
    ----------
    x : str
        The reflink or refname, which might look like:
        '<a refstr=TORRES_ET_AL__2008 href=https://ui.adsabs.harvard.edu/abs/2008ApJ...677.1324T/abstract target=ref> Torres et al. 2008 </a>'

    Returns
    -------
    name : str
        The human-friendly string for this reference. Be careful, this might not
        necessarily by 100% unique, particularly if folks publish multiple planets
        in the same year.
    """
    try:
        url = x.split("href=")[1].split(" ")[0]
        name = x.split(">")[1].split("<")[0].replace("&amp;", "+").strip()
        return name  # , url
    except AttributeError:
        return ""


class ExoplanetsPSCP(PredefinedPopulation):
    """
    Exoplanets from the NASA Exoplanet Archive.

    This uses the Planetary Systems Composite Parameters table,
    which merges together data from different papers for the
    same planet, which can sometimes result in parameters
    that aren't necessarily consistent with each other.
    """

    label = "ExoplanetsPSCP"

    # define columns to ingest with and without errors/limits
    _raw_columns_with_errors_and_limits = [
        "pl_orbper",
        "pl_orbsmax",
        "pl_rade",
        "pl_bmasse",
        "pl_msinie",
        "pl_dens",
        "pl_orbeccen",
        "pl_orbincl",
        "pl_orblper",
        "pl_tranmid",
        "pl_trandep",
        "pl_trandur",
        "pl_imppar",
        "pl_ratdor",
        "pl_rvamp",
        "pl_projobliq",
        "pl_trueobliq",
        "pl_insol",
        "pl_eqt",
        "st_teff",
        "st_rad",
        "st_mass",
        "st_met",
        "st_lum",
        "st_logg",
        "st_age",
        "st_dens",
        "st_vsin",
        "st_rotp",
        "st_radv",
    ]
    _raw_columns_with_errors = [
        "sy_dist",
        "sy_plx",
        "sy_bmag",
        "sy_vmag",
        "sy_jmag",
        "sy_hmag",
        "sy_kmag",
        "sy_umag",
        "sy_gmag",
        "sy_rmag",
        "sy_imag",
        "sy_zmag",
        "sy_w1mag",
        "sy_w2mag",
        "sy_w3mag",
        "sy_w4mag",
        "sy_gaiamag",
        "sy_tmag",
        "sy_kepmag",
        "sy_icmag",
    ]
    _raw_columns_without_errors = [
        "pl_name",
        "hostname",
        "pl_letter",
        "ra",
        "dec",
        "sy_pmra",
        "sy_pmdec",
        "gaia_id",
        "sy_snum",
        "sy_pnum",
        "discoverymethod",
        "disc_year",
        "disc_refname",
        "disc_facility",
        "rv_flag",
        "pul_flag",
        "ptv_flag",
        "tran_flag",
        "ast_flag",
        "obm_flag",
        "micro_flag",
        "etv_flag",
        "ima_flag",
        "dkin_flag",
        "pl_controv_flag",
        "ttv_flag",
        "st_spectype",
        "pl_nespec",
        "pl_ntranspec",
        "default_flag",
    ]
    _downloader = ExoplanetArchiveDownloader("pscomppars")

    def __init__(self, remake=False, **plotkw):
        """
        Initialize a population of all known exoplanets
        using the NASA Exoplanet Archive `pscomppars` table.

        This generates an Exoplanets population containing
        planets discovered through any method. It tries to
        be clever about loading cached files so that it
        will work quickly and not try to redownload everything
        from scratch every time it's called.

        Parameters
        ----------
        remake : bool
            Should the population be remade from raw files,
            even if a recent version already exists?
        **plotkw : dict
            All other keywords will go toward defining
            default plotting styles, for example like
            'alpha', 'marker', 'zorder', 'color', ...
        """

        # load standard table(s) or ingest from raw data
        PredefinedPopulation.__init__(self, remake=remake, **plotkw)
        # self._add_references_as_indices()

    def _add_references_as_indices(self):
        """
        Add all keys that look like references as table indices for faster lookup.
        """
        for k in self.standard.colnames:
            if "_reference" in k:
                self._make_sure_index_exists(k)

    def _create_standardized(self, raw):
        """
        Create a standardized table to make sure that at
        least a few necessary columns are populated.

        This standardization step aims to make it easier to
        compare exoplanet populations from different sources.
        It translates columns from a raw table into a standard
        set of expected column names, doing some calculations
        if necessary.

        Parameters
        ----------
        raw : astropy.table.QTable
            A raw unstandardized table from `ExoplanetArchiveDownloader`.

        Returns
        -------
        standard : astropy.table.QTable
            A standardize table of exoplanet properties.
        """

        # define the table from which we're deriving everying
        r = raw

        # the new standardized table
        s = QTable()

        def strip_even_if_masked(x):
            """
            Flexibly clean strings by removing trailing spaces.
            """
            try:
                return x.strip()
            except AttributeError:
                return None

        def attach_unit(x, unit=None):
            """
            Flexibly give units to quantities that should have them.
            """
            if isinstance(x[0], str):
                return [strip_even_if_masked(a) for a in x]
            if unit is None:
                return x
            else:
                return x * unit

        def populate_one_or_more_columns(k_new, k_original, unit=None):
            """
            A wrapper to help smoothly grab quantities from raw
            table and ingest them into the standardized table.

            For columns with

            Parameters
            ----------
            k_original : string
                The column name in the original, raw table.
            k_new : string
                The column name in the new, standardized table.
            unit : None, astropy.units.Unit
                A unit to attach to the quantity, if necessary.
            """

            if k_original not in r.colnames:
                warnings.warn(f"⚠️ No {k_original} found!")
                return

            if k_original in self._raw_columns_without_errors:
                # easy, just record the column itself with no error
                s[k_new] = attach_unit(r[k_original], unit)
                print(f"📕 populated {k_original} > {k_new}")

            elif k_original in self._raw_columns_with_errors:
                # record the column itself
                s[k_new] = attach_unit(r[k_original], unit)

                # record the upper and lower errorbars
                s[f"{k_new}_uncertainty_upper"] = attach_unit(
                    r[f"{k_original}err1"], unit
                )
                s[f"{k_new}_uncertainty_lower"] = attach_unit(
                    r[f"{k_original}err2"], unit
                )
                print(f"📏 populated {k_new} and errors with {k_original}")

            elif k_original in self._raw_columns_with_errors_and_limits:
                # from playing with table, I think lim = +1 is upper limit, -1 is lower limit

                # record the upper and lower errorbars
                limit_flag = r[f"{k_original}lim"]

                # populate lower limits
                is_lower = limit_flag == -1
                if np.sum(is_lower) > 0:
                    s[f"{k_new}_lower_limit"] = attach_unit(r[f"{k_original}"])
                    s[f"{k_new}_lower_limit"][is_lower == False] = np.nan

                # populate upper limits
                is_upper = limit_flag == +1
                if np.sum(is_upper) > 0:
                    s[f"{k_new}_upper_limit"] = attach_unit(r[f"{k_original}"])
                    s[f"{k_new}_upper_limit"][is_upper == False] = np.nan

                # populate bounded constraints
                is_bounded = limit_flag == 0
                if np.sum(is_bounded) > 0:
                    s[k_new] = attach_unit(r[k_original], unit)
                    s[f"{k_new}_uncertainty_upper"] = attach_unit(
                        r[f"{k_original}err1"], unit
                    )
                    s[f"{k_new}_uncertainty_lower"] = attach_unit(
                        r[f"{k_original}err2"], unit
                    )
                    s[k_new][is_bounded == False] = np.nan
                    s[f"{k_new}_uncertainty_upper"][is_bounded == False] = np.nan
                    s[f"{k_new}_uncertainty_lower"][is_bounded == False] = np.nan

                    print(f"👇 populated {k_original} > {k_new} and errors and limits ")

            else:
                # easy, just record the column itself with no error
                s[k_new] = attach_unit(r[k_original], unit)
                print(f"🙋 populated {k_original} > {k_new} , but not 100% sure...")

            # keep track of reference for measurements
            if (k_original in self._raw_columns_with_errors_and_limits) or (
                k_original in self._raw_columns_with_errors
            ):
                try:
                    s[f"{k_new}_reference"] = self._ingest_references(r, k_original)
                    print(
                        f"⚠️ ingested reference information for {k_original} > {k_new}"
                    )
                    # s.add_index(f"{k_new}_reference")
                except (KeyError, AssertionError):
                    print(
                        f"⚠️ no reference information found for {k_original} > {k_new}"
                    )

        # basic reference information
        populate_one_or_more_columns("name", "pl_name")
        populate_one_or_more_columns("hostname", "hostname")
        populate_one_or_more_columns("letter", "pl_letter")
        populate_one_or_more_columns("gaia_id", "gaia_id")
        populate_one_or_more_columns("number_of_stars", "sy_snum")
        populate_one_or_more_columns("number_of_planets", "sy_pnum")

        # what's the history?
        populate_one_or_more_columns("discovery_method", "discoverymethod")
        populate_one_or_more_columns("discovery_year", "disc_year")
        populate_one_or_more_columns("discovery_publication", "disc_refname")
        s["discovery_publication"] = [
            parse_reflink(x) for x in s["discovery_publication"]
        ]
        populate_one_or_more_columns("discovery_facility", "disc_facility")

        # what are the host positions and kinematics?
        populate_one_or_more_columns("ra", "ra", u.deg)
        populate_one_or_more_columns("dec", "dec", u.deg)
        populate_one_or_more_columns("pmra", "sy_pmra", u.mas / u.year)
        populate_one_or_more_columns("pmdec", "sy_pmdec", u.mas / u.year)
        populate_one_or_more_columns("systemic_rv", "st_radv", u.km / u.s)
        populate_one_or_more_columns("distance", "sy_dist", u.pc)

        # what are some useful flags?
        populate_one_or_more_columns("detected_in_rv", "rv_flag")
        populate_one_or_more_columns("detected_in_pulsar", "pul_flag")
        populate_one_or_more_columns("detected_in_pulsation_timing", "ptv_flag")
        populate_one_or_more_columns("detected_in_transit", "tran_flag")
        populate_one_or_more_columns("detected_in_astrometry", "ast_flag")
        populate_one_or_more_columns(
            "detected_in_orbital_brightness_modulations", "obm_flag"
        )
        populate_one_or_more_columns("detected_in_microlensing", "micro_flag")
        populate_one_or_more_columns(
            "detected_in_eclipse_timing_variations", "etv_flag"
        )
        populate_one_or_more_columns("detected_in_imaging", "ima_flag")
        populate_one_or_more_columns("detected_in_disk_kinematics", "dkin_flag")

        populate_one_or_more_columns("is_controversial", "pl_controv_flag")
        populate_one_or_more_columns("shows_ttv", "ttv_flag")

        # (these might need merging from the composite table?!)
        populate_one_or_more_columns("stellar_spectral_type", "st_spectype")
        populate_one_or_more_columns("stellar_teff", "st_teff", u.K)
        populate_one_or_more_columns("stellar_radius", "st_rad", u.Rsun)
        populate_one_or_more_columns("stellar_mass", "st_mass", u.Msun)
        populate_one_or_more_columns("stellar_age", "st_age", u.Gyr)
        populate_one_or_more_columns("stellar_metallicity", "st_met")
        populate_one_or_more_columns("stellar_luminosity", "st_lum")
        populate_one_or_more_columns("stellar_logg", "st_logg")
        populate_one_or_more_columns("stellar_density", "st_dens")
        populate_one_or_more_columns("stellar_vsini", "st_vsin")
        populate_one_or_more_columns("stellar_rotation_period", "st_rotp")

        # what are the stellar magnitudes?
        bands = [
            "u",
            "g",
            "r",
            "i",
            "z",
            "V",
            "B",
            "IC",
            "J",
            "H",
            "K",
            "W1",
            "W2",
            "W3",
            "W4",
            "gaia",
            "T",
            "kep",
        ]
        for b in bands:
            populate_one_or_more_columns(f"magnitude_{b}", f"sy_{b.lower()}mag", u.mag)

        # what are some basic orbital parameters
        populate_one_or_more_columns("period", "pl_orbper", u.day)
        populate_one_or_more_columns("semimajoraxis", "pl_orbsmax", u.AU)
        populate_one_or_more_columns("eccentricity", "pl_orbeccen")
        populate_one_or_more_columns("argument_of_periastron", "pl_orblper", u.deg)
        populate_one_or_more_columns("inclination", "pl_orbincl", u.deg)

        # what are the planet properties? (check Jupiter isn't better?!)
        populate_one_or_more_columns("radius", "pl_rade", u.Rearth)
        populate_one_or_more_columns("mass", "pl_bmasse", u.Mearth)
        populate_one_or_more_columns("density", "pl_dens", u.g / u.cm**3)
        populate_one_or_more_columns(
            "insolation",
            "pl_insol",
            (u.Lsun / 4 / np.pi / (1 * u.AU) ** 2).to("W/m**2"),
        )
        populate_one_or_more_columns("teq", "pl_eqt", u.K)

        # does it have transmission + emission spec?
        populate_one_or_more_columns(
            "number_of_transmission_measurements", "pl_ntranspec"
        )
        populate_one_or_more_columns("number_of_emission_measurements", "pl_nespec")

        # what are the (often) transit-derived properties?
        populate_one_or_more_columns("transit_midpoint", "pl_tranmid", u.day)
        populate_one_or_more_columns("transit_duration", "pl_trandur", u.hour)
        populate_one_or_more_columns("transit_depth", "pl_trandep", 0.01)
        populate_one_or_more_columns("transit_scaled_semimajoraxis", "pl_ratdor")
        populate_one_or_more_columns("transit_impact_parameter", "pl_imppar")
        populate_one_or_more_columns("rv_semiamplitude", "pl_rvamp", u.m / u.s)
        populate_one_or_more_columns("projected_obliquity", "pl_projobliq", u.deg)
        populate_one_or_more_columns("obliquity", "pl_trueobliq", u.deg)

        # is this the default parameter set?
        populate_one_or_more_columns("default_parameter_set", "default_flag")

        # sort these planets by their names
        s.sort("name")

        # set the fill_value for any numeric columns to nans
        for k in s.colnames:
            if not isinstance(s[k][0], str):
                s[k].fill_value = np.nan

        # fill in all the masked elements to make an unmasked array with nans
        standard = s.filled()

        #
        trimmed = self._trim_bad_data(standard)

        # return that standardized table
        return trimmed

    def _trim_bad_data(self, standard):
        """
        Mask bad data from the '.standard' table.

        Planetary Systems Composite Parameters contains lots of quantities
        that aren't quite reliable empircal measurements that we want to use.
        This trims theoretical calculations and/or suspicious measurements
        out of the population, masking them and replacing their central values
        with np.nan.

        Parameters
        ----------
        standard : astropy.table.QTable
            The filled table that will populate the `.standard` table.

        Returns
        -------
        trimmed : astropy.table.QTable
            The .standard table with bad data tidied away.
        """

        print("✂️ trimming weird data")

        # make a copy of the table
        trimmed = copy.deepcopy(standard)

        # find the columns that should have uncertainties
        columns_with_uncertainties = np.unique(
            [
                x.split("_uncertainty")[0]
                for x in standard.colnames
                if "uncertainty" in x
            ]
        )
        for k in columns_with_uncertainties:

            # remove quantities with no or nonsense uncertainties
            ok = np.isfinite(standard[k])
            N_has_value = np.sum(ok)
            for w in ["lower", "upper"]:
                ok *= np.isfinite(standard[f"{k}_uncertainty_{w}"])
                # this comment is a potentially harmless or potentially troublesome KLUDGE!
                # see https://github.com/zkbt/exoatlas/issues/58
                # ok *= standard[f"{k}_uncertainty_{w}"] != 0
                print(
                    f"🚨😳🔔‼️ some values with zero-uncertainty might have snuck through for {k} 🚨😳🔔‼️🚨😳🔔‼️"
                )
            N_has_uncertainty = np.sum(ok)
            bad = ok == False
            print(
                f"""{k:>20}: {N_has_value:>4} values and {N_has_uncertainty:>4} with uncertainties ({N_has_uncertainty/N_has_value:>6.1%})"""
            )
            trimmed[k][bad] = np.nan

        print("values without reasonable uncertainties have been trimmed")
        return trimmed

    def _ingest_references(self, r, k):
        """
        Try to get the reference for a particular measurement.

        For the Planetary Systems Composite Parameters table,
        where references have merged, each measured quantity has
        its own reference encoded as a '*reflink' column name.

        Parameters
        ----------
        r : astropy.table.QTable
            The raw table being ingested.
        k : str
            The key for which a reference should be checked.
        """

        # define a column name for this quantity
        k_reference = f"{k}_reflink"

        # get a bibcode for this quantity from each row
        list_of_references = [parse_reflink(x) for x in r[k_reference]]

        return list_of_references


class ExoplanetsPS(ExoplanetsPSCP):
    """
    Exoplanets from the NASA Exoplanet Archive.

    This uses the Planetary Systems table,
    which includes one row per planet per reference,
    so it's a massive dataset providing multiple options
    for most planets, although any one particular reference
    might be missing the complete set of parameters needed
    to fully characterize any one system.
    """

    label = "ExoplanetsPS"
    _downloader = ExoplanetArchiveDownloader("ps")

    def __init__(self, *args, **kwargs):
        """
        Initialize a population of all known exoplanets
        using the NASA Exoplanet Archive `ps` table.

        This contains many more rows than there are planets,
        because there's one row for each reference (for each planet).

        This also sets all the `*_reference` columns to serve as
        table indices, to make them easy to extract for updating
        individual planet references.
        """
        ExoplanetsPSCP.__init__(self, *args, **kwargs)
        self._add_references_as_indices()

    def _ingest_references(self, r, k):
        """
        Try to get the reference for a particular measurement.

        For the Planetary Systems table, where references are
        necessarily consistent across parameter sets, there
        are only three reference columns:
            'pl_refname' for the planet and orbit
            'st_refname' for the stellar properties
            'sy_refname' for the system positions + motions

        Parameters
        ----------
        r : astropy.table.QTable
            The raw table being ingested.
        k : str
            The key for which a reference should be checked.
        """

        # is this key 'pl', 'st', or 'sy'?
        prefix = k.split("_")[0]
        assert prefix in ["pl", "st", "sy"]

        # define a column name for this quantity
        k_reference = f"{prefix}_refname"

        # get a bibcode for this quantity from each row
        list_of_references = [parse_reflink(x) for x in r[k_reference]]

        return list_of_references


class Exoplanets(ExoplanetsPSCP):
    """
    Exoplanets from the NASA Exoplanet Archive.

    This uses the table of Planetary Systems Composite Parameters
    (one row per planet, references merged) but it secretly also
    loads the table of Planetary Systems (one row per reference,
    with multiple rows per planet), so that we can choose to use
    different references for any particular planet.
    """

    label = "Exoplanets"

    def __init__(self, remake=False, **plotkw):
        """
        Initialize a population of all known exoplanet from
        the NASA Exoplanet Archive `pscomppars` and `ps` tables.

        The standard population comes only from `pscomppars`.
        However, it should secretly also download the `ps` table,
        which would enable the user to choose to replace
        parameters for a particular planet by choosing a
        particular reference.

        Parameters
        ----------
        remake : bool
            Should the population be remade from raw files,
            even if a recent version already exists?
        **plotkw : dict
            All other keywords will go toward defining
            default plotting styles, for example like
            'alpha', 'marker', 'zorder', 'color', ...
        """

        # load standard table(s) or ingest from raw data
        PredefinedPopulation.__init__(self, remake=remake, **plotkw)

    def __getitem__(self, key):
        """
        Create a subpopulation of planets by indexing, slicing, or masking.

        `Exoplanets` needs its own version of this (beyond the one inherited from
        the overarching `Population` class) because it secretly contains a hidden
        `ExoplanetsPS` table that should also be trimmed down to only relevant
        rows every time that we select a subset.

        Parameters
        ----------
        key : int, list, array, slice
            An index, slice, or mask to select a subset of the population.

        Returns
        -------
        subset : Population
            A subset of the population, as set by the index, slice, or mask.
        """

        # first call the normal `Population` indexing/slicing/masking
        subset = super().__getitem__(key)

        # then trim the .ps population to only those (host)names in the subset
        if hasattr(self, "individual_references"):
            subset.individual_references = self.individual_references[
                list(np.unique(subset.tidyname()))
            ]
            subset.individual_references.label = "Individual References"

        return subset

    def load_individual_references(self):
        """
        Populate an internal `.individual_references` population
        with individual references for planets.

        The `pscomppars` table from which the Exoplanets population
        is defined merged together multiple different references.
        Sometimes we might want be able to see and/or choose the
        values associated with one particular reference. This
        function loads a separate Exoplanet-like population
        containing those individual references and stores it
        in the `.individual_references` attribute.
        """

        # skip slowly loading the data if it already exists
        try:
            return self.individual_references
        except AttributeError:

            print("(loading individual references may take up to a few minutes)")

            # load the internal population with data for individual references
            self.individual_references = ExoplanetsPS()
            # then trim the .ps population to only those (host)names in the subset
            self.individual_references = self.individual_references[
                list(np.unique(self.tidyname()))
            ]

            self.individual_references.label = "Individual References"
            return self.individual_references

    def plot_measurements_from_individual_references(self, x="mass", y="radius"):
        """
        Plot all available measurements from individual references.

        Parameters
        ----------
        x : str
            The quantity to plot on the x-axis.
        y : str
            The quantity to plot on the y-axis.
        """
        for p, c in zip([self.individual_references, self], ["orchid", "black"]):
            plt.errorbar(
                p.get(x),
                p.get(y),
                p.get_uncertainty_lowerupper(y),
                p.get_uncertainty_lowerupper(x),
                linewidth=0,
                elinewidth=1,
                marker="o",
                color=c,
                label=p.label,
            )
        plt.xlabel(x)
        plt.ylabel(y)
        plt.legend()

    def display_individual_references(
        self,
        keys=["mass", "radius", "stellar_mass", "stellar_radius", "stellar_teff"],
        planets=None,
    ):
        """
        Summarize the available measurements for a particular quantity.

        Parameters
        ----------
        planets : str, list
            One planet to investigate, as a string,
            or multiple planets to investigate, as a list of strings,
            or, if None, show the entire population.
        keys : str, list
            One quantity to investigate, as a string,
            or multiple quantities as a list of strings.

        Returns
        -------
        table : Table
            An abbreviated table summarizing the references and the
            measurements they each provide, for the requested key(s).
            The table will be sorted by planet name, to facilitate
            comparing values within a particular system.
        """

        # extract the rows relevant to the requested planets
        if planets == None:
            pop = self
        else:
            pop = self[planets]

        # replace a single key with a list of keys
        if isinstance(keys, str):
            keys = [keys]

        # construct a list of columsn to appear in the final table
        columns_to_include = ["name"]

        # add all columns relevant to the requested key(s)
        available_colnames = pop.standard.colnames
        for k in keys:
            if f"{k}_reference" in available_colnames:
                for suffix in suffixes:
                    if f"{k}{suffix}" in available_colnames:
                        columns_to_include += [f"{k}{suffix}"]

        # construct a subset table with just those columns, first row as actual, others as all available options
        table = vstack(
            [
                pop.standard[columns_to_include],
                pop.individual_references.standard[columns_to_include],
            ]
        )
        table.add_column(["🏷️"] * len(table), index=0, name="is_being_used")
        table["is_being_used"][: len(pop)] = "🗺️"
        table.sort("name")

        return table

    def update_reference(self, references, keys=None, planets=None, verbose=False):
        """
        Modify the internal standardized data table to use
        measurements from a particular reference.

        Each planet may have multiple measurements available in the
        NASA Exoplanet Archive. If you want to switch from the default
        parameters that come from the merged Planetary Systems
        Composite Parameters table to those coming from a particular
        individual reference in the Planetary Systems table, you might...

        1.  Run `.load_individual_references()` to make sure that the
            data for individual references are loaded and stored in the
            `.individul_references` attribute.
        2.  Run `.display_individual_references()` to see what references
            have useful data you might want to use for a particular planet.
        3.  Run `.update_reference()` to use a reference to update the
            parameters in the standardized population table.

        Running this method will permanently modify the `.standard` table
        inside the `Exoplanets` population. The only way to be sure to
        undo the effect is to create a new `Exoplanets` population.

        This will only update finite values; it won't overwrite existing
        values with nan. To do so, we either need to add a keyword option
        here, or you'll need to use `.update_values`.

        TO-DO: This feels slow and inefficient; I'm sure there's a faster
        and/or classier way to do this, but let's call this good for now!

        Parameters
        ----------
        references : str, list
            One reference to try to use, as a string,
            or multiple references to try to use, as a list of strings.
            References will be updated in order, so later references
            in a list will overwrite earlier ones.
        keys : str, list
            One quantity to try to update, as a string,
            or multiple quantities as a list of strings.
        planets : str, list
            One planet to try to update, as a string,
            or multiple planets to update, as a list of strings,
            or, if None, consider the entire population.
        verbose : bool
            Whether all changes should be displayed.
        """

        warning_message = f"""
        Your request to update references with...
            references = {references}
            keys = {keys}
            planets = {planets}
        ...failed somehow. Sorry!
        """
        # make sure references are a list, even with just one element
        if isinstance(references, str):
            references = [references]

        # make sure keys are a list, which might be all columns with references
        if keys == None:
            keys_with_reference = [
                x.split("_reference")[0]
                for x in self.standard.colnames
                if "_reference" in x
            ]
            keys = [k for k in keys_with_reference if k in self.standard.colnames]
        elif isinstance(keys, str):
            keys = [keys]

        # make sure planets are a list, which might be all planets in population
        if planets == None:
            planets_to_index = slice(None, None, None)
        elif isinstance(planets, str):
            planets_to_index = [clean(planets).lower()]
        else:
            planets_to_index = [clean(x).lower() for x in planets]

        # make sure something happens
        nothing_happened = True

        # loop over columns
        for k in keys:

            # extract a table of just these references, or give up
            try:
                # print(f"{k}_reference", references)
                these_references = QTable(
                    self.individual_references.standard.loc[
                        f"{k}_reference", references
                    ]
                )
            except KeyError:
                # if verbose:
                #    print(f'No reference(s) "{references}" found for key "{k}"')
                continue

            # extract a (possibly even smaller) table of just these planets, or give up
            try:
                these_planets = QTable(
                    these_references.loc["tidyname", planets_to_index]
                )
            except KeyError:
                if verbose:
                    print(
                        'No planet(s) "{planets_to_index}" found for reference "{references}"; moving on.'
                    )
                continue

            # find where there is a real new measurement for this key
            is_finite = np.isfinite(these_planets[f"{k}"])
            # if verbose:
            #    print(
            #        k,
            #        planets_to_index,
            #        is_finite,
            #        type(is_finite).__name__,
            #        (type(these_planets).__name__),
            #    )
            new = copy.deepcopy(these_planets[is_finite])
            new.add_index("tidyname")

            # if verbose:
            #    print("The new properties to be adopted:")
            #    display(new)

            # construct a new table of the references we want to update for this key
            if len(new) > 0:
                for tidyname in np.unique(new["tidyname"]):
                    if verbose:
                        print(f"'{k}' for '{tidyname}'")
                    old_for_this_planet = self.standard.loc["tidyname", tidyname]
                    new_for_this_planet = new.loc["tidyname", tidyname]
                    # if verbose:
                    # print("old", old_for_this_planet["tidyname"])
                    # print("new", new_for_this_planet["tidyname"])
                    assert (
                        old_for_this_planet["tidyname"]
                        == new_for_this_planet["tidyname"]
                    )

                    keys_to_display = ["tidyname"]
                    if verbose:
                        old_for_display = copy.deepcopy(old_for_this_planet)

                    for s in suffixes:
                        try:
                            if verbose:
                                old_value = old_for_this_planet[f"{k}{s}"]
                                new_value = new_for_this_planet[f"{k}{s}"]
                                # print(f"{k}{s}: {old_value} > {new_value}")
                            old_for_this_planet[f"{k}{s}"] = new_for_this_planet[
                                f"{k}{s}"
                            ]
                            keys_to_display.append(f"{k}{s}")
                        except (IndexError, KeyError):
                            pass
                    if verbose:
                        change_summary = vstack(
                            [
                                old_for_display[keys_to_display],
                                new_for_this_planet[keys_to_display],
                            ]
                        )
                        change_summary.add_column(
                            ["old", "new"], index=0, name="version"
                        )
                        display(change_summary)

                nothing_happened = False

        if nothing_happened:
            warnings.warn(warning_message)
