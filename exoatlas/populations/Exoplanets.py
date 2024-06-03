# exoplanet population of all "confirmed" exoplanets from exoplanet archive
from ..imports import *
from .Population import PredefinedPopulation
from .downloaders import *

__all__ = ["Exoplanets", "ExoplanetsComposite"]


# apply some kludges to correct bad planet properties
class Exoplanets(PredefinedPopulation):
    def __init__(self, label="Exoplanets", remake=False, **plotkw):
        """
        Initialize a population of all known exoplanets
        using the NASA Exoplanet Archive `ps` table.

        This generates an Exoplanets population containing
        planets discovered through any method. It tries to
        be clever about loading cached files so that it
        will work quickly and not try to redownload everything
        from scratch every time it's called.

        Parameters
        ----------
        label : string
            A default label to be associated with this population.
            This will show up in the legend on population
            comparison plots.
        remake : bool
            Should the population be remade from raw files,
            even if a recent version already exists?
        **plotkw : dict
            All other keywords will go toward defining
            default plotting styles, for example like
            'alpha', 'marker', 'zorder', 'color', ...
        """
        # set up the population
        self._downloader = planetary_systems_downloader
        PredefinedPopulation.__init__(self, label=label, remake=remake, **plotkw)

    def load_raw(self, remake=False):
        """
        Load the raw table of data.

        Parameters
        ----------
        remake : bool
            Should the raw table be re-made/re-downloaded,
            even if a recent one exists?

        Returns
        -------
        raw : astropy.table.Table
            A raw, untrimmed, unstandardized table.
        """

        # load (or download) the table of exoplanet properties
        raw = self._downloader.get(remake=remake)

        # for debugging, hang on to the raw table as a hidden attribute
        self._raw = raw

        return raw

    def trim_raw(self, raw):
        """
        Trim the raw table down to ones with reasonable values.

        Parameters
        ----------
        raw : astropy.table.Table
            A raw, untrimmed, unstandardized table.

        Returns
        -------
        trimmed : astropy.table.Table
            A raw, trimmed, unstandardized table.
        """

        masks = {}

        ok = np.ones(len(raw)).astype(bool)
        for k in masks:
            ok *= masks[k]
            N = sum(ok == True)
            self.speak(f"{N} planets pass the `{k}` filter")

        # trim down the table to just those that are OK
        trimmed = raw[ok]

        # for debugging, hang onto the trimmed table as a hidden attribute
        self._trimmed = trimmed

        # report how many points got trimmed away
        ntotal = len(raw)
        ntrimmed = len(trimmed)
        self.speak(f"trimmed down to {ntrimmed}/{ntotal} rows")

        return trimmed

    def create_standard(self, trimmed):
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
        trimmed : astropy.table.Table
            A raw, trimmed, unstandardized table.

        Returns
        -------
        standard : astropy.table.Table
            A standardize table of exoplanet properties.
        """

        # define the table from which we're deriving everying
        r = trimmed

        # the new standardized table
        s = Table()

        def strip_even_if_masked(x):
            try:
                return x.strip()
            except AttributeError:
                return None

        def attach_unit(x, unit=None):
            if isinstance(x[0], str):
                return [strip_even_if_masked(a) for a in x]
            if unit is None:
                return x
            else:
                return x * unit

        def populate_one_or_more_columns(k_new, k_original, unit=None):
            """
            A wrapper to help smoothly grab quantities with limits.

            Parameters
            ----------
            k_original : string
                The column name in the original, raw table.
            k_new : string
                The column name in the new, standardized table.
            unit : None, astropy.units.Unit
                A unit to attach to the quantity, if necessary.
            """

            if k_original in r.meta["columns-without-errors"]:
                # easy, just record the column itself with no error
                s[k_new] = attach_unit(r[k_original], unit)
                print(f"populated {k_new} with {k_original}")

            elif k_original in r.meta["columns-with-errors"]:
                # record the column itself
                s[k_new] = attach_unit(r[k_original], unit)

                # record the upper and lower errorbars
                s[f"{k_new}_uncertainty_upper"] = attach_unit(
                    r[f"{k_original}err1"], unit
                )
                s[f"{k_new}_uncertainty_lower"] = attach_unit(
                    r[f"{k_original}err2"], unit
                )
                print(f"populated {k_new} and errors with {k_original}")

            elif k_original in r.meta["columns-with-errors-and-limits"]:
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

                    print(f"populated {k_new} and errors and limits with {k_original}")

            else:
                # easy, just record the column itself with no error
                s[k_new] = attach_unit(r[k_original], unit)
                print(f"populated {k_new} with {k_original}, but not 100% sure...")

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
        populate_one_or_more_columns("discovery_reference", "disc_refname")
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
            populate_one_or_more_columns(f"{b}mag", f"sy_{b.lower()}mag", u.mag)

        # what are some basic orbital parameters
        populate_one_or_more_columns("period", "pl_orbper", u.day)
        populate_one_or_more_columns("semimajoraxis", "pl_orbsmax", u.AU)
        populate_one_or_more_columns("eccentricity", "pl_orbeccen")
        populate_one_or_more_columns("omega", "pl_orblper", u.deg)
        populate_one_or_more_columns("inclination", "pl_orbincl", u.deg)

        # what are the planet properties? (check Jupiter isn't better?!)
        populate_one_or_more_columns("radius", "pl_rade", u.Rearth)

        # kludge for dealing with different masses from different tables?
        try:
            # ps?
            populate_one_or_more_columns("mass", "pl_masse", u.Mearth)
            populate_one_or_more_columns("msini", "pl_msinie", u.Mearth)
        except KeyError:
            # pscomp?
            populate_one_or_more_columns("mass", "pl_bmasse", u.Mearth)
            x = "mass"
            provenance = r["pl_bmassprov"]
            is_measured_mass = (provenance == "Mass") | (provenance == "Msin(i)/sin(i)")
            for k in [
                x,
                f"{x}_uncertainty_upper",
                f"{x}_uncertainty_lower",
                f"{x}_lower_limit",
                f"{x}_upper_limit",
            ]:
                s[k][is_measured_mass == False].mask = True

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
        populate_one_or_more_columns("transit_ar", "pl_ratdor")
        populate_one_or_more_columns("transit_b", "pl_imppar")
        populate_one_or_more_columns("transit_depth", "pl_trandep", 0.01)
        populate_one_or_more_columns("transit_duration", "pl_trandur", u.hour)
        populate_one_or_more_columns("rv_semiamplitude", "pl_rvamp", u.m / u.s)
        populate_one_or_more_columns("projected_obliquity", "pl_projobliq", u.deg)
        populate_one_or_more_columns("obliquity", "pl_trueobliq", u.deg)

        populate_one_or_more_columns("default_parameter_set", "default_flag")

        # trim to default parameter set (be more careful here!!!)
        s = s[s["default_parameter_set"] == 1]

        # sort these planets by their names
        s.sort("name")

        # set the fill_value for any numeric columns to nans
        for k in s.colnames:
            if not isinstance(s[k][0], str):
                s[k].fill_value = np.nan

        # fill in all the masked elements to make an unmasked array with nans
        standard = s.filled()

        # return that standardized table
        return standard


class ExoplanetsComposite(Exoplanets):
    def __init__(self, label="ExoplanetsComposite", remake=False, **plotkw):
        """
        Initialize a population of all known exoplanets
        using the NASA Exoplanet Archive `ps` table.

        This generates an Exoplanets population containing
        planets discovered through any method. It tries to
        be clever about loading cached files so that it
        will work quickly and not try to redownload everything
        from scratch every time it's called.

        Parameters
        ----------
        label : string
            A default label to be associated with this population.
            This will show up in the legend on population
            comparison plots.
        remake : bool
            Should the population be remade from raw files,
            even if a recent version already exists?
        **plotkw : dict
            All other keywords will go toward defining
            default plotting styles, for example like
            'alpha', 'marker', 'zorder', 'color', ...
        """
        # set up the population
        self._downloader = composite_planetary_systems_downloader
        PredefinedPopulation.__init__(self, label=label, remake=remake, **plotkw)
