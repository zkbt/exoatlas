# exoplanet population of all "confirmed" exoplanets from exoplanet archive
from ..imports import *
from .Population import PredefinedPopulation
from .downloaders import merged_exoplanets

__all__ = ["Exoplanets"]

# apply some kludges to correct bad planet properties
class Exoplanets(PredefinedPopulation):
    def __init__(self, label="All Exoplanets", remake=False, **plotkw):
        """
        Initialize a population of all known exoplanets,
        from a table downloaded from the NASA Exoplanet Archive.
        """
        # set up the population
        PredefinedPopulation.__init__(self, label=label, remake=remake, **plotkw)

    def load_raw(self, remake=False):
        """
        Load the raw table of data from the NASA Exoplanet Archive.
        """

        # load (or download) the table of composite exoplanet properties
        raw = merged_exoplanets.get(remake=remake)

        # for debugging, hang on to the raw table as a hidden attribute
        self._raw = raw

        return raw

    def trim_raw(self, raw):
        """
        Trim the raw table down to ones with reasonable values.
        """

        masks = {}

        ok = np.ones(len(raw)).astype(np.bool)
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
        Create a standardized table, pulling at least the necessary columns
        from the raw table and potentially including others too.
        """

        # define the table from which we're deriving everying
        t = trimmed

        s = Table()

        # what's the name of the planet?
        s["name"] = [x.strip() for x in t["pl_name"]]
        s["hostname"] = [x.strip() for x in t["hostname"]]

        # badpos = (t['ra'] ==0.0)*(t['dec'] == 0.0)
        s["ra"] = t["ra"] * u.deg  # t.MaskedColumn(t['ra'], mask=badpos)
        s["dec"] = t["dec"] * u.deg  # t.MaskedColumn(t['dec'], mask=badpos)

        def merge(k):
            """
            For some particular column from the raw table, take all the data
            we can from the confirmed planets table, and then
            try to fill in from the composite planets table.
            """

            # make an array from the confirmed table
            x = t[k].copy()

            # try to fill in the bad ones from the composite table
            bad = x.mask
            n = len(x)
            nmissing = sum(bad)
            self.speak(f"{nmissing}/{n} missing from confirmed")

            x[bad] = t["f" + k][bad]
            self.speak(f"{nmissing}/{n} missing from confirmed and composite")

            return x

        # what's the period of the planet?
        s["period"] = t["pl_orbper"] * u.day
        s["semimajoraxis"] = t["pl_orbsmax"] * u.AU
        s["e"] = t["pl_orbeccen"]
        s["omega"] = t["pl_orblper"] * u.deg
        s["inclination"] = t["pl_orbincl"] * u.deg

        # what are the observed transit properties?
        s["transit_epoch"] = t["pl_tranmid"] * u.day
        s["transit_duration"] = t["pl_trandur"] * u.day
        s["transit_depth"] = t["pl_trandep"] / 100.0

        # what are the basic stellar properties?
        s["stellar_teff"] = merge("st_teff") * u.K
        s["stellar_radius"] = merge("st_rad") * u.Rsun
        s["stellar_mass"] = merge("st_mass") * u.Msun

        # add stellar age
        s["stellar_age"] = merge("st_age") * u.Gyr
        upper = t["st_ageerr1"].data
        lower = t["st_ageerr2"].data
        bad = upper.mask | lower.mask
        upper[bad] = np.inf
        lower[bad] = np.inf
        s["stellar_age_uncertainty_upper"] = upper * u.Gyr
        s["stellar_age_uncertainty_lower"] = lower * u.Gyr

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
            s[f"{b}mag"] = t[f"sy_{b.lower()}mag"]

        # what is the planet radius?
        s["radius"] = t["pl_rade"] * u.Rearth

        # pull out the radius uncertainties
        upper = t["pl_radeerr1"].data
        lower = t["pl_radeerr2"].data
        bad = upper.mask | lower.mask
        upper[bad] = np.inf
        lower[bad] = np.inf
        s["radius_uncertainty_upper"] = upper * u.Rearth
        s["radius_uncertainty_lower"] = lower * u.Rearth

        # what are the (often) transit-derived properties?
        s["transit_ar"] = t["pl_ratdor"]
        s["transit_b"] = t["pl_imppar"]

        # KLUDGE?
        # s['rv_semiamplitude'] =  t['pl_rvamp'] #t.MaskedColumn(t['K'], mask=t['K']==0.0)

        s["mass"] = t["pl_masse"] * u.Mearth

        # pull out the mass uncertainties
        upper = t["pl_masseerr1"].data
        lower = t["pl_masseerr2"].data
        bad = upper.mask | lower.mask
        upper[bad] = np.inf
        lower[bad] = np.inf
        s["mass_uncertainty_upper"] = upper * u.Mearth
        s["mass_uncertainty_lower"] = lower * u.Mearth

        # keep track of msini (for non-transiting planets)
        s["msini"] = t["pl_msinie"] * u.Mearth

        # how far away is the star?
        s["distance"] = merge("sy_dist") * u.pc
        s["distance_uncertainty_upper"] = merge("sy_disterr1") * u.pc
        s["distance_uncertainty_lower"] = merge("sy_disterr2") * u.pc

        # what facility discovered the planet?
        s["discoverer"] = t["disc_facility"]
        s["method"] = t["discoverymethod"]

        s["has_transit"] = t["tran_flag"]
        s["has_RV"] = t["rv_flag"]

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
