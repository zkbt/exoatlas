# exoplanet population of all "confirmed" exoplanets from exoplanet archive
from ..imports import *
from .Population import PredefinedPopulation
from .downloaders import toi_merged

__all__ = ["TOI"]

# url = 'https://exofop.ipac.caltech.edu/tess/download_ctoi.php?sort=ctoi&output=pipe'
# initial_filename = directories['data'] + 'TOI-from-exofop.psv'
#
# def downloadLatest():
#    print('downloading the latest list of TOI candidates from ExoFOP')
#    request.urlretrieve(url, initial_filename)


class TOI(PredefinedPopulation):
    def __init__(self, label="TOI", remake=False, **kw):
        """
        Initialize a population of TESS Objects of Interest
        from a table downloaded from the NASA Exoplanet Archive.
        """
        # set up the population
        PredefinedPopulation.__init__(self, label=label, remake=remake, **kw)
        self.color = "crimson"

    def load_raw(self, remake=False):
        """
        Load the raw table of TOI data from the NASA Exoplanet Archive.
        """

        # load (or download) the table of composite exoplanet properties
        raw = toi_merged.get(remake=remake)

        # for debugging, hang on to the raw table as a hidden attribute
        self._raw = raw

        return raw

    def trim_raw(self, raw):
        """
        Trim the raw table down to ones with reasonable values.
        """

        masks = {}

        with np.errstate(invalid="ignore"):

            # is this a relatively cool star?
            masks["cool"] = raw["Stellar Eff Temp (K)"] < 20000

            # is this not a single-transit candidate (without period)?
            masks["notsingle"] = raw["Period (days)"] > 0

        ok = np.ones(len(raw)).astype(np.bool)
        for k in masks:
            ok *= masks[k]
            N = sum(ok == True)
            self.speak(f"{N} planets pass the `{k}` filter")

        # trim down the table to just those that are OK
        trimmed = raw[ok]

        # for debugging, hang onto the trimmed table as a hidden attribute
        self._trimmed = trimmed

        ntotal = len(raw)
        ntrimmed = len(trimmed)
        self.speak(f"trimmed down to {ntrimmed}/{ntotal} rows")

        return trimmed

    def create_standard(self, trimmed):
        """
        Create a standardized table, pulling at least the necessary columns
        from the raw table and potentially including others too.
        """

        t = trimmed
        n = len(t)
        s = Table()

        # use the TOI as the name
        s["name"] = [f"TOI{toi:.2f}" for toi in t["TOI"]]
        s["hostname"] = [f"TOI{toi:.0f}" for toi in t["TOI"]]

        s["TIC ID"] = t["TIC ID"]  # maybe remove, carefully?
        s["tic"] = s["TIC ID"]

        s["distance"] = t["Stellar Distance (pc)"] * u.pc
        s["distance_uncertainty"] = t["Stellar Distance (pc) err"] * u.pc
        s["discoverer"] = "TESS"

        s["period"] = t["Period (days)"] * u.day
        s["period_uncertainty"] = t["Period (days) err"] * u.day

        s["e"] = np.nan
        s["omega"] = np.nan * u.deg

        s["transit_epoch"] = t["Epoch (BJD)"] * u.day
        # s['transit_epoch_uncertainty'] = t['Epoch (BJD) err']*u.day

        s["transit_duration"] = t["Duration (hours)"] * u.hour
        s["transit_duration_uncertainty"] = t["Duration (hours) err"] * u.hour

        s["transit_depth"] = t["Depth (ppm)"] / 1e6
        s["transit_depth_uncertainty"] = t["Depth (ppm) err"] / 1e6

        s["stellar_teff"] = t["Stellar Eff Temp (K)"] * u.K
        s["stellar_teff_uncertainty"] = t["Stellar Eff Temp (K) err"] * u.K

        s["stellar_radius"] = t["Stellar Radius (R_Sun)"] * u.Rsun
        s["stellar_radius_uncertainty"] = t["Stellar Radius (R_Sun) err"] * u.Rsun

        s["radius"] = t["Planet Radius (R_Earth)"] * u.Rearth
        s["radius_uncertainty"] = t["Planet Radius (R_Earth) err"] * u.Rearth

        s["mass"] = np.nan * u.Mearth
        s["mass_uncertainty"] = np.nan * u.Mearth

        s["semimajoraxis"] = np.nan * u.AU

        # KLUDGE! FIXME!
        s["inclination"] = 90 * u.deg
        s["transit_ar"] = np.nan
        s["transit_b"] = np.nan

        # pull out some magnitudes
        s["Jmag"] = t["TIC Jmag"]
        s["Hmag"] = t["TIC Hmag"]
        s["Kmag"] = t["TIC Kmag"]
        s["Vmag"] = t["TIC Vmag"]
        s["Gmag"] = t["TIC GAIAmag"]
        s["GBPmag"] = t["TIC gaiabp"]
        s["GRPmag"] = t["TIC gaiarp"]
        s["Tmag"] = t["TIC Tmag"]
        s["Bmag"] = t["TIC Bmag"]

        s["stellar_mass"] = t["TIC mass"] * u.Msun
        s["stellar_mass_uncertainty"] = t["TIC e_mass"] * u.Msun

        s["tic_luminosity"] = t["TIC lum"] * u.Lsun
        s["tic_luminosity_uncertainty"] = t["TIC e_lum"] * u.Lsun

        s["ra"] = t["TIC ra"] * u.deg
        s["dec"] = t["TIC dec"] * u.deg

        othercolumns = {
            "ACWG": "ACWG",
            "Comments": "Comments",
            #'Date TOI Edited (UTC)': 'Date TOI Edited (UTC)',
            "Master": "Priority-Master",
            #'Planet Num': 'Planet Number',
            #'Planet SNR': 'Planet SNR',
            "SG1A": "Priority-SG1A",
            "SG1B": "Priority-SG1B",
            "SG2": "Priority-SG2",
            "SG3": "Priority-SG3",
            "SG4": "Priority-SG4",
            "SG5": "Priority-SG5",
            "Sectors": "Sectors",
            "Source": "Source",
            "TESS Disposition": "TESS Disposition",
            "TFOPWG Disposition": "TFOPWG Disposition",
            "TIC GAIA": "TIC GAIA",
            "TIC ID": "TIC ID",
            "TIC KIC": "TIC KIC",
            "TIC PARflag": "TIC PARflag",
            "TIC TWOMASS": "TIC TWOMASS",
            "TIC contratio": "TIC contratio",
            "TIC disposition": "TIC disposition",
        }
        for c in othercolumns:
            k = othercolumns[c]
            s[k] = t[c]

        s.sort("name")
        standard = s.filled()
        return standard

    @property
    def is_falsepositive(self):
        return self.standard["TFOPWG Disposition"] == "FP"

    @property
    def is_knownplanet(self):
        tess = self.standard["TESS Disposition"] == "KP"
        tfop = self.standard["TFOPWG Disposition"] == "KP"
        return tess | tfop

    @property
    def is_candidateplanet(self):
        return self.standard["TFOPWG Disposition"] == "CP"


"""class TransitingExoplanetsSubset(TOI):
    def __init__(self, label, color='black', zorder=0):

        # set the label
        self.label=label
        self.color=color
        self.zorder=zorder
        try:
            # first try to load this population
            Talker.__init__(self)
            self.load_standard()
        except IOError:
            # if that fails, recreate it from the confirmed population
            KOI.__init__(self)
            self.label=label
            self.selectSubsample()
        self.respond_to_color=True

    def selectSubsample(self):
        tr = self.toRemove()
        self.speak('removing {0} rows'.format(np.sum(tr)))
        self.removeRows(tr)
        self.speak('leaving {0} rows'.format(self.n))
        self.save_standard()

"""
"""class UnconfirmedKepler(TransitingExoplanetsSubset):
    def __init__(self):
        TransitingExoplanetsSubset.__init__(self, label="Kepler (candidates)", color='gray', zorder=-1e6)
        self.respond_to_color=True

    def toRemove(self):
        isconfirmed = self.standard['disposition'] == 'CONFIRMED'
        isjunk = self.distance == 10.0
        return isconfirmed | isjunk"""
