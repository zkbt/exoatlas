from .imports import *
from .population import Population


initial_filename = directories["data"] + "TESSsimulations.tsv"


class PredictedTESS(PredefinedPopulation):
    """TESS population object contains a simulated TESS planet yield
    from 200,000 two-minute cadence postage stamps,
    as calculated by Peter Sullivan et al. (2015)"""

    def __init__(self, **kwargs):
        """Initialize a population of simulated TESS planets."""
        Population.__init__(self, label="Predicted TESS", **kwargs)
        self.color = "darkorange"
        self.zorder = 1

    def loadFromScratch(self):
        self.table = ascii.read(initial_filename)
        self.speak("loaded TESS simulated population from {0}".format(initial_filename))

    def trim_raw(self):
        self.trimmed = self.table

    def create_standardardized(self):
        t = self.trimmed
        s = Table()
        s["name"] = ["tess{0:04}i".format(i) for i in range(len(t))]
        s["kepid"] = None
        s["period"] = t["P"]
        s["stellar_teff"] = t["stellar_teff"]
        s["stellar_radius"] = t["Rstar"]
        s["Jmag"] = t["Jmag"]
        s["Vmag"] = t["Vmag"]

        s["radius"] = t["Rplanet"]
        s["radius_uncertainty_lower"] = np.zeros_like(t["Rplanet"]) + 0.0001
        s["radius_uncertainty_upper"] = np.zeros_like(t["Rplanet"]) + 0.0001

        s["mass"] = np.nan + np.zeros_like(t["Rplanet"])
        s["mass_uncertainty_lower"] = np.nan + np.zeros_like(t["Rplanet"])
        s["mass_uncertainty_upper"] = np.nan + np.zeros_like(t["Rplanet"])

        s["teq"] = 280.0 * t["S/SEarth"] ** 0.25
        s["transit_ar"] = 0.5 * (s["stellar_teff"] / s["teq"]) ** 2
        s["rv_semiamplitude"] = t["K"]
        s["radius_ratio"] = (
            t["Rplanet"] * craftroom.units.Rearth / (t["Rstar"] * craftroom.units.Rsun)
        )
        s["distance"] = 10 * 10 ** (0.2 * t["Dist"]) + np.nan
        s["ra"] = t["RA"]
        s["dec"] = t["Dec"]

        self.standard = s
