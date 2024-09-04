from ...imports import *
from ..predefined import PredefinedPopulation
from .sbdb_downloader import *


class SolarSystemMinorPlanets(PredefinedPopulation):
    label = "Solar System Minor Planets"

    def __init__(self, minimum_diameter=100 * u.km, **kw):
        self.minimum_diameter = minimum_diameter
        self._downloader = SBDBDownloader(minimum_diameter=minimum_diameter)
        PredefinedPopulation.__init__(self, **kw)

        self.color = "midnightblue"
        self.zorder = 1
        self.s = 20
        self.respond_to_color = False
        self.exact = True
        self.marker = "s"

    @property
    def fileprefix(self):
        """
        Define a fileprefix for this population, to be used
        for setting the filename of the standardized population.
        """
        return f'{clean(self.label)}-{self.minimum_diameter.to_value("km")}km'

    def create_standardardized(self, raw):
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
            A raw unstandardized table from `SBDBDownloader`.

        Returns
        -------
        standard : astropy.table.QTable
            A standardize table of exoplanet properties.
        """

        # define the table from which we're deriving everying
        r = raw

        # the new standardized table
        s = QTable()

        # names
        for k in ["name", "full_name", "class"]:
            s[k] = raw[k]

        s["eccentricity"] = raw["e"]
        s["semimajoraxis"] = raw["a"] * u.AU
        s["perihelion"] = raw["q"] * u.AU
        s["inclination"] = raw["i"] * u.deg
        s["Omega"] = raw["om"] * u.deg
        s["omega"] = raw["w"] * u.deg
        s["mean_anomaly"] = raw["ma"] * u.deg
        s["time_of_perihelion"] = raw["tp"] * u.day
        s["period"] = raw["per"] * u.day
        s["H"] = raw["H"]
        s["G"] = raw["G"]
        s["radius"] = (0.5 * raw["diameter"] * u.km).to(u.Rearth)
        s["mass"] = (raw["GM"] * u.km**3 / u.s**2 / con.G).to(u.Mearth)
        s["density"] = raw["density"] * u.g / u.cm**3
        s["rotational_period"] = raw["rot_per"] * u.hour
        s["albedo"] = raw["albedo"]
        s["hostname"] = "Sun"
        s["stellar_teff"] = 5780 * u.K
        s["stellar_radius"] = 1 * u.Rsun
        s["stellar_mass"] = 1 * u.Msun
        s["stellar_age"] = 4.5 * u.Gyr

        # for k in ["mass", "radius"]:
        #    for side in ["lower", "upper"]:
        #        s[f"{k}_uncertainty_{side}"] = 0 * s[k].unit

        # sort these minor planets by their names (with numbers)
        s.sort("full_name")

        # set the fill_value for any numeric columns to nans
        for k in s.colnames:
            if not isinstance(s[k][0], str):
                s[k].fill_value = np.nan

        # fill in all the masked elements to make an unmasked array with nans
        standard = s.filled()

        return standard
