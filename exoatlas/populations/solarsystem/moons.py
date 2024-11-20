from ...imports import *
from ..predefined import PredefinedPopulation
from .planets import SolarSystem

__all__ = ["SolarSystemMoons"]

initial_filename = os.path.join(
    code_directory, "populations/data/solarsystem/jpl-ssd-moons.txt"
)


class SolarSystemMoons(PredefinedPopulation):
    """
    The Solar System Moons, very crudely.
    All orbital parameters are of the planet around its star.
    """

    label = "Solar System Moons"

    # the data in the table probably don't need to be updated
    _expiration = np.inf * u.day

    def __init__(self, **kwargs):
        """
        Initialize a population of Solar System (main) planets.
        """
        PredefinedPopulation.__init__(self, **kwargs)
        self.color = "mediumblue"
        self.zorder = None
        self.s = 20
        self.respond_to_color = False
        self.exact = True
        self.marker = "s"

    def _download_raw_data(self, **kw):
        """
        Load the raw table of data for the Solar System.
        """

        # load the table of Solar System planets
        raw = ascii.read(initial_filename)

        # for debugging, hang on to the raw table as a hidden attribute
        self._raw = raw

        # a table of unstandardized planet properties
        return raw

    def create_standardardized(self, raw):
        """
        Create a standardized table of planet properties.
        """

        # start from the raw trimmed table
        t = raw

        s = QTable()
        s["name"] = t["satellite"]
        s["planet"] = t["planet"]
        s["code"] = t["code"]
        s["mass"] = (t["GM"] * u.km**3 / u.s**2 / con.G).to(u.Mearth)
        s["radius"] = (t["radius"] * u.km).to(u.Rearth)
        s["density"] = t["density"] * u.g / u.cm**3
        for k in ["lower", "upper"]:
            s[f"mass_uncertainty_{k}"] = (
                t["GM_uncertainty"] * u.km**3 / u.s**2 / con.G
            ).to(u.Mearth)
            s[f"radius_uncertainty_{k}"] = (t["radius_uncertainty"] * u.km).to(u.Rearth)
            s[f"density_uncertainty_{k}"] = t["density_uncertainty"] * u.g / u.cm**3

        # get stellar properties from the Solar System
        solar = SolarSystem()
        for k in ["mass", "radius", "teff", "age"]:
            s[f"stellar_{k}"] = solar.get(f"stellar_{k}")[0]
        s["hostname"] = "Sun"

        # get orbit properties from the Solar System planet hosts
        planets = solar.standard["name", "period", "semimajoraxis"]
        planets.columns[0].name = "planet"

        merged = join(planets, s, keys="planet")

        return merged
