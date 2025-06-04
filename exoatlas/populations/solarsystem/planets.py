from ...imports import *
from ..predefined import PredefinedPopulation

__all__ = ["SolarSystem", "SolarSystemDwarfPlanets"]


class SolarSystem(PredefinedPopulation):
    """The Solar System, very crudely."""

    _raw_mass_unit = u.Unit(1e24 * u.kg)
    _filename = "jpl-ssd-planets.txt"

    label = "Solar System"
    # the data in the table probably don't need to be updated
    _expiration = np.inf * u.day

    def __init__(self, **kwargs):
        """
        Initialize a population of Solar System (main) planets.
        """
        PredefinedPopulation.__init__(self, **kwargs)
        self.color = "cornflowerblue"
        self.zorder = 1e10
        self.s = 64
        self.respond_to_color = False
        self.exact = True
        self.marker = "s"

    def _download_raw_data(self, **kw):
        """
        Load the raw table of data for the Solar System.
        """

        # load the table of Solar System planets
        initial_filename = os.path.join(
            code_directory, f"populations/data/solarsystem/{self._filename}"
        )
        raw = ascii.read(initial_filename)

        # for debugging, hang on to the raw table as a hidden attribute
        self._raw = raw

        # a table of unstandardized planet properties
        return raw

    def _create_standardized(self, trimmed):
        """
        Create a standardized table of planet properties.
        """

        # start from the trimmed table
        t = QTable(trimmed)

        # create an empty standardized table
        s = QTable()
        s["name"] = t["name"]
        s["hostname"] = "Sun"

        # set up the Sun
        s["stellar_teff"] = 5780 * u.K
        s["stellar_radius"] = 1.0 * u.Rsun
        s["stellar_mass"] = 1.0 * u.Msun
        s["stellar_luminosity"] = 1 * u.Lsun
        s["stellar_age"] = 4.57 * u.Gyr

        # store the period, by default, in day
        s["period"] = (t["period"] * u.year).to(u.day)

        def add_with_uncertainty(
            their_name, our_name, their_unit=u.Unit(), our_unit=u.Unit()
        ):
            s[our_name] = (t[their_name] * their_unit).to(our_unit)
            try:
                s[f"{our_name}_uncertainty_{k}"] = (
                    t[f"{their_name}_uncertainty"] * their_unit
                ).to(our_unit)
            except KeyError:
                pass

        add_with_uncertainty("equatorial_radius", "radius_equatorial", u.km, u.Rearth)
        add_with_uncertainty("mean_radius", "radius", u.km, u.Rearth)
        add_with_uncertainty("mass", "mass", self._raw_mass_unit, u.Mearth)
        add_with_uncertainty("density", "density", u.g / u.cm**3, u.g / u.cm**3)
        add_with_uncertainty("rotational_period", "rotational_period", u.day, u.day)
        add_with_uncertainty("period", "period", u.year, u.day)
        add_with_uncertainty("absolute_magnitude", "absolute_magnitude", u.mag, u.mag)
        add_with_uncertainty("albedo_geometric", "albedo_geometric")
        add_with_uncertainty(
            "surface_gravity", "surface_gravity", u.m / u.s**2, u.m / u.s**2
        )
        add_with_uncertainty(
            "escape_velocity", "escape_velocity", u.km / u.s, u.km / u.s
        )
        try:
            add_with_uncertainty("eccentricity", "eccentricity")
            add_with_uncertainty(
                "argument_of_periastron", "argument_of_periastron", u.deg, u.deg
            )
            s["inclination"] = (90 - t["inclination"]) * u.deg

        except KeyError:
            s["eccentricity"] = np.nan
            s["argument_of_periastron"] = np.nan * u.deg
            s["inclination"] = 90 * u.deg

        # use Kepler's (actual) 3rd Law to get semimajor axis
        semimajoraxis = (s["period"].to(u.year).value) ** (2.0 / 3.0) * u.AU
        s["semimajoraxis"] = semimajoraxis
        # s["transit_scaled_semimajoraxis"] = (semimajoraxis / u.Rsun).decompose().value

        # equilibrium temperature assuming uniform redistribution + 0 albedo
        # s["teq"] = s["stellar_teff"] * (0.25 * 1.0 / s["transit_scaled_semimajoraxis"] ** 2) ** 0.25

        # some other planet parameters we might not need for the solar system
        # s["rv_semiamplitude"] = np.nan * u.m / u.s
        s["radius_ratio"] = (s["radius"] / s["stellar_radius"]).decompose().value
        s["distance"] = np.nan * u.pc
        s["ra"] = np.nan * 0.0 * u.deg
        s["dec"] = np.nan * 0.0 * u.deg
        # s["discovery_facility"] = "humans"
        # s["transit_midpoint"] = np.nan * u.day
        # s["transit_duration"] = np.nan * u.day
        # s["transit_depth"] = (s["radius"] / s["stellar_radius"]).decompose() ** 2
        # s["transit_impact_parameter"] = 0.0
        # s["inclination"] = 90 * u.deg

        s

        self.standard = s
        return s


class SolarSystemDwarfPlanets(SolarSystem):
    _raw_mass_unit = u.Unit(1e18 * u.kg)
    _filename = "jpl-ssd-dwarf-planets.txt"

    label = "Solar System Dwarf Planets"

    def __init__(self, **kwargs):
        """
        Initialize a population of Solar System (dwarf) planets.
        """
        PredefinedPopulation.__init__(self, **kwargs)
        self.color = "royalblue"
        self.s = 32
        self.respond_to_color = False
        self.exact = True
        self.marker = "s"
