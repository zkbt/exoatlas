from .setup_tests import *

from exoatlas import *


def test_planet_calculations():
    """
    Do the planet-related population calculations work?
    """
    e = TransitingExoplanets()
    for k in [
        "semimajoraxis",
        "scaled_semimajoraxis",
        "eccentricity",
        "argument_of_periastron",
        "transit_impact_parameter",
    ]:
        assert np.any(e.get(k, distribution=False) != 0)
        assert np.any(e.get(k, distribution=True) != 0)
        assert np.any(e.get_uncertainty(k) != 0)


def test_stellar_calculations():
    """
    Do the planet-related population calculations work?
    """
    e = TransitingExoplanets()
    for k in [
        "stellar_luminosity",
        "stellar_luminosity_from_radius_and_teff",
        "distance_modulus",
    ]:
        assert np.any(e.get(k, distribution=False) != 0)
        assert np.any(e.get(k, distribution=True) != 0)
        assert np.any(e.get_uncertainty(k) != 0)
