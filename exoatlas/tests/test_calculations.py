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
        "impact_parameter",
    ]:
        e.get(k, distribution=False)
        e.get(k, distribution=True)
        e.get_uncertainty(k)
