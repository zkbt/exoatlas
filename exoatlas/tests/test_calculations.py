from .setup_tests import *

from exoatlas import *


def test_planet_calculations():
    """
    Do the planet-related population calculations work?
    """
    e = TransitingExoplanets()
    calculations = """semimajoraxis_from_period,
semimajoraxis_from_transit_scaled_semimajoraxis,
semimajoraxis,
scaled_semimajoraxis_from_semimajoraxis,
scaled_semimajoraxis,
eccentricity,
argument_of_periastron,
transit_impact_parameter_from_inclination,
transit_impact_parameter,
insolation,
relative_insolation,
log_relative_insolation,
relative_cumulative_xuv_insolation,
teq,
planet_luminosity,
transit_depth_from_radii,
scaled_radius_from_radii,
scaled_radius,
transit_duration_from_orbit,
transit_duration,
surface_gravity,
density,
escape_velocity,
orbital_velocity,
impact_velocity,
escape_parameter,
scale_height""".split(
        ",\n"
    )
    for k in tqdm(calculations):
        print(f"testing planet calculations for {k}")
        assert np.any(e.get(k, distribution=False) != 0)
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
        assert np.any(e.get_uncertainty(k) != 0)
