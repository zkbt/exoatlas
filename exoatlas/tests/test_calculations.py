from .setup_tests import *

from exoatlas import *


def test_planet_calculations():
    """
    Do the planetary-related population calculations work?
    """
    calculations = """
semimajoraxis_from_period,
semimajoraxis,
scaled_semimajoraxis,
eccentricity,
argument_of_periastron,
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
scale_height
""".strip().split(
        ",\n"
    )
    test_calculations(calculations=calculations)


def test_stellar_calculations():
    """
    Do the stellar-related population calculations work?
    """
    calculations = """
stellar_luminosity_from_radius_and_teff,
stellar_luminosity,
distance_modulus
""".strip().split(
        ",\n"
    )
    test_calculations(calculations=calculations)


def test_observability_calculations():
    """
    Do the observability-related population calculations work?
    """
    calculations = """
angular_separation,
transmission_signal,
emission_signal,
reflection_signal,
stellar_brightness,
stellar_brightness_in_telescope_units,
depth_uncertainty,
depth_snr,
emission_snr,
reflection_snr,
transmission_snr
""".strip().split(
        ",\n"
    )
    test_calculations(calculations=calculations)


def test_calculations(calculations=[]):
    e = TransitingExoplanets()[::100]
    s = SolarSystem()
    for k in calculations:
        for p in [e, s]:
            print(f'trying to calculate "{k}" for {p}')
            p.get(k)
            p.get_uncertainty(k)
