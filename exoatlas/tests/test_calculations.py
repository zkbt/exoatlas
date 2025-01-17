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


def test_add_columns():
    s = SolarSystem()
    new_column_name = "is_inhabited"
    new_column_data = (s.name() == "Earth") * 1
    new_column_uncertainty = (s.name() != "Earth") * 0.01
    s.add_column(
        name=new_column_name, data=new_column_data, uncertainty=new_column_uncertainty
    )
    s.is_inhabited()
    s.is_inhabited_uncertainty()


def test_add_calculations():
    s = SolarSystem()

    def f(self, distribution=False):
        """
        Surface Area (m)

        Calculate the surface area of a planet,
        based on its radius.
        """
        return 4 * np.pi * self.radius(distribution=distribution)

    s.add_calculation(name="surface_area", function=f)
    s.surface_area()
    s.surface_area_uncertainty()
