# FIXME! right now it's not obvious except to a pro user what items
# are attributes and which are methods which need to be called as
# functions. We should try to make that more transparent/easy!

core_descriptions = {
    "name": "name of the planet/star/object",
    "ra": "Right Ascension of the system",
    "dec": "Declination of the system",
    "distance": "distance to the system",
}

core_stellar_descriptions = {
    "stellar_teff": "stellar effective temperature",
    "stellar_mass": "stellar mass",
    "stellar_radius": "stellar radius",
}

derived_stellar_descriptions = {
    "stellar_luminosity": "bolometric luminosity of the star",
    "distance_modulus": "apparent magnitude - absolute magnitude",
}


# transit_columns = [
#    "period",
#    "semimajoraxis",
#    "eccentricity",
#    "omega",
#    "inclination",
#    "transit_midpoint",
#    "transit_duration",
#    "transit_depth",
#    "stellar_teff",
#    "stellar_mass",
#    "stellar_radius",
#    "radius",
#    "mass",
#    "transit_ar",
#    "transit_b",
# ]

core_planet_descriptions = {
    "discovery_facility": "telescope/project that found this planet",
    "period": "orbital period of the planet",
    "impact_parameter": "impact parameter b",
    "eccentricity": "eccentricity",
    "omega": "argument of periastron",
    "radius": "planet radius",
    "mass": "planet mass",
}

core_transit_descriptions = {
    "transit_midpoint": "a transit midpoint",
    "transit_duration": "duration of the transit",
    "transit_depth": "fraction of starlight the planet blocks",
    "transit_ar": "(transit-derived) scaled orbital distance a/R*",
    "transit_b": "(transit-derived) impact parameter b",
}


# calculated_columns = [
#    "a_over_rs",
#    "b",
#    "insolation",
#    "relative_insolation",
#    "log_relative_insolation",
#    "teq",
#    "planet_luminosity",
#    "density",
#    "surface_gravity",
#    "distance_modulus",
#    "escape_velocity",
#    "escape_parameter",
#    "stellar_luminosity",
# ]

derived_planet_descriptions = {
    "inclination": "orbital inclination",
    "semimajoraxis": "the semimajor axis of the planet's orbit",
    "a_over_rs": "scaled orbital distance a/R*",
    "density": "density of the planet",
    "insolation": "bolometric energy flux the planet receives from its star",
    "relative_insolation": "insolation relative to Earth",
    "log_relative_insolation": "log10(insolation relative to Earth)",
    "teq": "equilibrium temperature of the planet (assuming 0 albedo)",
    "planet_luminosity": "power emitted by the planet (assuming 0 albedo)",
    "surface_gravity": "surface gravitational acceleration of the planet",
    "scale_height": "scale height of an H2-rich atmosphere",
    "escape_velocity": "escape velocity of the planet",
    "escape_parameter": "ratio of gravitational potential to thermal energy for an H atom",
    "transmission_signal": "transit depth of one scale height of atmosphere (a function of mean molecular weight)",
    "emission_signal": "thermal-emission eclipse depth (a function of wavelength)",
    "reflection_signal": "reflected-light eclipse depth (a function of albedo)",
}

column_descriptions = (
    core_descriptions
    | core_stellar_descriptions
    | derived_stellar_descriptions
    | core_planet_descriptions
    | core_transit_descriptions
)
derived_descriptions = (
    derived_stellar_descriptions
    | derived_planet_descriptions
    | derived_observability_descriptions
)
all_descriptions = core_basic_descriptions | derived_descriptions

basic_columns = list(core_basic_descriptions.keys())
core_columns = list(core_descriptions.keys())
derived_columns = list(derived_descriptions.keys())
all_columns = list(all_descriptions.keys())


def describe_columns():
    """
    Describe some of the common columns you might want to
    access in an exoplanet population.
    """
    N = max([len(x) for x in column_descriptions])
    f = "{:>" + str(N) + "} = {}"
    for k, v in column_descriptions.items():
        print(f.format(k, v))
