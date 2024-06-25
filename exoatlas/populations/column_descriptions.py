# FIXME! right now it's not obvious except to a pro user what items
# are attributes and which are methods which need to be called as
# functions. We should try to make that more transparent/easy!

core_stellar_descriptions = {
    "name": "name of the planet/star/object",
    "ra": "Right Ascension of the system",
    "dec": "Declination of the system",
    "distance": "distance to the system",
    "stellar_teff": "stellar effective temperature",
    "stellar_mass": "stellar mass",
    "stellar_radius": "stellar radius",
    "stellar_luminosity": "luminosity of the star",
}

derived_stellar_descriptions = {
    "stellar_brightness": "photon flux from the star at Earth (a function of wavelength)",
    "distance_modulus": "apparent magnitude - absolute magnitude",
}

core_planet_descriptions = {
    "discovery_facility": "telescope/project that found this planet",
    "period": "orbital period of the planet",
    "b": "impact parameter b",
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

derived_planet_descriptions = {
    "semimajoraxis": "the semimajor axis of the planet's orbit",
    "a_over_rs": "scaled orbital distance a/R*",
    "density": "density of the planet",
    "insolation": "bolometric energy flux the planet receives from its star",
    "relative_insolation": "insolation relative to Earth",
    "teq": "equilibrium temperature of the planet (assuming 0 albedo)",
    "surface_gravity": "surface gravity of the planet",
    "scale_height": "scale height of an H2-rich atmosphere",
    "escape_velocity": "escape velocity of the planet",
    "escape_parameter": "ratio of gravitational potential to thermal energy for an H atom",
    "transmission_signal": "transit depth of one scale height of atmosphere (a function of mean molecular weight)",
    "emission_signal": "thermal-emission eclipse depth (a function of wavelength)",
    "reflection_signal": "reflected-light eclipse depth (a function of albedo)",
}

column_descriptions = (
    core_stellar_descriptions
    | derived_stellar_descriptions
    | core_planet_descriptions
    | core_transit_descriptions
    | derived_planet_descriptions
)


def describe_columns():
    """
    Describe some of the common columns you might want to
    access in an exoplanet population.
    """
    N = max([len(x) for x in column_descriptions])
    f = "{:>" + str(N) + "} = {}"
    for k, v in column_descriptions.items():
        print(f.format(k, v))
