column_descriptions = {
    "name": "name of the planet",
    "ra": "Right Ascension of the system",
    "dec": "Declination of the system",
    "distance": "distance to the system",
    "distance_modulus": "apparent magnitude - absolute magnitude",
    "discovery_facility": "telescope/project that found this planet",
    "stellar_teff": "stellar effective temperature",
    "stellar_mass": "stellar mass",
    "stellar_radius": "stellar radius",
    "stellar_luminosity": "luminosity of the star",
    "stellar_brightness": "photon flux from the star at Earth (a function of wavelength)",
    "period": "orbital period of the planet",
    "semimajoraxis": "the semimajor axis of the planet's orbit",
    "a_over_rs": "scaled orbital distance a/R*",
    "b": "impact parameter b",
    "eccentricity": "eccentricity",
    "omega": "argument of periastron",
    "radius": "planet radius",
    "mass": "planet mass",
    "density": "density of the planet",
    "insolation": "bolometric energy flux the planet receives from its star",
    "relative_insolation": "insolation relative to Earth",
    "teq": "equilibrium temperature of the planet (assuming 0 albedo)",
    "surface_gravity": "surface gravity of the planet",
    "scale_height": "scale height of an H2-rich atmosphere",
    "escape_velocity": "escape velocity of the planet",
    "escape_parameter": "ratio of gravitational potential to thermal energy for an H atom",
    "transit_midpoint": "a transit midpoint",
    "transit_duration": "duration of the transit",
    "transit_depth": "fraction of starlight the planet blocks",
    "transit_ar": "(transit-derived) scaled orbital distance a/R*",
    "transit_b": "(transit-derived) impact parameter b",
    "transmission_signal": "transit depth of one scale height of atmosphere",
    "emission_signal": "thermal-emission eclipse depth (a function of wavelength)",
    "reflection_signal": "reflected-light eclipse depth (for an albedo of 1)",
}


def describe_columns():
    """
    Describe some of the common columns you might want to
    access in an exoplanet population.
    """
    N = max([len(x) for x in column_descriptions])
    f = "{:>" + str(N) + "} = {}"
    for k, v in column_descriptions.items():
        print(f.format(k, v))
