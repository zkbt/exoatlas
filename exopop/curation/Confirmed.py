def correct(pop):

    """ I went through all the planets shown in the Dressing et al. (2015)
        mass-radius diagram to double the values that are present in the
        NASA Exoplanet Archive."""

    # parameters taken from the updates in Demory et al. (2015, *submitted*)
    pop.correct('55 Cnce',
                planet_radius=1.92,
                planet_radius_upper=0.08,
                planet_radius_lower=-0.08,
                a_over_rs=3.53,
                planet_mass=8.08,
                planet_mass_upper=0.31,
                planet_mass_lower=-0.31,
                b=0.36)

    # from the updates in Anglada-Escude et al. (2013)
    pop.correct('GJ 1214b',
                teff=3252.0,
                stellar_radius=0.211,
                planet_radius=2.72,
                planet_radius_upper=0.24,
                planet_radius_lower=-0.24,
                a_over_rs=14.62,
                rv_semiamplitude=10.9,
                planet_mass=6.19,
                planet_mass_upper=0.91,
                planet_mass_lower=-0.91,
                b=0.2,
                stellar_distance=14.55,
                stellar_distance_upper=0.13,
                stellar_distance_lower=-0.13)

    # from Knutson et al. (2014)
    pop.correct('HD 97658b',
                a_over_rs=26.24)

    # from Haywood et al. (2014), allowing eccentricity as they do in abstract
    pop.correct('CoRoT-7b',
                rv_semiamplitude=3.42,
                planet_mass=4.73,
                planet_mass_upper=0.95,
                planet_mass_lower=-0.95)

    # from Motalebi et al. (2015)
    pop.correct('HD 219134b',
                a_over_rs=10.53,
                planet_mass=4.46,
                planet_mass_upper=0.47,
                planet_mass_lower=-0.47)

    # from Carter et al.
    pop.correct('Kepler-36b',
                a_over_rs = 15.24)

    # the archive parameters look reasonable for
    # HIP 116454b
    # Kepler-10b
    # Kepler-10c
    # Kepler-93b
    #
