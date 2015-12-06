import numpy as np
def correct(pop):

    """ I went through all the planets shown in the Dressing et al. (2015)
        mass-radius diagram to double the values that are present in the
        NASA Exoplanet Archive."""

    # parameters taken from the updates in Demory et al. (2015, *submitted*)
    '''pop.correct('55 Cnce',
                planet_radius=1.92,
                planet_radius_upper=0.08,
                planet_radius_lower=-0.08,
                a_over_r=3.53,
                planet_mass=8.08,
                planet_mass_upper=0.31,
                planet_mass_lower=-0.31,
                b=0.36)'''

    # from Gillon et al. (2012) (and some from Demory et al.)
    pop.correct('55 Cnce',
                planet_radius=2.17,
                planet_radius_upper=0.1,
                planet_radius_lower=-0.1,
                a_over_r=3.53,
                planet_mass=8.08,
                planet_mass_upper=0.31,
                planet_mass_lower=-0.31,
                b=0.459)

    # from von Braun et al. (2012)
    pop.correct('GJ 436b',
                teff=3416.0)

    # from the updates in Anglada-Escude et al. (2013)
    pop.correct('GJ 1214b',
                teff=3252.0,
                stellar_radius=0.211,
                planet_radius=2.72,
                planet_radius_upper=0.24,
                planet_radius_lower=-0.24,
                a_over_r=14.62,
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
                a_over_r=26.24)

    # from Haywood et al. (2014), allowing eccentricity as they do in abstract
    pop.correct('CoRoT-7b',
                rv_semiamplitude=3.42,
                planet_mass=4.73,
                planet_mass_upper=0.95,
                planet_mass_lower=-0.95)

    # from Motalebi et al. (2015)
    pop.correct('HD 219134b',
                a_over_r=10.53,
                planet_mass=4.46,
                planet_mass_upper=0.47,
                planet_mass_lower=-0.47)

    # from Carter et al.
    pop.correct('Kepler-36b',
                a_over_r = 15.24)

    # the archive parameters look reasonable for
    # HIP 116454b
    # Kepler-10b
    # Kepler-10c
    # Kepler-93b

    """ Then, I corrected things that popped up as weird non-physical outliers
        in various slices of parameter space, trying to track down real values
        whenever possible."""


    # Kepler-11's radius was listed as 0.07!
    pop.correct('Kepler-11b',
                stellar_radius=1.053)
    pop.correct('Kepler-11c',
                stellar_radius=1.053)
    pop.correct('Kepler-11d',
                stellar_radius=1.053)
    pop.correct('Kepler-11e',
                stellar_radius=1.053)
    pop.correct('Kepler-11f',
                stellar_radius=1.053)
    pop.correct('Kepler-11g',
                stellar_radius=1.053)

    # googling the ones that have messed up Teffs
    pop.correct('Qatar-1b',
                teff=4861.)
    pop.correct('Kepler-450b',
                teff=6215.0)
    pop.correct('Kepler-450c',
                teff=6215.0)
    pop.correct('Kepler-450d',
                teff=6215.0)

    pop.correct('Kepler-449b',
                teff=5588.0)
    pop.correct('Kepler-449c',
                teff=5588.0)

    pop.correct('Kepler-414b',
                teff=5523.0)
    pop.correct('Kepler-414c',
                teff=5523.0)

    pop.correct('Kepler-415b',
                teff=4523.0)
    pop.correct('Kepler-415c',
                teff=4523.0)

    pop.correct('Kepler-416b',
                teff=5670.0)
    pop.correct('Kepler-416c',
                teff=5670.0)

    pop.correct('Kepler-417b',
                teff=5376.0)
    pop.correct('Kepler-417c',
                teff=5376.0)

    # from discovery paper
    pop.correct('WASP-100b',
                stellar_radius=2.0)

    # supposedly HIPPARCOS, via exoplanets.org
    pop.correct('HD80606b',
                stellar_distance=58,
                stellar_distance_upper=30,
                stellar_distance_lower=-14.8)

    # from discovery paper
    pop.correct('Kepler-128b',
                planet_mass=np.nan)
    pop.correct('Kepler-128c',
                planet_mass=np.nan)
