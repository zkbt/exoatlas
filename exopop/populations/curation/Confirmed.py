import numpy as np
import astropy.units as u

def correct(pop):


    # Grimm et al. masses
    pop.correct('TRAPPIST-1b',      planet_mass=  1.017,
                                    planet_mass_upper=0.154,
                                    planet_mass_lower= -0.143)

    pop.correct('TRAPPIST-1c',      planet_mass= 1.156,
                                    planet_mass_upper=0.142,
                                    planet_mass_lower= -0.131)

    pop.correct('TRAPPIST-1d',      planet_mass= 0.297,
                                    planet_mass_upper=0.039,
                                    planet_mass_lower= -0.035)


    pop.correct('TRAPPIST-1e',      planet_mass= 0.772,
                                    planet_mass_upper=0.079,
                                    planet_mass_lower= -0.075)


    pop.correct('TRAPPIST-1f',      planet_mass= 0.934,
                                    planet_mass_upper=0.080,
                                    planet_mass_lower=-0.078)
    pop.correct('TRAPPIST-1g',      planet_mass=1.148,
                                    planet_mass_upper=0.098,
                                    planet_mass_lower=-0.095)
    pop.correct('TRAPPIST-1h',      planet_mass=0.331,
                                    planet_mass_upper=0.056,
                                    planet_mass_lower=-0.049)

    pop.correct('GJ 1132b',
            a_over_r=16.54,
            planet_mass=1.62,
            planet_mass_upper=0.55,
            planet_mass_lower=-0.55,
			planet_radius=1.130,
			planet_radius_upper=0.056,
			planet_radius_lower=-0.056,
			stellar_teff=3270.0, b=0.38,
			stellar_radius=0.2105,
			stellar_mass=0.181,
            transit_epoch=2457184.55786)

    pop.correct('30 Ari Bb',
            planet_radius=np.nan)


    pop.correct('HD 17156b',
            transit_duration=0.1338, a_over_r=23.2, b=0.4)

    pop.correct('HD 17156b',
            transit_duration=0.1338, a_over_r=23.2, b=0.477)

    pop.correct('HD 80606b',
            transit_duration=0.4958, a_over_r=94, b=0.75)

    pop.correct('GJ 1132b',
                transit_duration=47.0/60.0/24.0)

    pop.correct('GJ 1214b',
                transit_duration=0.03620)

    pop.correct('GJ 436b',
                transit_duration=0.04227)

    pop.correct('HD 219134b',
            a_over_r=1/(0.03876*u.au/0.778/u.Rsun).decompose().value,
            planet_mass=4.74,
            planet_mass_upper=0.19,
            planet_mass_lower=-0.19)

    pop.correct('HD 219134c',
            a_over_r=1/(0.06530*u.au/0.778/u.Rsun).decompose().value,
            planet_mass=4.36,
            planet_mass_upper=0.22,
            planet_mass_lower=-0.22)

    pop.correct('LHS 1140b',
            transit_duration=2.0/24.0,
            a_over_r=101.0,
            b=0.16)

    pop.correct('GJ 3470b',
            transit_duration=0.07992)

    pop.correct('Kepler-51b',
            transit_duration=0.24154)

    pop.correct('Kepler-51c',
            transit_duration=0.1171)

    pop.correct('Kepler-51d',
            transit_duration=0.351992)

    pop.correct('Kepler-42b',
            transit_duration=0.022329)

    pop.correct('Kepler-42c',
            transit_duration=0.02)

    pop.correct('Kepler-42d',
            transit_duration=0.016554)

    pop.correct('WASP-43b',
            transit_duration=0.0483, a_over_r=4.8, b=0.66)

    pop.correct('K2-18b', transit_duration=0.1117)

    """I corrected the famous ones, which have weird choices in the archive."""
    pop.correct('HD 209458b',
                planet_mass=235.0,
                planet_mass_upper=19.0,
                planet_mass_lower=-19.0,
                stellar_teff=6065.0,
                a_over_r=8.9,
                b=0.5)


    pop.correct('HD 189733b',
                a_over_r=8.9,
                planet_radius=13.630,
                planet_radius_upper=0.269,
                planet_radius_lower=-0.269,
                planet_mass=369.303,
                planet_mass_upper=18.433,
                planet_mass_lower=-18.433,
                transit_epoch=2454279.436714)

    """ I went through all the planets shown in the Dressing et al. (2015)
        mass-radius diagram to double check the values that are present in the
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
                stellar_teff=3416.0)

    # from the updates in Anglada-Escude et al. (2013)
    pop.correct('GJ 1214b',
                stellar_teff=3252.0,
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

    # googling the ones that have messed up stellar_teffs
    pop.correct('Qatar-1b',
                stellar_teff=4861.)
    pop.correct('Kepler-450b',
                stellar_teff=6215.0)
    pop.correct('Kepler-450c',
                stellar_teff=6215.0)
    pop.correct('Kepler-450d',
                stellar_teff=6215.0)

    pop.correct('Kepler-449b',
                stellar_teff=5588.0)
    pop.correct('Kepler-449c',
                stellar_teff=5588.0)

    pop.correct('Kepler-414b',
                stellar_teff=5523.0)
    pop.correct('Kepler-414c',
                stellar_teff=5523.0)

    pop.correct('Kepler-415b',
                stellar_teff=4523.0)
    pop.correct('Kepler-415c',
                stellar_teff=4523.0)

    pop.correct('Kepler-416b',
                stellar_teff=5670.0)
    pop.correct('Kepler-416c',
                stellar_teff=5670.0)

    pop.correct('Kepler-417b',
                stellar_teff=5376.0)
    pop.correct('Kepler-417c',
                stellar_teff=5376.0)

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
