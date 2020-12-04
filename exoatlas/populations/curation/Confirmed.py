import numpy as np
import astropy.units as u

def update_planet(pop):


    # Grimm et al. masses
    pop.update_planet('TRAPPIST-1b',      mass=  1.017,
                                    mass_uncertainty_upper=0.154,
                                    mass_uncertainty_lower= -0.143)

    pop.update_planet('TRAPPIST-1c',      mass= 1.156,
                                    mass_uncertainty_upper=0.142,
                                    mass_uncertainty_lower= -0.131)

    pop.update_planet('TRAPPIST-1d',      mass= 0.297,
                                    mass_uncertainty_upper=0.039,
                                    mass_uncertainty_lower= -0.035)


    pop.update_planet('TRAPPIST-1e',      mass= 0.772,
                                    mass_uncertainty_upper=0.079,
                                    mass_uncertainty_lower= -0.075)


    pop.update_planet('TRAPPIST-1f',      mass= 0.934,
                                    mass_uncertainty_upper=0.080,
                                    mass_uncertainty_lower=-0.078)
    pop.update_planet('TRAPPIST-1g',      mass=1.148,
                                    mass_uncertainty_upper=0.098,
                                    mass_uncertainty_lower=-0.095)
    pop.update_planet('TRAPPIST-1h',      mass=0.331,
                                    mass_uncertainty_upper=0.056,
                                    mass_uncertainty_lower=-0.049)

    pop.update_planet('GJ 1132b',
            a_over_r=16.54,
            mass=1.62,
            mass_uncertainty_upper=0.55,
            mass_uncertainty_lower=-0.55,
			radius=1.130,
			radius_uncertainty_upper=0.056,
			radius_uncertainty_lower=-0.056,
			stellar_teff=3270.0, b=0.38,
			stellar_radius=0.2105,
			stellar_mass=0.181,
            transit_epoch=2457184.55786)

    pop.update_planet('30 Ari Bb',
            radius=np.nan)


    pop.update_planet('HD 17156b',
            transit_duration=0.1338, a_over_r=23.2, b=0.4)

    pop.update_planet('HD 17156b',
            transit_duration=0.1338, a_over_r=23.2, b=0.477)

    pop.update_planet('HD 80606b',
            transit_duration=0.4958, a_over_r=94, b=0.75)

    pop.update_planet('GJ 1132b',
                transit_duration=47.0/60.0/24.0)

    pop.update_planet('GJ 1214b',
                transit_duration=0.03620)

    pop.update_planet('GJ 436b',
                transit_duration=0.04227)

    pop.update_planet('HD 219134b',
            a_over_r=1/(0.03876*u.au/0.778/u.Rsun).decompose().value,
            mass=4.74,
            mass_uncertainty_upper=0.19,
            mass_uncertainty_lower=-0.19)

    pop.update_planet('HD 219134c',
            a_over_r=1/(0.06530*u.au/0.778/u.Rsun).decompose().value,
            mass=4.36,
            mass_uncertainty_upper=0.22,
            mass_uncertainty_lower=-0.22)

    pop.update_planet('LHS 1140b',
            transit_duration=2.0/24.0,
            a_over_r=101.0,
            b=0.16)

    pop.update_planet('GJ 3470b',
            transit_duration=0.07992)

    pop.update_planet('Kepler-51b',
            transit_duration=0.24154)

    pop.update_planet('Kepler-51c',
            transit_duration=0.1171)

    pop.update_planet('Kepler-51d',
            transit_duration=0.351992)

    pop.update_planet('Kepler-42b',
            transit_duration=0.022329)

    pop.update_planet('Kepler-42c',
            transit_duration=0.02)

    pop.update_planet('Kepler-42d',
            transit_duration=0.016554)

    pop.update_planet('WASP-43b',
            transit_duration=0.0483, a_over_r=4.8, b=0.66)

    pop.update_planet('K2-18b', transit_duration=0.1117)

    """I corrected the famous ones, which have weird choices in the archive."""
    pop.update_planet('HD 209458b',
                mass=235.0,
                mass_uncertainty_upper=19.0,
                mass_uncertainty_lower=-19.0,
                stellar_teff=6065.0,
                a_over_r=8.9,
                b=0.5)


    pop.update_planet('HD 189733b',
                a_over_r=8.9,
                radius=13.630,
                radius_uncertainty_upper=0.269,
                radius_uncertainty_lower=-0.269,
                mass=369.303,
                mass_uncertainty_upper=18.433,
                mass_uncertainty_lower=-18.433,
                transit_epoch=2454279.436714)

    """ I went through all the planets shown in the Dressing et al. (2015)
        mass-radius diagram to double check the values that are present in the
        NASA Exoplanet Archive."""

    # parameters taken from the updates in Demory et al. (2015, *submitted*)
    '''pop.update_planet('55 Cnce',
                radius=1.92,
                radius_uncertainty_upper=0.08,
                radius_uncertainty_lower=-0.08,
                a_over_r=3.53,
                mass=8.08,
                mass_uncertainty_upper=0.31,
                mass_uncertainty_lower=-0.31,
                b=0.36)'''



    # from Gillon et al. (2012) (and some from Demory et al.)
    pop.update_planet('55 Cnce',
                radius=2.17,
                radius_uncertainty_upper=0.1,
                radius_uncertainty_lower=-0.1,
                a_over_r=3.53,
                mass=8.08,
                mass_uncertainty_upper=0.31,
                mass_uncertainty_lower=-0.31,
                b=0.459)

    # from von Braun et al. (2012)
    pop.update_planet('GJ 436b',
                stellar_teff=3416.0)

    # from the updates in Anglada-Escude et al. (2013)
    pop.update_planet('GJ 1214b',
                stellar_teff=3252.0,
                stellar_radius=0.211,
                radius=2.72,
                radius_uncertainty_upper=0.24,
                radius_uncertainty_lower=-0.24,
                a_over_r=14.62,
                rv_semiamplitude=10.9,
                mass=6.19,
                mass_uncertainty_upper=0.91,
                mass_uncertainty_lower=-0.91,
                b=0.2,
                distance=14.55,
                distance_uncertainty_upper=0.13,
                distance_uncertainty_lower=-0.13)

    # from Knutson et al. (2014)
    pop.update_planet('HD 97658b',
                a_over_r=26.24)

    # from Haywood et al. (2014), allowing eccentricity as they do in abstract
    pop.update_planet('CoRoT-7b',
                rv_semiamplitude=3.42,
                mass=4.73,
                mass_uncertainty_upper=0.95,
                mass_uncertainty_lower=-0.95)

    # from Carter et al.
    pop.update_planet('Kepler-36b',
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
    pop.update_planet('Kepler-11b',
                stellar_radius=1.053)
    pop.update_planet('Kepler-11c',
                stellar_radius=1.053)
    pop.update_planet('Kepler-11d',
                stellar_radius=1.053)
    pop.update_planet('Kepler-11e',
                stellar_radius=1.053)
    pop.update_planet('Kepler-11f',
                stellar_radius=1.053)
    pop.update_planet('Kepler-11g',
                stellar_radius=1.053)

    # googling the ones that have messed up stellar_teffs
    pop.update_planet('Qatar-1b',
                stellar_teff=4861.)
    pop.update_planet('Kepler-450b',
                stellar_teff=6215.0)
    pop.update_planet('Kepler-450c',
                stellar_teff=6215.0)
    pop.update_planet('Kepler-450d',
                stellar_teff=6215.0)

    pop.update_planet('Kepler-449b',
                stellar_teff=5588.0)
    pop.update_planet('Kepler-449c',
                stellar_teff=5588.0)

    pop.update_planet('Kepler-414b',
                stellar_teff=5523.0)
    pop.update_planet('Kepler-414c',
                stellar_teff=5523.0)

    pop.update_planet('Kepler-415b',
                stellar_teff=4523.0)
    pop.update_planet('Kepler-415c',
                stellar_teff=4523.0)

    pop.update_planet('Kepler-416b',
                stellar_teff=5670.0)
    pop.update_planet('Kepler-416c',
                stellar_teff=5670.0)

    pop.update_planet('Kepler-417b',
                stellar_teff=5376.0)
    pop.update_planet('Kepler-417c',
                stellar_teff=5376.0)

    # from discovery paper
    pop.update_planet('WASP-100b',
                stellar_radius=2.0)

    # supposedly HIPPARCOS, via exoplanets.org
    pop.update_planet('HD80606b',
                distance=58,
                distance_uncertainty_upper=30,
                distance_uncertainty_lower=-14.8)

    # from discovery paper
    pop.update_planet('Kepler-128b',
                mass=np.nan)
    pop.update_planet('Kepler-128c',
                mass=np.nan)
