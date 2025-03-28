"""
Define some wrappers to help generate a collection of different populations
that might be useful to compare with one another. This is mostly a tool

"""

from ..imports import *
from .exoplanets import *
from .solarsystem import *


def get_transiting_and_nontransiting_exoplanets():
    e = Exoplanets()
    transit = e[e.detected_in_transit() != 0]
    transit.label = "Transiting Exoplanets"
    nontransit = e[e.detected_in_transit() == 0]
    nontransit.label = "Non-transiting Exoplanets"
    return dict(nontransit=nontransit, transit=transit)


def get_exoplanets_by_method(methods="all", include_solar_system=True):
    """
    Create a dictionary that contains a collection of exoplanet
    populations, separated by their discovery method.

    Returns
    -------

    pops : dict
        A dictionary of planet populations,
        grouped by discovery method.
    """

    # create a population of exoplanets
    e = Exoplanets()

    # use all the methods
    if methods == "all":
        methods = np.unique(e.discovery_method())

    # create a dictionary of populations,
    pops = {}
    for m in methods:
        pops[m] = e[e.discovery_method() == m]
        pops[m].label = m
        pops[m].color = None

    if include_solar_system:
        s = SolarSystem()
        s.color = "black"
        pops["Solar System"] = s

    return pops


def get_solar_system_objects():
    """
    Create a dictionary that contains a collection of
    Solar System populations, grouped by category.

    Returns
    -------
    pops : dict
        A dictionary of Solar System populations,
        grouped by (JPL/SSD) category.
    """
    pops = dict(
        major=SolarSystem(),
        dwarf=SolarSystemDwarfPlanets(),
        moons=SolarSystemMoons(),
        minor=SolarSystemMinorPlanets(),
    )
    return pops


def get_transiting_exoplanets_by_mass_precision(sigma=5):
    """
    Create a dictionary that contains a collection of exoplanet
    populations, separated by whether their masses are good or bad.

    Parameters
    ----------
    sigma : float
        At how many sigma should the mass be
        distinguishable from 0? The maximum
        fractional uncertainty on planets
        classified as 'good' will be 1/sigma.
    Returns
    -------
    pops : dict
        A dictionary of planet populations,
        grouped by mass precision.
    """
    # create a dictionary of populations,
    pops = dict(good=GoodMass(sigma=sigma), bad=BadMass(sigma=sigma))
    return pops


def get_exoplanets_by_teff():
    """
    Create a dictionary that contains a collection of exoplanet
    populations, separated by their stellar effective temperature.

    Returns
    -------
    pops : dict
        A dictionary of planet populations,
        grouped by stellar effective temperature.
    """

    e = Exoplanets()
    # spectral type boundaries from https://www.pas.rochester.edu/~emamajek/EEM_dwarf_UBVIJHK_colors_Teff.txt
    upper_temperature_limits = dict(
        Y=480,
        T=1312,
        L=2310,
        M=3890,
        K=5325,
        G=5960,
        F=7310,
        A=10050,
        B=31650,
        O=np.inf,
    )

    pops = {}
    for k, t in upper_temperature_limits.items():
        is_cool_enough = e.stellar_teff() < t * u.K
        try:
            pops[k] = e[is_cool_enough]
            pops[k].label = k
            pops[k].color = None
        except IndexError:
            pass
        e = e[is_cool_enough == False]

    # reverse order to get OBAFGKMLTY
    items = list(pops.items())
    items.reverse()
    return dict(items)


def get_all_planets():
    p = get_transiting_and_nontransiting_exoplanets() | get_solar_system_objects()
    p["transit"].color = "black"
    p["nontransit"].color = "coral"
    p["major"].annotate_planets = True
    p["dwarf"].annotate_planets = True
    return p
