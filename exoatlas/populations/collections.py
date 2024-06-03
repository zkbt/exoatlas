from ..imports import *
from .Exoplanets import Exoplanets
from .SolarSystem import SolarSystem


def get_exoplanets_by_method(methods="all", include_solar_system=True):
    """
    Create a dictionary that contains a collection of exoplanet
    populations, separated by their discovery method.

    Returns
    -------

    pops : dict
        A dictionary of planet populations, grouped by
        keys indicating which method was used to find them.
    """

    # create a population of exoplanets
    e = Exoplanets()

    # use all the methods
    if methods == "all":
        methods = np.unique(e.method)

    # create a dictionary of populations,
    pops = {}
    for m in methods:
        pops[m] = e[e.method == m]
        pops[m].label = m
        pops[m].color = None

    if include_solar_system:
        s = SolarSystem()
        s.color = "black"
        pops["Solar System"] = s

    return pops


