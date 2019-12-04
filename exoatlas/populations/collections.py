from ..imports import *
from .Exoplanets import Exoplanets
from .SolarSystem import SolarSystem

def get_exoplanets_by_method(methods='all'):
    '''
    Create a dictionary that contains a collection of exoplanet
    populations, separated by their discovery method.

    Returns
    -------

    '''

    # create a population of exoplanets
    e = Exoplanets()

    # use all the methodsx
    if methods == 'all':
        methods = np.unique(e.method)

    # create a dictionary of populations,
    pops = {}
    for m in methods:
        pops[m] = e[e.method == m]
        pops[m].label = m
        pops[m].color = None


    s = SolarSystem()
    s.color = 'black'
    pops['Solar System'] = s

    return pops
