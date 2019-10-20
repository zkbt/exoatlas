from exopop.imports import *
from exopop.populations import *

def test_population():
    '''
    Can we make a population from scratch from a table?
    '''

    fake = Table({x:[0]*3 for x in necessary_columns}, masked=True)
    p = Population(standard=fake, label='fake')
    p.validate_columns()
    return p

def test_transitingexoplanets():
    '''
    Can we make a population of confirmed transiting exoplanets?
    '''
    p = TransitingExoplanets()
    p.validate_columns()
    return p

def test_exoplanets():
    '''
    Can we make a population of confirmed exoplanets?
    '''
    p = Exoplanets()
    p.validate_columns()
    return p

def test_indexing():
    '''
    Can we make a population of confirmed transiting exoplanets?
    '''
    p = TransitingExoplanets()

    # try different subsets
    a = p[[]]
    b = p[5]
    c = p[0:10]
    d = p[p.stellar_radius < 1.0]
    e = p['GJ 1214b']
    f = p[['GJ 1214b', 'GJ 1132b', 'LHS 1140b']]
    g = p[p.discoverer == 'Kepler']

    return a, b, c, d, e

if __name__ == '__main__':
    outputs = {k.split('_')[-1]:v()
               for k, v in locals().items()
               if 'test_' in k}
