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

def test_solarsystem():
    '''
    Can we make a population of Solar System planets?
    '''
    p = SolarSystem()
    p.validate_columns()
    return p

def test_transitingexoplanets():
    '''
    Can we make a population of confirmed transiting exoplanets?
    '''
    p = TransitingExoplanets(skip_update=True)
    p.validate_columns()
    return p

def test_exoplanets():
    '''
    Can we make a population of confirmed exoplanets?
    '''
    p = Exoplanets(skip_update=True)
    p.validate_columns()
    return p

def test_kepler():
    '''
    Can we make a population of confirmed Kepler planets?
    '''
    p = Kepler(skip_update=True)
    p.validate_columns()
    return p


def test_nonkepler():
    '''
    Can we make a population of confirmed non-Kepler planets?
    '''
    p = NonKepler(skip_update=True)
    p.validate_columns()
    return p

def test_tess():
    '''
    Can we make a population of confirmed TESS planets?
    '''
    p = TESS(skip_update=True)
    p.validate_columns()
    return p

def test_nontess():
    '''
    Can we make a population of confirmed TESS planets?
    '''
    p = NonTESS(skip_update=True)
    p.validate_columns()
    return p

def test_space():
    '''
    Can we make a population of confirmed space-discovered planets?
    '''
    p = Space(skip_update=True)
    p.validate_columns()
    return p

def test_ground():
    '''
    Can we make a population of confirmed ground-discovered planets?
    '''
    p = Ground(skip_update=True)
    p.validate_columns()
    return p

def test_goodmass():
    '''
    Can we make a population of confirmed planets with good masses?
    '''
    p = GoodMass(skip_update=True)
    p.validate_columns()
    return p

def test_badmass():
    '''
    Can we make a population of confirmed planets with bad masses?
    '''
    p = BadMass(skip_update=True)
    p.validate_columns()
    return p

def test_indexing():
    '''
    Can we make a population of confirmed transiting exoplanets?
    '''
    p = TransitingExoplanets(skip_update=True)

    # try different subsets
    a = p[[]]
    b = p[5]
    c = p[0:10]
    d = p[p.stellar_radius < 1.0*u.Rsun]
    e = p['GJ 1214b']
    f = p[['GJ 1214b', 'LHS 1140b']] #'GJ 1132b',
    g = p[p.discoverer == 'Kepler']

    return a, b, c, d, e

def test_attributes():
    '''
    Can we make a population of Solar System planets?
    '''
    p = SolarSystem()


    print(p.color)
    p.alpha = 0.5
    assert(p.plotkw['alpha'] == 0.5)


    for k in necessary_columns:
        getattr(p, k)

    p.b
    p.transit_duration
    p.insolation
    p.stellar_luminosity
    p.a_over_rs
    p.teq
    p.mu
    p.surface_gravity
    p.escape_velocity
    p.escape_parameter
    p.scale_height
    return p

if __name__ == '__main__':
    outputs = {k.split('_')[-1]:v()
               for k, v in locals().items()
               if 'test_' in k}
