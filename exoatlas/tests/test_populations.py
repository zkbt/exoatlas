from exoatlas.imports import *
from exoatlas.populations import *

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
    with mock.patch('builtins.input', return_value=""):
        p = TransitingExoplanets()
    p.validate_columns()
    return p

def test_exoplanets():
    '''
    Can we make a population of confirmed exoplanets?
    '''
    with mock.patch('builtins.input', return_value=""):
        p = Exoplanets()
    p.validate_columns()
    return p

def test_subsets():
    '''
    Can we make a population of confirmed Kepler planets?
    '''
    with mock.patch('builtins.input', return_value=""):
        for x in [Kepler,
                  NonKepler,
                  TESS,
                  NonTESS,
                  Space,
                  Ground,
                  GoodMass,
                  BadMass]:
            p = x()
            p.validate_columns()

        with pytest.raises(NotImplementedError):
            q = ExoplanetsSubset()

def test_tess():
    '''
    Can we make a population of confirmed TESS planets?
    '''
    with mock.patch('builtins.input', return_value=""):
        p = TESS()
    p.validate_columns()
    return p



def test_indexing():
    '''
    Can we make a population of confirmed transiting exoplanets?
    '''
    with mock.patch('builtins.input', return_value=""):
        p = TransitingExoplanets()

    # try different subsets
    a = p[[]]
    b = p[5]
    c = p[0:10]
    with np.errstate(invalid='ignore'):
        d = p[p.stellar_radius < 1.0*u.Rsun]
    e = p['GJ 1214b']
    f = p[['GJ 1214b', 'LHS 1140b']] #'GJ 1132b',
    g = p[p.discoverer == 'Kepler']
    h = p['TRAPPIST-1b']
    i = p['TRAPPIST-1']
    j = e + h
    k = f - e
    l = p.create_subset_by_name('GJ1214b')
    m = p.create_subset_by_hostname('GJ1214')
    coordinates = SkyCoord(e.ra, e.dec)
    n = p.create_subset_by_position(coordinates)
    return a, b, c, d, e, f, g, h, i, j, k, l, m, n



def test_attributes():
    '''
    Can we make a population of Solar System planets?
    '''
    p = SolarSystem()
    #for p in [SolarSystem(), TransitingExoplanets()]:


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
    p.density
    return p

if __name__ == '__main__': # pragma: no cover
    outputs = {k.split('_')[-1]:v()
               for k, v in locals().items()
               if 'test_' in k}
