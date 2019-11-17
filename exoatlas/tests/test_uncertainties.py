from exoatlas.imports import *
from exoatlas.populations import *

def test_uncertainties():
    '''
    Can we pull out c
    '''

    p = SolarSystem()

    uncertainty = p.uncertainty('planet_radius')
    assert(np.all(uncertainty == 0*u.Rearth))

    upper, lower = p.uncertainty_lowerupper('planet_radius')
    assert(np.all(lower == 0*u.Rearth))
    assert(np.all(upper == 0*u.Rearth))

    bad = p.uncertainty('stellar_distance')
    assert(np.all(np.isnan(bad)))
    return p

if __name__ == '__main__':
    outputs = {k.split('_')[-1]:v()
               for k, v in locals().items()
               if 'test_' in k}
