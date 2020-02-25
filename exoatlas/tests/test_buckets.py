from exoatlas.imports import *
from exoatlas import *

def test_buckets():
    '''
    Make sure that our photon-counting tools plottables work,
    with their various possible telescope units.
    '''

    with mock.patch('builtins.input', return_value=""):
        t = TransitingExoplanets()
    DepthBrightness().build([t])
    BubblePanel(StellarBrightness, Depth).build([t]);
    BubblePanel(StellarBrightness(5*u.micron), Depth).build([t]);
    BubblePanel(StellarBrightness(5*u.micron, telescope='JWST'), Depth).build([t])
    for k in telescope_units:
        BubblePanel(StellarBrightness(telescope=k), Depth).build([t])

if __name__ == '__main__':
    outputs = {k.split('_')[-1]:v()
               for k, v in locals().items()
               if 'test_' in k}
