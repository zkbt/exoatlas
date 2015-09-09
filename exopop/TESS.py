from imports import *
from Population import Population


initial_filename = directories['data'] + 'TESSsimulations.tsv'

class TESS(Population):
    '''TESS population object contains a simulated TESS planet yield
        from 200,000 two-minute cadence postage stamps,
        as calculated by Peter Sullivan et al. (2015)'''

    def __init__(self, **kwargs):
        '''Initialize a population of simulated TESS planets.'''
        Population.__init__(self, label='Predicted TESS', **kwargs)

    def loadFromScratch(self):

        self.table = astropy.io.ascii.read(initial_filename)
        self.speak('loaded TESS simulated population from {0}'.format(initial_filename))

    def trimRaw(self):
        self.trimmed = self.table

    def createStandard(self):
        t = self.trimmed
        s = astropy.table.Table()
        s['name'] = ['tess{0:04}i'.format(i) for i in range(len(t))]
        s['kepid'] = None
        s['period'] = t['P']
        s['teff'] = t['Teff']
        s['stellar_radius'] = t['Rstar']
        s['J'] = t['Jmag']
        s['V'] = t['Vmag']

        s['planet_radius'] = t['Rplanet']

        s['teq'] = 280.0*t['S/SEarth']**0.25
        s['a_over_r'] = 0.5*(s['teff']/s['teq'])**2
        s['rv_semiamplitude'] = t['K']
        s['radius_ratio'] = t['Rplanet']*zachopy.units.Rearth/(t['Rstar']*zachopy.units.Rsun)
        s['stellar_distance'] = 10*10**(0.2*t['Dist']) + np.nan
        s['ra'] = t['RA']
        s['dec'] = t['Dec']

        self.standard = s
