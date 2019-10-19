from .imports import *
from .Population import Population


initial_filename = directories['data'] + 'TESSsimulations.tsv'

class PredictedTESS(Population):
    '''TESS population object contains a simulated TESS planet yield
        from 200,000 two-minute cadence postage stamps,
        as calculated by Peter Sullivan et al. (2015)'''

    def __init__(self, **kwargs):
        '''Initialize a population of simulated TESS planets.'''
        Population.__init__(self, label='Predicted TESS', **kwargs)
        self.color = 'darkorange'
        self.zorder = 1

    def loadFromScratch(self):

        self.table = ascii.read(initial_filename)
        self.speak('loaded TESS simulated population from {0}'.format(initial_filename))

    def trimRaw(self):
        self.trimmed = self.table

    def createStandard(self):
        t = self.trimmed
        s = Table()
        s['name'] = ['tess{0:04}i'.format(i) for i in range(len(t))]
        s['kepid'] = None
        s['period'] = t['P']
        s['teff'] = t['Teff']
        s['stellar_radius'] = t['Rstar']
        s['J'] = t['Jmag']
        s['V'] = t['Vmag']

        s['planet_radius'] = t['Rplanet']
        s['planet_radius_lower'] = np.zeros_like(t['Rplanet']) + 0.0001
        s['planet_radius_upper'] = np.zeros_like(t['Rplanet']) + 0.0001

        s['planet_mass'] = np.nan + np.zeros_like(t['Rplanet'])
        s['planet_mass_lower'] = np.nan + np.zeros_like(t['Rplanet'])
        s['planet_mass_upper'] = np.nan + np.zeros_like(t['Rplanet'])

        s['teq'] = 280.0*t['S/SEarth']**0.25
        s['a_over_r'] = 0.5*(s['teff']/s['teq'])**2
        s['rv_semiamplitude'] = t['K']
        s['radius_ratio'] = t['Rplanet']*craftroom.units.Rearth/(t['Rstar']*craftroom.units.Rsun)
        s['stellar_distance'] = 10*10**(0.2*t['Dist']) + np.nan
        s['ra'] = t['RA']
        s['dec'] = t['Dec']

        self.standard = s
