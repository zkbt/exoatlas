from .imports import *
from .Population import Population
import astropy.units as u

initial_filename = directories['data'] + 'solarsystem.txt'

class SolarSystem(PredefinedPopulation):
    '''The Solar System, very crudely.'''

    def __init__(self, **kwargs):
        '''Initialize a population of simulated TESS planets.'''
        Population.__init__(self, label='Solar System', **kwargs)
        self.color = 'cornflowerblue'
        self.zorder = 1e10

    def loadFromScratch(self):

        self.table = ascii.read(initial_filename)
        self.speak('loaded TESS simulated population from {0}'.format(initial_filename))

    def trim_raw(self):
        self.trimmed = self.table

    def create_standard(self):
        t = self.trimmed
        s = Table()
        s['name'] = t['name']
        s['kepid'] = None

        s['stellar_teff'] = 5780

        s['period'] = (t['period']*u.year/u.day).decompose()

        s['stellar_radius'] = 1.0
        s['stellar_mass'] = 1.0

        s['J'] = None
        s['V'] = None

        s['planet_radius'] = (t['radius']*u.km/u.Rearth).decompose().value
        s['planet_radius_lower'] = 0.0
        s['planet_radius_upper'] = 0.0
        assert((s['planet_radius'] < 20).all())
        s['planet_mass'] = (t['mass']*u.kg/u.Mearth).decompose().value
        s['planet_mass_lower'] = 0.0
        s['planet_mass_upper'] = 0.0

        flux = (s['period'])**(-4.0/3.0)
        semimajoraxis = t['period']**(2.0/3.0)
        s['a_over_r'] = (semimajoraxis*u.au/(1*u.Rsun)).decompose().value

        s['teq'] = s['stellar_teff']*(0.25*1.0/s['a_over_r']**2)**0.25

        s['rv_semiamplitude'] = None
        s['radius_ratio'] = (s['planet_radius']*u.Rearth/(s['stellar_radius']*u.Rsun)).decompose().value
        s['stellar_distance'] = 0.0
        s['ra'] = 0.0
        s['dec'] = 0.0

        self.standard = s

    @property
    def distance(self):
        return np.nan*np.zeros_like(self.standard['stellar_distance'])
