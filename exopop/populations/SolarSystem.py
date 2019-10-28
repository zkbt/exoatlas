from ..imports import *
from .Population import PredefinedPopulation
import astropy.units as u

initial_filename = directories['data'] + 'solarsystem.txt'

class SolarSystem(PredefinedPopulation):
    '''The Solar System, very crudely.'''

    # the data in the table probably don't need to be updated
    expiration = np.inf

    def __init__(self, **kwargs):
        '''
        Initialize a population of simulated TESS planets.
        '''
        PredefinedPopulation.__init__(self, label='Solar System', **kwargs)
        self.color = 'cornflowerblue'
        self.zorder = 1e10

    def load_raw(self):
        '''
        Load the raw table of data from the NASA Exoplanet Archive.
        '''

        # load the table of Solar System planets
        raw = ascii.read(initial_filename)

        # for debugging, hang on to the raw table as a hidden attribute
        self._raw = raw

        # a table of unstandardized planet properties
        return raw

    def trim_raw(self, raw):
        '''
        Trim bad/unnecessary rows out of a raw table of planet properties.
        '''

        # no trimming necessary
        trimmed = raw

        # for debugging, hang onto the trimmed table
        self._trimmed = trimmed

        # a trimmed table
        return self._trimmed


    def create_standard(self, trimmed):
        '''
        Create a standardized table of planet properties.
        It must at least contain the columns in
        `necessary_columns`.
        '''

        # start from the trimmed table
        t = trimmed

        # create an empty standardized table
        s = Table()
        s['name'] = t['name']

        # set up the Sun
        s['stellar_teff'] = 5780*u.K

        s['stellar_radius'] = 1.0*u.Rsun
        s['stellar_mass'] = 1.0*u.Msun

        # store the period, by default, in day
        s['period'] = (t['period']*u.year).to(u.day)

        # hide the stellar magnitudes
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
        s['discoverer'] = 'humans'
        s['transit_epoch'] = 0.0
        s['transit_duration'] = 0.0
        s['transit_depth'] = 0.0
        s['b'] = 0.0

        self.standard = s
        return s
    @property
    def distance(self):
        return np.nan*np.zeros_like(self.standard['stellar_distance'])
