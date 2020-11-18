from ..imports import *
from .Population import PredefinedPopulation
import astropy.units as u

__all__ = ['SolarSystem']

initial_filename = resource_filename(__name__,
                                     'data/solarsystem/solarsystem.txt')

class SolarSystem(PredefinedPopulation):
    '''The Solar System, very crudely.'''

    # the data in the table probably don't need to be updated
    expiration = np.inf

    def __init__(self, **kwargs):
        '''
        Initialize a population of Solar System (main) planets.
        '''
        PredefinedPopulation.__init__(self, label='Solar System', **kwargs)
        self.color = 'cornflowerblue'
        self.zorder = 1e10
        self.s = 80
        self.ink = False
        self.exact = True

    def load_raw(self):
        '''
        Load the raw table of data for the Solar System.
        '''

        # load the table of Solar System planets
        raw = ascii.read(initial_filename)

        # for debugging, hang on to the raw table as a hidden attribute
        self._raw = raw

        # a table of unstandardized planet properties
        return raw

    def create_standard(self, trimmed):
        '''
        Create a standardized table of planet properties.
        It must at least contain the columns in
        `attribute_columns`.
        '''

        # start from the trimmed table
        t = trimmed

        # create an empty standardized table
        s = Table()
        s['name'] = t['name']
        s['hostname'] = 'Sun'

        # set up the Sun
        s['stellar_teff'] = 5780*u.K

        s['stellar_radius'] = 1.0*u.Rsun
        s['stellar_mass'] = 1.0*u.Msun
        s['stellar_age'] = 4.5*u.Gyr
        
        # store the period, by default, in day
        s['period'] = (t['period']*u.year).to(u.day)

        # store an eccentricity and longitude of periastron
        s['e'] = np.nan
        s['omega'] = np.nan*u.deg

        # hide the stellar magnitudes
        s['Jmag'] = np.nan
        s['Vmag'] = np.nan

        # pull out the radius and mass
        s['radius'] = (t['radius']*u.km).to(u.Rearth)
        s['radius_uncertainty_lower'] = 0.0*u.Rearth
        s['radius_uncertainty_upper'] = 0.0*u.Rearth
        assert((s['radius'] < 20*u.Rearth).all())
        s['mass'] = (t['mass']*u.kg).to(u.Mearth)
        s['mass_uncertainty_lower'] = 0.0*u.Mearth
        s['mass_uncertainty_upper'] = 0.0*u.Mearth

        # use Kepler's (actual) 3rd Law to get semimajor axis
        semimajor_axis = (s['period'].to(u.year).value)**(2.0/3.0)*u.AU
        s['semimajoraxis'] = semimajor_axis
        s['transit_ar'] = (semimajor_axis/u.Rsun).decompose().value

        # equilibrium temperature assuming uniform redistribution + 0 albedo
        s['teq'] = s['stellar_teff']*(0.25*1.0/s['transit_ar']**2)**0.25

        # some other planet parameters we might not need for the solar system
        s['rv_semiamplitude'] = np.nan*u.m/u.s
        s['radius_ratio'] = (s['radius'].quantity/s['stellar_radius'].quantity).decompose().value
        s['distance'] = np.nan*u.pc
        s['ra'] = 0.0*u.deg
        s['dec'] = 0.0*u.deg
        s['discoverer'] = 'humans'
        s['transit_epoch'] = np.nan*u.day
        s['transit_duration'] = np.nan*u.day
        s['transit_depth'] = (s['radius'].quantity/s['stellar_radius'].quantity).decompose()**2
        s['transit_b'] = 0.0
        s['inclination'] = 90*u.deg

        self.standard = s
        return s
