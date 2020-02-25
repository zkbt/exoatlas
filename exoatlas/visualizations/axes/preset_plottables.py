from ...imports import *
from .plottable import *
from ...telescopes import *

class Flux(PlottableAxis):
    source = 'relative_insolation'
    label = 'Bolometric Flux Received (relative to Earth)'
    scale = 'log'
    lim = [6e4, 2e-4]

class Radius(PlottableAxis):
    source = 'radius'
    label = 'Planet Radius (Earth radii)'
    scale = 'log'
    lim = [0.3, 30]

    def value_lowerupper(self):
        return self.panel.pop.uncertainty_lowerupper('radius')

class Mass(PlottableAxis):
    source = 'mass'
    label = 'Planet Mass\n(Earth masses)'
    scale = 'log'
    lim = [0.03, 3000]

    def value_lowerupper(self):
        return self.panel.pop.uncertainty_lowerupper('mass')

class SemimajorAxis(PlottableAxis):
    source = 'semimajoraxis'
    label = 'Semimajor Axis (AU)\n'
    scale = 'log'
    lim = [0.001, 1000]

class KludgedMass(PlottableAxis):
    source = 'kludge_mass'
    label = 'Planet Mass or msini\n(Earth masses)'
    scale = 'log'
    lim = [0.03, 4000]

class StellarTeff(PlottableAxis):
    source = 'stellar_teff'
    label = 'Stellar Effective Temperature (K)'
    scale = 'linear'
    lim = [2000, 7000]

class Distance(PlottableAxis):
    source = 'distance'
    label = 'Distance\n(parsecs)'
    scale = 'log'
    lim = [5,1000]

class EscapeVelocity(PlottableAxis):
    source='escape_velocity'
    label='Escape\nVelocity\n(km/s)'
    scale = 'log'
    lim=[2, 500]

class Escape(PlottableAxis):
    source = 'escape_parameter'
    label = '$\lambda = E_{grav}/E_{thermal}$'
    scale = 'log'
    lim = [None, None]

    def value_lowerupper(self):
        '''
        Calculate a rough estimate of the uncertainty on the escape
        velocity, simply by propagating the mass and radius uncertainties.
        This could probably be done more cleverly by deriving it from
        the planetary surface gravity uncertainty, which as fewer hidden
        covariances because it can be derived directly from transit
        and RV observables.
        '''

        m = self.panel.pop.mass
        sigma_m = self.panel.pop.uncertainty('mass')

        r = self.panel.pop.radius
        sigma_r = self.panel.pop.uncertainty('radius')

        dlnm = sigma_m/m
        dlnr = sigma_r/r

        dlne = np.sqrt((dlnm/2)**2 + (dlnr/2)**2)
        sigma = dlne*self.panel.pop.escape_velocity

        return sigma, sigma


class Density(PlottableAxis):
    source = 'density'
    label = 'Planet Density\n(g/cm$^3$)'
    scale = 'log'
    lim = [0.01, 100]

class StellarRadius(PlottableAxis):
    source = 'stellar_radius'
    label = 'Stellar Radius\n(solar radii)'
    scale = 'linear'
    lim = [0.0, 2.0]


class Period(PlottableAxis):
    source = 'period'
    label = 'Period (days)\n'
    scale = 'log'
    lim = [0.15, 365]

class Jmag(PlottableAxis):
    source = 'Jmag'
    label = 'J (magnitude)\n'
    scale = 'linear'
    lim = [3.5, 14.5]

class Depth(PlottableAxis):
    source = 'transit_depth'
    label = 'Transit Depth\nof Planet'
    scale = 'log'
    lim = [2e-6, 2e-1]

class Transmission(Depth):
    source = 'transmission_signal'
    label = 'Transit Depth\nof 1 Scale Height\n of H$_2$-rich Planet'

class Reflection(Depth):
    source = 'reflection_signal'
    label = 'Eclipse Depth\nin Reflected Light\n(100% albedo)'

class Emission(Depth):
    source = 'emission_signal'

    def __init__(self, wavelength=5*u.micron, **kw):
        '''
        Initialize for a particular wavelength, because the
        eclipse depth will depend on the thermal emission
        spectrum of the planet and star.
        '''
        PlottableAxis.__init__(self, **kw)
        self.wavelength = wavelength
        self.label = f'Eclipse Depth\nin Thermal Emission\nat $\lambda={self.wavelength.to(u.micron).value}\mu m$'

    def value(self):
        return self.panel.pop.emission_signal(self.wavelength)

class StellarBrightness(PlottableAxis):
    source='stellar_brightness'
    label='Stellar Brightness at Earth\n(photons/s/m$^2$/nm)'
    scale='log'
    lim=[1, 1e5]

    def __init__(self, panel=None,
                       orientation=None,
                       wavelength=None,
                       telescope=None,
                       R=20,
                       dt=1*u.hr,
                       **kw):
        '''
        Initialize the StellarBrightness plottable.
        It depends on wavelength, due to the thermal
        emission Planck spectrum.

        Parameters
        ----------
        wavelength : astropy.units.quantity.Quantity
            The wavelength at which we want the
            stellar brightness to be computed.

        telescope : None, str
            The telescope unit in which to express the
            stellar brightness. Options include:
                'JWST'
                'Hubble'

        R : float
            The spectral resolution at which the
            telescope will be binned.
            (Ignored if telescope is None.)

        dt : astropy.units.quantity.Quantity
            The time over which the telescope exposes.
            (Ignored if telescope is None.)
        '''

        # initialize the basic plottable axis
        PlottableAxis.__init__(self, panel=panel, orientation=orientation, **kw)

        # keep track of the telescope (if any)
        self.telescope = telescope

        # select the appropriate photon unit
        if self.telescope is None:
            if wavelength is None:
                wavelength = 1.0*u.micron
            self.wavelength = wavelength
            self.telescope_unit = u.s*u.m**2*u.micron
            self.unit = u.Unit(('ph s^-1 m^-2 micron^-1'))
            unit_string = 'photons/s/m$^2$/$\mu$m'
            self.lim = [1e2, 1e8]
        else:
            self.telescope_unit = define_telescope_unit_by_name(self.telescope,
                                wavelength=wavelength, R=R, dt=dt, **kw)
            self.unit = photon_unit/self.telescope_unit
            self.wavelength = self.telescope_unit.wavelength
            # if self.telescope == 'JWST':
            #     self.unit = photon_unit/define_JWST_unit(wavelength=wavelength,
            #                                              R=R, dt=dt)
            # if self.telescope == 'HST':
            #     self.unit = photon_unit/define_HST_unit(wavelength=wavelength,
            #                                             R=R, dt=dt)
            # unit_string = f'{photon_unit}/{self.telescope}/{dt}/R={R}'
            unit_string = self.unit.to_string()
            self.lim = [1e-3, 1e3]


        w = self.wavelength.to(u.micron).value
        self.label = f'Stellar Brightness\nat Earth at $\lambda={w}\mu$m\n({unit_string})'

    def value(self):
        return self.panel.pop.stellar_brightness(self.wavelength).to(self.unit)
