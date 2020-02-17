from ...imports import *
from .plottable import *

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
    lim = [0.3, 12]

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

    def __init__(self, wavelength=5*u.micron,
                       JWST=False,
                       **kw):
        '''
        Initialize for a particular wavelength.
        '''
        PlottableAxis.__init__(self, **kw)
        self.wavelength = wavelength
        w = self.wavelength.to(u.micron).value

        self.JWST = JWST
        if JWST:
            unit = 'Gigaphotons/JWST/hr/R=20'
            self.lim = [.5e-1, 5e3]
        else:
            unit = 'photons/s/m$^2$/nm'
        self.label = f'Stellar Brightness\nat Earth at $\lambda={w}\mu m$\n({unit})'

    def value(self):
        if self.JWST:
            return self.panel.pop.stellar_brightness_JWST(self.wavelength)
        else:
            return self.panel.pop.stellar_brightness(self.wavelength)
