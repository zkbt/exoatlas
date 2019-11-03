from ..imports import *
from .panels import BubblePanel, ErrorPanel

import craftroom.units as u

class DistanceRadius(BubblePanel):
    title = ''
    xlabel = 'Distance\n(parsecs)'
    ylabel = 'Planet Radius (Earth radii)'
    xscale = 'log'
    yscale = 'log'
    xlim = [3,1000]
    ylim = [0.4, 40]

    def __init__(self, lightyears=False, **kw):
        BubblePanel.__init__(self, **kw)

        self.lightyears = lightyears
        if lightyears:
            self.xlabel = "Distance from Earth (lightyears)"
            self.xlim = np.array([3,1000])*3.26

        self.normalization = 0.0001#self.unnormalizedsize()[self.pop.standard['name'] == 'GJ 1214b']

    def unnormalizedsize(self):
        # set the symbol size to be transit depth
        #return (self.pop.planet_radius/self.pop.stellar_radius)**2
        return self.pop.transit_depth

    def get_sizes(self):
        return self.unnormalizedsize()/self.normalization

    @property
    def x(self):
        if self.lightyears:
            return self.pop.stellar_distance.to(u.ly)
        else:
            return self.pop.stellar_distance

    @property
    def y(self):
        return self.pop.planet_radius

    def build(self, **kw):
        for key in self.pops.keys():
            self.plot(key, **kw)

class EscapeRadius(DistanceRadius):
    xlabel = '$\lambda = E_{grav}/E_{thermal}$'
    xscale = 'log'
    xlim = [None, None]

    @property
    def x(self):
        return self.pop.escape_parameter

class FluxRadius(DistanceRadius):
    xlabel = 'Flux Received\n(relative to Earth)'
    xscale = 'log'
    xlim = [5e4, 0.5]

    @property
    def x(self):
        return (self.pop.insolation/self.pop.earth_insolation)
Teq = FluxRadius

class PlanetDensityRadius(DistanceRadius):

    title = ''
    xlabel = 'Planet Density\n(g/cm$^3$)'
    xscale = 'log'
    xlim = [0.3, 12]

    @property
    def x(self):
        mass = self.pop.planet_mass
        volume = 4/3*np.pi*(self.pop.planet_radius)**3
        return (mass/volume).to('g/cm**3')

class StellarRadius(DistanceRadius):
    xlabel = 'Stellar Radius\n(solar radii)'
    xscale = 'linear'
    xlim = [0.1, 1.1]

    @property
    def x(self):
        return self.pop.stellar_radius

class JRadius(DistanceRadius):
    xlabel = 'J (magnitude)\n'
    xscale = 'linear'
    xlim = [3.5, 14.5]

    @property
    def x(self):
        return self.pop.J

class PeriodRadius(DistanceRadius):
    xlabel = 'Period (days)\n'
    xscale = 'log'
    xlim = [0.15, 365]

    @property
    def x(self):
        return self.pop.period

class MassRadius(ErrorPanel):
    
    title = ''
    xlabel = 'Planet Mass (Earth masses)'
    ylabel = 'Planet Radius (Earth radii)'
    xscale = 'log'
    yscale = 'log'
    xlim = [0.03, 3000]
    ylim = [0.3, 25]

    @property
    def x(self):
        return self.pop.planet_mass

    @property
    def y(self):
        return self.pop.planet_radius

    @property
    def x_lowerupper(self):
        return self.pop.uncertainty_lowerupper('planet_mass')

    @property
    def y_lowerupper(self):
        return self.pop.uncertainty_lowerupper('planet_radius')
