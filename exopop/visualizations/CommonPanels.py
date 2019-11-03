from ..imports import *
from .panels import BubblePanel, ErrorPanel
from ..models import plot_both_seager


class Scatter(BubblePanel):
    def __init__(self, **kw):

        BubblePanel.__init__(self, **kw)

        for k in ['xsource', 'ysource',
                  'xlabel', 'ylabel',
                  'xscale', 'yscale',
                  'xlim', 'ylim']:
            if k in kw:
                vars(self)[k] = kw[k]

    @property
    def x(self):
        return getattr(self.pop, self.xsource)

    @property
    def y(self):
        return getattr(self.pop, self.ysource)

class FluxRadius(Scatter):
    xsource='relative_insolation'
    ysource='planet_radius'
    xlabel='Flux Received\n(relative to Earth)'
    ylabel='Planet Radius (Earth radii)'
    xscale='log'
    yscale='log'
    xlim=[5e4,5e-4]
    ylim=[.25, 25]

class FluxEscape(FluxRadius):
    ysource='escape_velocity'
    ylabel='Escape Velocity (m/s)'
    ylim=[None, None]

class DistanceRadius(FluxRadius):
    xsource = 'stellar_distance'
    xlabel = 'Distance\n(parsecs)'
    xlim = [3,1000]

    '''        if lightyears:
                self.xlabel = "Distance from Earth (lightyears)"
                self.xlim = np.array([3,1000])*3.26
    '''
class EscapeRadius(FluxRadius):
    xsource = 'escape_parameter'
    xlabel = '$\lambda = E_{grav}/E_{thermal}$'
    xlim = [None, None]

class PlanetDensityRadius(FluxRadius):
    xsource = 'planet_density'
    xlabel = 'Planet Density\n(g/cm$^3$)'
    xlim = [0.3, 12]

class StellarRadius(FluxRadius):
    xsource = 'stellar_radius'
    xlabel = 'Stellar Radius\n(solar radii)'
    xscale = 'linear'
    xlim = [0.08, 1.5]

class JRadius(FluxRadius):
    xsource = 'Jmag'
    xlabel = 'J (magnitude)\n'
    xscale = 'linear'
    xlim = [3.5, 14.5]

class PeriodRadius(FluxRadius):
    xsource = 'period'
    xlabel = 'Period (days)\n'
    xscale = 'log'
    xlim = [0.15, 365]

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

    def plot_both_seager(self):
        plt.sca(self.ax)
        plot_both_seager()
