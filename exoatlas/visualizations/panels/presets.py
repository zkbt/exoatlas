from ...imports import *
from ...models import plot_both_seager
from .BubblePanel import BubblePanel
from .ErrorPanel import ErrorPanel

# FIXME - perhaps define some common x-axes and y-axes, so that we
# can define new classes simply from inheriting from them in the
# correct order. (We'll need to look up and practice multiple inherits.)

class FluxRadius(BubblePanel):
    xsource='relative_insolation'
    xlabel='Bolometric Flux Received (relative to Earth)'
    xscale='log'
    xlim=[6e4, 2e-4]

    ysource='radius'
    ylabel='Planet Radius (Earth radii)'
    yscale='log'
    ylim = [0.3, 30]

    def add_teqaxis(self):
        '''
        Add an extra axis along the bottom of this panel,
        quoting the equilibrium temperature associated
        with a particular bolometric flux received.
        '''

        ax_temp = self.ax.twiny()
        mn, mx = self.ax.get_xlim()
        ax_temp.set_xscale('log')
        teqearth = (5780*u.K/np.sqrt(2*u.AU/u.Rsun)).to(u.K)
        ax_temp.set_xlim(teqearth*mn**.25, teqearth*mx**.25)
        extra_axis_color = 'gray'

        ax_temp.set_xlabel('Planetary Equilibrium Temperature (for zero albedo + efficient circulation)', color=extra_axis_color)
        ax_temp.xaxis.set_ticks_position('bottom') # set the position of the second x-axis to bottom
        ax_temp.xaxis.set_label_position('bottom') # set the position of the second x-axis to bottom
        ax_temp.spines['bottom'].set_position(('outward', 40))
        ax_temp.tick_params(axis='x', colors=extra_axis_color, which='both')
        ax_temp.spines['bottom'].set_edgecolor(extra_axis_color)
        ax_temp.xaxis.set_major_formatter(FormatStrFormatter('%dK'))
        plt.sca(self.ax)

    def plot_hz(self, color='cornflowerblue', alpha=0.25, linewidth=0, **kw):
        '''
        Add a bar that indicates an approximate habitable zone.
        (Estimated very roughly by eye from Kopparapu et al.)
        '''
        plt.axvspan(1.4, 0.4, color=color, alpha=alpha, linewidth=linewidth, **kw)

class FluxTeff(FluxRadius):
    ysource = 'stellar_teff'
    ylabel = 'Stellar Effective Temperature (K)'
    yscale = 'linear'
    ylim = [2000, 7000]


class DistanceRadius(FluxRadius):
    xsource = 'distance'
    xlabel = 'Distance\n(parsecs)'
    xlim = [5,1000]



class EscapeRadius(FluxRadius):
    xsource = 'escape_parameter'
    xlabel = '$\lambda = E_{grav}/E_{thermal}$'
    xlim = [None, None]

class PlanetDensityRadius(FluxRadius):
    xsource = 'density'
    xlabel = 'Planet Density\n(g/cm$^3$)'
    xlim = [0.3, 12]

class StellarRadius(FluxRadius):
    xsource = 'stellar_radius'
    xlabel = 'Stellar Radius\n(solar radii)'
    xscale = 'linear'
    xlim = [0.0, 2.0]

class DepthRadius(FluxRadius):
    xsource = 'transit_depth'
    xlabel = 'Transit Depth\nof Planet'
    xscale = 'log'
    xlim = [2e-6, 2e-1]


class TransmissionRadius(DepthRadius):
    xsource = 'transmission_signal'
    xlabel = 'Transit Depth\nof 1 Scale Height\n of H$_2$-rich Planet'

class ReflectionRadius(DepthRadius):
    xsource = 'reflection_signal'
    xlabel = 'Eclipse Depth\nin Reflected Light\n(100% albedo)'

class EmissionRadius(DepthRadius):
    xsource = 'emission_signal'

    def __init__(self, wavelength=5*u.micron, **kw):
        '''
        Initialize for a particular wavelength, because the
        eclipse depth will depend on the thermal emission
        spectrum of the planet and star.
        '''
        BubblePanel.__init__(self, **kw)
        self.wavelength = wavelength
        self.xlabel = f'Eclipse Depth\nin Thermal Emission\nat $\lambda={self.wavelength.to(u.micron).value}\mu m$'

    @property
    def x(self):
        return self.pop.emission_signal(self.wavelength)


class DistanceBrightness(DistanceRadius):
    ysource='stellar_brightness'
    ylabel='Stellar Brightness at Earth\n(photons/s/m$^2$/nm)'
    yscale='log'
    ylim=[1, 1e5]

    def __init__(self, wavelength=5*u.micron,
                       JWST=False,
                       **kw):
        '''
        Initialize for a particular wavelength.
        '''
        BubblePanel.__init__(self, **kw)
        self.wavelength = wavelength
        w = self.wavelength.to(u.micron).value

        self.JWST = JWST
        if JWST:
            unit = 'Gigaphotons/JWST/hr/R=20'
            self.ylim = [.5e-1, 5e3]
        else:
            unit = 'photons/s/m$^2$/nm'
        self.ylabel = f'Stellar Brightness\nat Earth at $\lambda={w}\mu m$\n({unit})'

    @property
    def y(self):
        if self.JWST:
            return self.pop.stellar_brightness_JWST(self.wavelength)
        else:
            return self.pop.stellar_brightness(self.wavelength)


class DepthBrightness(DistanceBrightness):
    xsource = 'transit_depth'
    xlabel = 'Transit Depth\nof Planet'
    xscale = 'log'
    xlim = [2e-6, 2e-1]

    def plot_sigma(self, color='black', linewidth=3, alpha=0.5, **kw):
        '''
        Plot the 1-sigma uncertainty from photon noise.
        '''

        w = self.wavelength
        photons = np.logspace(0, 14)*u.ph
        sigma = 1/np.sqrt(photons).decompose().value
        # FIXME (we should porbably move the unit handling to the plot panels?)
        photons_in_unit = photons.to(self.pop.photon_unit)/self.pop.JWST_transit_unit(w)
        plt.plot(sigma, photons_in_unit, color=color, linewidth=linewidth, alpha=alpha, **kw)

class TransmissionBrightness(DepthBrightness):
    xsource = 'transmission_signal'
    xlabel = 'Transit Depth\nof 1 Scale Height\n of H$_2$-rich Planet'
    xscale = 'log'

class ReflectionBrightness(DepthBrightness):
    xsource = 'reflection_signal'
    xlabel = 'Eclipse Depth\nin Reflected Light\n(100% albedo)'
    xscale = 'log'

class EmissionBrightness(DepthBrightness):
    xsource = 'emission_signal'
    xscale = 'log'

    def __init__(self, wavelength=5*u.micron,
                       JWST=False,
                       **kw):
        '''
        Initialize for a particular wavelength, because the
        eclipse depth will depend on the thermal emission
        spectrum of the planet and star.
        '''
        BubblePanel.__init__(self, **kw)
        self.wavelength = wavelength
        self.xlabel = f'Eclipse Depth\nin Thermal Emission\nat $\lambda={self.wavelength.to(u.micron).value}\mu m)$'

        w = self.wavelength.to(u.micron).value
        self.JWST = JWST
        if JWST:
            unit = 'Gigaphotons/JWST/hr/R=20'
            self.ylim = [.5e-1, 5e3]
        else:
            unit = 'photons/m$^2$/s/nm'
        self.ylabel = f'Stellar Brightness\nat Earth at $\lambda={w}\mu m$\n({unit})'

    @property
    def x(self):
        return self.pop.emission_signal(self.wavelength)

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
    xlabel = 'Planet Mass\n(Earth masses)'
    ylabel = 'Planet Radius (Earth radii)'
    xscale = 'log'
    yscale = 'log'
    xlim = [0.03, 3000]
    ylim = [0.3, 30]

    @property
    def x(self):
        return self.pop.mass

    @property
    def y(self):
        return self.pop.radius

    @property
    def x_lowerupper(self):
        return self.pop.uncertainty_lowerupper('mass')

    @property
    def y_lowerupper(self):
        return self.pop.uncertainty_lowerupper('radius')

    def plot_both_seager(self, **kw):
        plt.sca(self.ax)
        plot_both_seager(**kw)

    def plot_constant_density(self, densities=10.0**np.arange(-4,7)*u.g/u.cm**3,
                                    color='coral',
                                    alpha=0.5,
                                    zorder=-100,
                                    **kw):
        for density in densities:
            mass = np.logspace(-2, 3)*u.Mearth
            radius = (mass/density/4/np.pi*3)**(1/3)
            plt.plot(mass, radius.to(u.Rearth),
                             color=color,
                             alpha=alpha,
                             zorder=zorder,
                             **kw)

class FluxEscape(ErrorPanel):

    ysource='escape_velocity'
    ylabel='Escape\nVelocity\n(km/s)'
    yscale = 'log'
    ylim=[2, 500]

    xsource='relative_insolation'
    xlabel='Bolometric Flux Received (relative to Earth)'
    xscale='log'
    xlim=[6e4, 2e-4]

    @property
    def x_lowerupper(self):
        sigma = self.x*0.0
        return sigma, sigma

    @property
    def y_lowerupper(self):
        '''
        Calculate a rough estimate of the uncertainty on the escape
        velocity, simply by propagating the mass and radius uncertainties.
        This could probably be done more cleverly by deriving it from
        the planetary surface gravity uncertainty, which as fewer hidden
        covariances because it can be derived directly from transit
        and RV observables.
        '''

        dlnm = self.pop.uncertainty('mass')/self.pop.mass
        dlnr = self.pop.uncertainty('radius')/self.pop.radius


        dlne = np.sqrt((dlnm/2)**2 + (dlnr/2)**2)
        sigma = dlne*self.pop.escape_velocity

        return sigma, sigma

    def plot_constant_lambda(self, alpha=0.5, color='gray', **kw):
        '''
        Plot the escape velocity vs insolation for
        different constant values of the escape parameter.

        This assumes the exosphere is at the planet's equilibrium
        temperature, which is is *terrible* approximation. This
        also represented the *current* escape parameter, saying
        nothing about the history of XUV radiation the planet
        may have received.

        Parameters
        ----------
        **kw : dict
            Keyword parameters get passed along to plot.
        '''

        # label one of the lines
        plt.text(.005, 145, r'$E_{grav}/E_{thermal}$', rotation=-4.5, fontsize=8, color=color)

        teq = np.logspace(1, 4)*u.K
        earth_teq = (5780*u.K/np.sqrt(2*u.AU/u.Rsun)).to(u.K)

        relative_insolation = (teq/earth_teq)**4
        m = 1*u.M_p

        def escape_velocity(T, lam=1):
            # calculate the escape velocity for a particular escape parameter
            kT = con.k_B*T
            return np.sqrt(2*kT*lam/m).to('km/s')

        plt.sca(self.ax)
        max_insolation = self.xlim[1]
        max_teq = earth_teq*max_insolation**.25

        # loop over factors of 10 of lambda
        for lam in 10**np.arange(6):

            plt.plot(relative_insolation,
                     escape_velocity(teq, lam=lam),
                     color=color, alpha=alpha, **kw)


class SemimajorRadius(FluxRadius):
    xsource = 'semimajoraxis'
    xlabel = 'Semimajor Axis (AU)\n'
    xscale = 'log'
    xlim = [0.001, 1000]

class SemimajorMass(SemimajorRadius):
    ysource = 'kludge_mass'
    ylabel = 'Planet Mass\n(Earth masses)'
    yscale = 'log'
    ylim = [0.03, 4000]
