from ..imports import *
from .panels import BubblePanel, ErrorPanel
from ..models import plot_both_seager


"""class Scatter(BubblePanel):
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
"""
Scatter = BubblePanel

class FluxRadius(BubblePanel):
    xsource='relative_insolation'
    xlabel='Bolometric Flux Received (relative to Earth)'
    xscale='log'
    xlim=[6e4, 2e-4]

    ysource='planet_radius'
    ylabel='Planet Radius (Earth radii)'
    yscale='log'
    ylim = [0.3, 30]

    def add_teq_axis(self):

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
        plt.axvspan(1.4, 0.4, color=color, alpha=alpha, linewidth=linewidth, **kw)


class DistanceRadius(FluxRadius):
    xsource = 'stellar_distance'
    xlabel = 'Distance\n(parsecs)'
    xlim = [5,1000]



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
    xlim = [0.0, 2.0]

class DepthRadius(FluxRadius):
    xsource = 'transit_depth'
    xlabel = 'Transit Depth\nof Planet'
    xscale = 'log'
    xlim = [2e-6, 2e-1]


class TransmissionRadius(DepthRadius):
    xsource = 'transmissionsignal'
    xlabel = 'Transit Depth\nof 1 Scale Height\n of H$_2$-rich Planet'

class ReflectionRadius(DepthRadius):
    xsource = 'reflectionsignal'
    xlabel = 'Eclipse Depth\nin Reflected Light\n(100% albedo)'

class EmissionRadius(DepthRadius):
    xsource = 'emissionsignal'

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
        return self.pop.emissionsignal(self.wavelength)


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
    xsource = 'transmissionsignal'
    xlabel = 'Transit Depth\nof 1 Scale Height\n of H$_2$-rich Planet'
    xscale = 'log'

class ReflectionBrightness(DepthBrightness):
    xsource = 'reflectionsignal'
    xlabel = 'Eclipse Depth\nin Reflected Light\n(100% albedo)'
    xscale = 'log'

class EmissionBrightness(DepthBrightness):
    xsource = 'emissionsignal'
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
        return self.pop.emissionsignal(self.wavelength)

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

    def plot_both_seager(self, **kw):
        plt.sca(self.ax)
        plot_both_seager(**kw)

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

        dlnm = self.pop.uncertainty('planet_mass')/self.pop.planet_mass
        dlnr = self.pop.uncertainty('planet_radius')/self.pop.planet_radius


        dlne = np.sqrt((dlnm/2)**2 + (dlnr/2)**2)
        sigma = dlne*self.pop.escape_velocity

        return sigma, sigma

    def plot_constant_lambda(self, alpha=0.5, color='gray', **kw):
        '''
        Plot the escape velocity vs insolation for
        different constant values of the escape parameter.

        Parameters
        ----------
        **kw : dict
            Keyword parameters get passed along to plot.
        '''

        # label one of the lines
        plt.text(.005, 145, r'$E_{grav}/E_{thermal}$', rotation=-4.5, fontsize=8, color=color)

        teq = np.logspace(1, 4)*u.K
        teq_earth = (5780*u.K/np.sqrt(2*u.AU/u.Rsun)).to(u.K)

        relative_insolation = (teq/teq_earth)**4
        m = 1*u.M_p

        def escape_velocity(T, lam=1):
            # calculate the escape velocity for a particular escape parameter
            kT = con.k_B*T
            return np.sqrt(2*kT*lam/m).to('km/s')

        plt.sca(self.ax)
        max_insolation = self.xlim[1]
        max_teq = teq_earth*max_insolation**.25

        # loop over factors of 10 of lambda
        for lam in 10**np.arange(6):

            plt.plot(relative_insolation,
                     escape_velocity(teq, lam=lam),
                     color=color, alpha=alpha, **kw)
