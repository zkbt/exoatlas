from ...imports import *
from ...models import plot_both_seager
from ..axes.preset_plottables import *
from .Panel import *
from .BubblePanel import *
from .ErrorPanel import *


class FluxRadius(BubblePanel):
    xaxis = Flux
    yaxis = Radius

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

class FluxTeff(BubblePanel):
    xaxis = Flux
    yaxis = StellarTeff

class DistanceRadius(BubblePanel):
    xaxis = Distance
    yaxis = Radius

class EscapeRadius(BubblePanel):
    xaxis = Escape
    yaxis = Radius

class DensityRadius(BubblePanel):
    xaxis = Density
    yaxis = Radius

class StellarRadiusPlanetRadius(BubblePanel):
    xaxis = StellarRadius
    yaxis = Radius

class DepthRadius(BubblePanel):
    xaxis = Depth
    yaxis = Radius

class TransmissionRadius(BubblePanel):
    xaxis = Transmission
    yaxis = Radius

class ReflectionRadius(BubblePanel):
    xaxis = Reflection
    yaxis = Radius

class EmissionRadius(BubblePanel):
    xaxis = Emission
    yaxis = Radius

class DistanceBrightness(BubblePanel):
    xaxis = Distance
    yaxis = StellarBrightness

class DepthBrightness(BubblePanel):
    xaxis = Depth
    yaxis = StellarBrightness

    def plot_sigma(self, color='black', linewidth=3, alpha=0.5, **kw):
        '''
        Plot the 1-sigma uncertainty from photon noise.
        '''

        w = self.plottable['y'].wavelength
        photons = np.logspace(0, 14)*u.ph
        sigma = 1/np.sqrt(photons).decompose().value

        # what is "1" (e.g. JWST for one hour at R=20)?
        telescope_unit = self.plottable['y'].telescope_unit

        # how many photons do we collect with that one?
        photons_collected = photons/telescope_unit

        unit = self.plottable['y'].unit
        photons_in_unit = photons_collected.to(unit)


        plt.plot(sigma, photons_in_unit, color=color, linewidth=linewidth, alpha=alpha, **kw)

class TransmissionBrightness(DepthBrightness):
    xaxis = Transmission

class ReflectionBrightness(DepthBrightness):
    xaxis = Reflection

class EmissionBrightness(DepthBrightness):
    xaxis = Emission

class JRadius(BubblePanel):
    xaxis = Jmag
    yaxis = Radius

class PeriodRadius(BubblePanel):
    xaxis = Period
    yaxis = Radius

class SemimajorRadius(BubblePanel):
    xaxis = SemimajorAxis
    yaxis = Radius

class SemimajorMass(SemimajorRadius):
    xaxis = SemimajorAxis
    yaxis = KludgedMass

class MassRadius(ErrorPanel):
    xaxis = Mass
    yaxis = Radius

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
    xaxis = Flux
    yaxis = EscapeVelocity

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


# construct a list of predefined panels
predefined_panels = []
options = list(locals().items())
for k, v in options:
    try:
        assert(issubclass(v, Panel))
        assert(v != BubblePanel)
        assert(v != ErrorPanel)
        assert(v != Panel)
        predefined_panels.append(v)
    except (AssertionError, TypeError):
        pass
