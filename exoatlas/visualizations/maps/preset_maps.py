from ...imports import *
from ...models import plot_both_seager, make_hz
from ..plottables.preset_plottables import *
from .Map import *
from .BubbleMap import *
from .ErrorMap import *


class Flux_x_Radius(ErrorMap):
    xaxis = Flux
    yaxis = Radius

    def add_teqaxis(self, position=50):
        """
        Add an extra axis along the bottom of this map,
        quoting the equilibrium temperature associated
        with a particular bolometric flux received.
        """

        ax_temp = self.ax.twiny()
        mn, mx = self.ax.get_xlim()
        ax_temp.set_xscale("log")
        teqearth = (5780 * u.K / np.sqrt(2 * u.AU / u.Rsun)).to(u.K)
        ax_temp.set_xlim(teqearth * mn**0.25, teqearth * mx**0.25)
        extra_axis_color = "gray"

        ax_temp.set_xlabel(
            "Planetary Equilibrium Temperature (for zero albedo + efficient circulation)",
            color=extra_axis_color,
        )
        ax_temp.xaxis.set_ticks_position(
            "bottom"
        )  # set the position of the second x-axis to bottom
        ax_temp.xaxis.set_label_position(
            "bottom"
        )  # set the position of the second x-axis to bottom
        ax_temp.spines["bottom"].set_position(("outward", position))
        ax_temp.tick_params(axis="x", colors=extra_axis_color, which="both")
        ax_temp.spines["bottom"].set_edgecolor(extra_axis_color)
        ax_temp.xaxis.set_major_formatter(FormatStrFormatter("%dK"))
        plt.sca(self.ax)

    def plot_hz(
        self,
        Teff=[2600, 6000],
        inner="moist-greenhouse",
        outer="maximum-greenhouse",
        color="cornflowerblue",
        alpha=0.25,
        linewidth=0,
        **kw
    ):
        """
        Add a bar that indicates an approximate habitable zone.
        (Estimated very roughly by eye from Kopparapu et al.)
        """

        # make sure the plotting happens in this map
        plt.sca(self.ax)

        # define the functions
        S_inner = make_hz(inner)
        S_outer = make_hz(outer)

        if np.atleast_1d(Teff).size == 1:
            plt.axvspan(
                S_inner(Teff),
                S_outer(Teff),
                color=color,
                alpha=alpha,
                linewidth=linewidth,
                **kw
            )
        elif np.atleast_1d(Teff).size == 2:
            N = 10
            for Teff in np.linspace(Teff[0], Teff[1], N):
                plt.axvspan(
                    S_inner(Teff),
                    S_outer(Teff),
                    color=color,
                    alpha=alpha / N,
                    linewidth=0,
                    **kw
                )

        # Teff_inner = 5780
        # Teff_outer = 2600


class Flux_x_Teff(BubbleMap):
    xaxis = Flux
    yaxis = StellarTeff

    def plot_hz(
        self,
        inner="moist-greenhouse",
        outer="maximum-greenhouse",
        color="cornflowerblue",
        alpha=0.25,
        linewidth=0,
        zorder=-100,
        **kw
    ):
        """
        Add a bar that indicates an approximate habitable zone.
        (Estimated very roughly by eye from Kopparapu et al.)
        """

        # make sure the plotting happens in this map
        plt.sca(self.ax)

        # define the functions
        S_inner = make_hz(inner)
        S_outer = make_hz(outer)

        # define the temperature grid
        Teff = np.linspace(2000, 7500)

        # plot the swath
        plt.fill_betweenx(
            Teff,
            S_inner(Teff),
            S_outer(Teff),
            color=color,
            alpha=alpha,
            linewidth=linewidth,
            zorder=zorder,
            **kw
        )


class SemimajorAxis_x_StellarLuminosity(BubbleMap):
    xaxis = SemimajorAxis
    yaxis = StellarLuminosity

    def plot_hz(
        self,
        inner="moist-greenhouse",
        outer="maximum-greenhouse",
        color="cornflowerblue",
        alpha=0.25,
        linewidth=0,
        zorder=-100,
        **kw
    ):
        """
        Add a bar that indicates an approximate habitable zone.
        (Estimated very roughly by eye from Kopparapu et al.)
        """

        # make sure the plotting happens in this map
        plt.sca(self.ax)

        # define the functions
        S_inner = make_hz(inner)
        S_outer = make_hz(outer)

        # define the temperature grid
        Teff = np.linspace(2000, 7500)

        # plot the swath
        plt.fill_betweenx(
            Teff,
            S_inner(Teff),
            S_outer(Teff),
            color=color,
            alpha=alpha,
            linewidth=linewidth,
            zorder=zorder,
            **kw
        )


class Distance_x_Radius(BubbleMap):
    xaxis = Distance
    yaxis = Radius


class Distance_x_Teff(BubbleMap):
    xaxis = Distance
    yaxis = StellarTeff


class EscapeParameter_x_Radius(BubbleMap):
    xaxis = EscapeParameter
    yaxis = Radius


class Density_x_Radius(BubbleMap):
    xaxis = Density
    yaxis = Radius


class StellarRadius_x_PlanetRadius(ErrorMap):
    xaxis = StellarRadius
    yaxis = Radius


class Depth_x_Radius(BubbleMap):
    xaxis = Depth
    yaxis = Radius


class Transmission_x_Radius(BubbleMap):
    xaxis = Transmission
    yaxis = Radius


class Reflection_x_Radius(BubbleMap):
    xaxis = Reflection
    yaxis = Radius


class Emission_x_Radius(BubbleMap):
    xaxis = Emission
    yaxis = Radius


class Distance_x_Brightness(BubbleMap):
    xaxis = Distance
    yaxis = StellarBrightnessTelescope


class Depth_x_Brightness(BubbleMap):
    xaxis = Depth
    yaxis = StellarBrightnessTelescope

    def plot_sigma(self, color="black", linewidth=3, alpha=0.5, **kw):
        """
        Plot the 1-sigma uncertainty from photon noise,
        with sigma on the x axis and number of photons on y-axis,
        in the telescope units associated with the map.

        Parameters
        ----------
        All keywords will be passed to plt.plot
        """

        w = self.plottable["y"].wavelength
        photons = np.logspace(0, 14) * u.ph
        sigma = 1 / np.sqrt(photons).decompose().value

        # what is "1" (e.g. JWST for one hour at R=20)?
        telescope_unit = self.plottable["y"].telescope_unit

        # how many photons do we collect with that one?
        photons_collected = photons / telescope_unit

        unit = self.plottable["y"].unit
        photons_in_unit = photons_collected.to(unit)

        plt.plot(
            sigma, photons_in_unit, color=color, linewidth=linewidth, alpha=alpha, **kw
        )


class Transmission_x_Brightness(Depth_x_Brightness):
    xaxis = Transmission


class Reflection_x_Brightness(Depth_x_Brightness):
    xaxis = Reflection


class Emission_x_Brightness(Depth_x_Brightness):
    xaxis = Emission


class Period_x_Radius(BubbleMap):
    xaxis = Period
    yaxis = Radius


class SemimajorAxis_x_Radius(BubbleMap):
    xaxis = SemimajorAxis
    yaxis = Radius


class SemimajorAxis_x_Mass(SemimajorAxis_x_Radius):
    xaxis = SemimajorAxis
    yaxis = KludgedMass


class Mass_x_Radius(ErrorMap):
    xaxis = Mass
    yaxis = Radius

    def plot_both_seager(self, **kw):
        plt.sca(self.ax)
        plot_both_seager(**kw)

    def plot_constant_density(
        self,
        densities=10.0 ** np.arange(-4, 7) * u.g / u.cm**3,
        color="coral",
        alpha=0.5,
        zorder=-100,
        **kw
    ):
        for density in densities:
            mass = np.logspace(-2, 3) * u.Mearth
            radius = (mass / density / 4 / np.pi * 3) ** (1 / 3)
            plt.plot(
                mass, radius.to(u.Rearth), color=color, alpha=alpha, zorder=zorder, **kw
            )


class Mass_x_EscapeVelocity(ErrorMap):
    xaxis = Mass
    yaxis = EscapeVelocity


class Flux_x_EscapeVelocity(ErrorMap):
    xaxis = Flux
    yaxis = EscapeVelocity

    def plot_constant_lambda(
        self, alpha=0.5, color="gray", x=0.01, y=100, rotation=-4.5, zorder=-100, **kw
    ):
        """
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
        """

        # label one of the lines
        plt.text(
            x,
            y,
            r"$E_{grav}/E_{thermal}$",
            rotation=rotation,
            fontsize=8,
            color=color,
        )

        teq = np.logspace(1, 4) * u.K
        earth_teq = (5780 * u.K / np.sqrt(2 * u.AU / u.Rsun)).to(u.K)

        relative_insolation = (teq / earth_teq) ** 4
        m = 1 * u.M_p

        def escape_velocity(T, lam=1):
            # calculate the escape velocity for a particular escape parameter
            kT = con.k_B * T
            return np.sqrt(2 * kT * lam / m).to("km/s")

        plt.sca(self.ax)
        max_insolation = self.xlim[1]
        max_teq = earth_teq * max_insolation**0.25

        # loop over factors of 10 of lambda
        for lam in 10 ** np.arange(6):
            plt.plot(
                relative_insolation,
                escape_velocity(teq, lam=lam),
                color=color,
                alpha=alpha,
                zorder=zorder,
                **kw
            )

    def plot_jeans_shoreline(
        self,
        reference_relative_insolation=0.43337238,
        reference_escape_velocity=5.02270418 * u.km / u.s,
        alpha=0.5,
        color="gray",
        x=0.01,
        y=100,
        **kw
    ):
        relative_insolation = np.logspace(-5, 5)
        # v ~ T**(1/2) ~ I**(1/8)
        escape_velocity = (
            reference_escape_velocity
            * (relative_insolation / reference_relative_insolation) ** 0.125
        )
        plt.plot(relative_insolation, escape_velocity, alpha=alpha, color=color)

        plt.text(
            reference_relative_insolation,
            reference_escape_velocity,
            r"$v_{thermal} \propto v_{escape}$    ",
            fontsize=8,
            color=color,
            alpha=alpha,
            ha="right",
            **kw
        )

    def plot_shoreline(
        self,
        reference_relative_insolation=0.43337238,
        reference_escape_velocity=5.02270418 * u.km / u.s,
        alpha=0.5,
        color="blue",
        x=0.01,
        y=100,
        **kw
    ):
        relative_insolation = np.logspace(-5, 5)
        escape_velocity = (
            reference_escape_velocity
            * (relative_insolation / reference_relative_insolation) ** 0.25
        )
        plt.plot(relative_insolation, escape_velocity, alpha=alpha, color=color)

        plt.text(
            reference_relative_insolation,
            reference_escape_velocity,
            r"   $I \propto v_{escape}^4$",
            fontsize=8,
            color=color,
            alpha=alpha,
            **kw
        )


class EscapeVelocity_x_Flux(BubbleMap):
    xaxis = EscapeVelocity
    yaxis = Flux

    def plot_jeans_shoreline(
        self,
        reference_relative_insolation=0.43337238,
        reference_escape_velocity=5.02270418 * u.km / u.s,
        alpha=0.5,
        color="gray",
        x=0.01,
        y=100,
        **kw
    ):
        relative_insolation = np.logspace(-5, 5)
        # v ~ T**(1/2) ~ I**(1/8)
        escape_velocity = (
            reference_escape_velocity
            * (relative_insolation / reference_relative_insolation) ** 0.125
        )
        plt.plot(
            escape_velocity,
            relative_insolation,
            alpha=alpha,
            color=color,
            label=r"$v_{thermal} \propto v_{escape}$",
        )

    def plot_shoreline(
        self,
        reference_relative_insolation=0.43337238,
        reference_escape_velocity=5.02270418 * u.km / u.s,
        alpha=0.5,
        color="blue",
        x=0.01,
        y=100,
        **kw
    ):
        relative_insolation = np.logspace(-5, 5)
        escape_velocity = (
            reference_escape_velocity
            * (relative_insolation / reference_relative_insolation) ** 0.25
        )
        plt.plot(
            escape_velocity,
            relative_insolation,
            alpha=alpha,
            color=color,
            label=r"$I \propto v_{escape}^4$",
        )


class EscapeVelocity_x_CumulativeXUVFlux(EscapeVelocity_x_Flux):
    yaxis = CumulativeXUVFlux


class CumulativeXUVFlux_x_EscapeVelocity(Flux_x_EscapeVelocity):
    xaxis = CumulativeXUVFlux


class ImpactVelocity_x_EscapeVelocity(Flux_x_EscapeVelocity):
    xaxis = ImpactVelocity

    def plot_shoreline(
        self,
        reference_impact_velocity=24.650357 * u.km / u.s,
        reference_escape_velocity=5.02270418 * u.km / u.s,
        alpha=0.5,
        color="brown",
        x=0.01,
        y=100,
        rotation=-4.5,
        **kw
    ):
        impact_velocity = np.logspace(0, 3) * u.km / u.s
        escape_velocity = (
            reference_escape_velocity
            * (impact_velocity / reference_impact_velocity) ** 0.25
        )
        plt.plot(impact_velocity, escape_velocity, alpha=alpha, color=color)

        plt.text(
            reference_impact_velocity,
            reference_escape_velocity,
            r"   $v_{impact} \propto {v_{escape}}$",
            fontsize=8,
            color=color,
            alpha=alpha,
            va="center",
            ha="left",
        )


class RA_x_Dec(BubbleMap):
    xaxis = RightAscension
    yaxis = Declination

    def add_monthaxis(self, ax=None, position=40):
        """
        Add an extra axis along the bottom of this map,
        quoting the month when things are best to observe
        """
        ax_month = ax.twiny()
        ax_month.set_xticks([24, 18, 12, 6, 0], ["Sep", "Jun", "Mar", "Dec", "Sep"])
        extra_axis_color = "gray"
        ax_month.xaxis.set_ticks_position("bottom")
        ax_month.xaxis.set_label_position("bottom")
        ax_month.spines["bottom"].set_position(("outward", position))
        ax_month.tick_params(axis="x", colors=extra_axis_color, which="both")
        ax_month.spines["bottom"].set_edgecolor(extra_axis_color)


preset_maps = {}
local_variables = dict(**locals())
for k, v in local_variables.items():
    try:
        assert issubclass(v, Map)
        preset_maps[k] = v
    except (AssertionError, TypeError):
        continue
for k in ["Map", "BubbleMap", "ErrorMap"]:
    preset_maps.pop(k)
