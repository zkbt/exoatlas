from .shoreline_maps import *
from ..animation import get_animation_writer
from ..plottables import *
from ...models import *


class Shoreline_SemimajorAxis_x_StellarLuminosity_x_Radius(ShorelineErrorMap):
    """ """

    xaxis = SemimajorAxis(lim=[1e-3, 1e1])
    yaxis = StellarLuminosity(lim=[1e-4, 1e1])
    sliceaxis = Radius(lim=[0.3, 2.0] * u.Rearth)

    def _convert_radius_to_relative_escape_velocity(self, radius):
        """
        Use a mass-radius relation to convert from
        log(radius) to log(relative escape velocity).
        """
        lnR = np.log(radius.to_value(u.Rearth))
        mass_radius_slope = 3.353
        mass_radius_intercept = 0.011
        lnM = mass_radius_slope * lnR + mass_radius_intercept
        lnV = 0.5 * (lnM - lnR)
        return np.exp(lnV)

    def _get_slice_relative_escape_velocities(self, N=100):
        """
        Generate log10(relative escape velocity) for this radius slice,
        including both one value for the "center" of the slice

        Parameters
        ----------
        N : int
            The number of samples to generate.

        Returns
        -------
        log10v : float
            A grid of log10(relative escape velocity) values
            that span the slice from limit to limit.
        """

        # uniform grid in R
        # R = np.linspace(*self.plottable['slice'].lim, num=N)

        # uniform grid in log R
        lnR = np.linspace(
            *np.log(self.plottable["slice"].lim.to_value(u.R_earth)), num=N
        )
        R = np.exp(lnR) * u.R_earth

        # use mass-radius relation to convert to escape velocity
        relative_escape_velocity = self._convert_radius_to_relative_escape_velocity(R)
        return relative_escape_velocity

    def plot_hz_earth(
        self, color="seagreen", alpha=0.75, ax=None, annotate=False, **kw
    ):
        """
        Plot the Earth equivalent habitable zone.

        Parameters
        ----------
        **kw : dict
            All keywords will be passed to `plt.plot`
        """

        plt.sca(ax or self.ax)

        # create grid of luminosities
        L = np.logspace(-4, 2) * u.Lsun

        # plot Earth-equivalent instellation distance
        a_ee = 1 * u.AU * (L / u.Lsun) ** 0.5
        plt.plot(a_ee, L, color=color, alpha=alpha, **kw)

        if annotate:
            plt.text(
                1.8,
                7,
                "habitable\nzone   ",
                va="top",
                ha="right",
                fontsize=7,
                color="seagreen",
            )

    def plot_hz_kopparapu(
        self,
        color="seagreen",
        alpha=0.25,
        linewidth=0,
        zorder=-100,
        ax=None,
        annotate=False,
        **kw,
    ):
        """
        Plot the Kopparapu et al. (2013) habitable zone.

        This includes the Kopparapu et al. (2013)
        HZ as a function of stellar Teff (translated
        from luminosity via the Mamajek table).

        It makes no use of the planet's mass,
        just the stellar effective temperature.

        Parameters
        ----------
        **kw : dict
            All keywords will be passed to `plt.fill_betweenx`
        """

        plt.sca(ax or self.ax)

        # create grid of luminosities
        L = np.logspace(-4, 2) * u.Lsun

        # define the HZ functions from Kopparapu
        f_inner = make_hz("moist-greenhouse")
        f_outer = make_hz("maximum-greenhouse")

        # translate from luminosity to Teff
        m = Mamajek()
        logT = m.tofrom("logT")("logL")(np.log10(L / u.Lsun))
        Teff = 10**logT

        # calculate habitable zone distance (in AU)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            a_inner = (L.to_value("Lsun") / f_inner(Teff)) ** 0.5
            a_outer = (L.to_value("Lsun") / f_outer(Teff)) ** 0.5

        # plot the swath for the HZ
        plt.fill_betweenx(
            L,
            a_inner,
            a_outer,
            color=color,
            alpha=alpha,
            linewidth=linewidth,
            zorder=zorder,
            **kw,
        )

        if annotate:
            plt.text(
                1.8,
                7,
                "habitable\nzone   ",
                va="top",
                ha="right",
                fontsize=7,
                color="seagreen",
            )

    def plot_stellar_radius(self, color="gray", alpha=0.5, **kw):
        """
        Plot the radius of the star, for reference.
        """

        # create grid of luminosities
        L = np.logspace(-4, 2) * u.Lsun

        # translate luminosity to radius
        m = Mamajek()
        logR = m.tofrom("logR")("logL")(np.log10(L / u.Lsun))
        radius = ((10**logR) * u.Rsun).to(u.AU)
        plt.plot(radius, L, color=color, alpha=alpha, **kw)

    def plot_shoreline_lines(self, N_samples=100, ax=None):
        """
        Plot shoreline locations.
        """

        plt.sca(ax or self.ax)

        # create grid of luminosities
        L = np.logspace(-4, 2) * u.Lsun

        best_parameters = self.best_parameters()
        sampled_parameters = self.sampled_parameters(N_samples)
        log_f_0 = best_parameters["log_f_0"]
        q = best_parameters["q"]
        p = best_parameters["p"]
        log_v = np.log10(self._get_slice_relative_escape_velocities(N_samples))

        log_v_center = np.median(log_v)
        log_a = 0.5 * ((1 - q) * np.log10(L / u.Lsun) - p * log_v_center - log_f_0)
        a_shoreline = 10**log_a * u.AU
        plt.plot(a_shoreline, L, color="gray", alpha=1, linestyle="--")

        for s in range(N_samples):
            log_f_0 = sampled_parameters["log_f_0"].values[s]
            q = sampled_parameters["q"].values[s]
            p = sampled_parameters["p"].values[s]

            log_a = 0.5 * ((1 - q) * np.log10(L / u.Lsun) - p * log_v[s] - log_f_0)
            a_shoreline = 10**log_a * u.AU
            plt.plot(a_shoreline, L, color="gray", alpha=0.1)

    def plot_shoreline_probability(
        self, N_samples=100, colormesh=True, contour=True, ax=None, annotate=True
    ):

        plt.sca(ax or self.ax)

        best_parameters = self.best_parameters()
        sampled_parameters = self.sampled_parameters(N_samples)
        log_v = np.log10(self._get_slice_relative_escape_velocities(N_samples))

        # set up 2D grid for background colormap
        log_a_1d = np.linspace(-3, 1, 1000)
        log_L_1d = np.linspace(-4, 1, 1000)
        log_a_2d, log_L_2d = np.meshgrid(log_a_1d, log_L_1d)

        # generate lots of 2D probability maps
        samples_of_log_P_2d = []
        for s in range(N_samples):

            # set up parameters for this sample iteration
            input_parameters = dict(
                p=sampled_parameters["p"].values[s],
                q=sampled_parameters["q"].values[s],
                log_f_0=sampled_parameters["log_f_0"].values[s],
                ln_w=sampled_parameters["ln_w"].values[s],
            )

            log_f_2d = log_L_2d - 2 * log_a_2d
            input_data = dict(log_v=log_v[s], log_f=log_f_2d, log_L=log_L_2d)

            # calculate the 2D probability and store it as part of
            this_P_2d = self.probability_of_atmosphere(**input_data, **input_parameters)
            samples_of_log_P_2d.append(this_P_2d)

        log_P_2d = np.mean(samples_of_log_P_2d, axis=0)

        if colormesh:
            self._shoreline_probability_imshow = plt.pcolormesh(
                10**log_a_2d,
                10**log_L_2d,
                log_P_2d,
                cmap=one2another("burlywood", "lightskyblue"),
                alpha=1,
                zorder=-1e9,
                rasterized=True,
            )

        if contour:
            self._shoreline_probability_contours = plt.contour(
                10**log_a_2d,
                10**log_L_2d,
                log_P_2d,
                levels=[0.05, 0.5, 0.95],
                linestyles="--",  # (0, (5,5)),
                alpha=0.5,
                colors=["gray", "black", "gray"],
            )

        if annotate:
            a_text = 10 ** np.interp(0.5, log_P_2d[-1, :], log_a_1d)
            plt.text(
                a_text,
                7,
                "   cosmic\n    shoreline",
                va="top",
                ha="left",
                fontsize=7,
                color="black",
                alpha=0.5,
            )

    def refine(self, **kw):
        self.plot_hz_earth(annotate=True)
        self.plot_hz_kopparapu()
        self.plot_shoreline_probability(annotate=True)


class ShorelineStandardMap(ShorelineErrorMap):

    linear_plottables = dict(
        v=RelativeEscapeVelocity, f=RelativeInstellation, L=RelativeStellarLuminosity
    )

    log_plottables = dict(
        v=LogRelativeEscapeVelocity,
        f=LogRelativeInstellation,
        L=LogRelativeStellarLuminosity,
    )

    def __init__(
        self, order="vfL", posterior=None, vlim=None, flim=None, Llim=None, **kw
    ):
        # decide which quantity is plotted on which axis
        x_var, y_var, slice_var = order

        # make sure they're all different
        assert sorted(order) == ["L", "f", "v"]

        # set default limits
        self.lims = dict(
            v=np.array([-1.5, 1.0]), f=np.array([-4.0, 4.0]), L=np.array([0.5, -3.5])
        )
        # adjust limits if necessary
        if vlim is not None:
            self.lims["v"] = np.array(vlim)
        if flim is not None:
            self.lims["f"] = np.array(flim)
        if Llim is not None:
            self.lims["L"] = np.array(Llim)

        # initialize the map
        super().__init__(
            xaxis=self.linear_plottables[x_var](lim=10 ** self.lims[x_var], **kw),
            yaxis=self.linear_plottables[y_var](lim=10 ** self.lims[y_var], **kw),
            sliceaxis=self.log_plottables[slice_var](lim=self.lims[slice_var], **kw),
            posterior=posterior,
            **kw,
        )

        self.x_var = x_var
        self.y_var = y_var
        self.slice_var = slice_var

        # set a buffer to extend the plots a little beyond the ranges
        self._axes_buffer = 1.5

    def plot_shoreline_probability(
        self,
        N_samples=100,
        colormesh=True,
        contour=True,
        ax=None,
    ):

        plt.sca(ax or self.ax)

        best_parameters = self.best_parameters()
        sampled_parameters = self.sampled_parameters(N_samples)

        # set up 1D grid for background colormap
        smooth_log_grids = {}
        for k in self.lims:
            smooth_log_grids[k] = np.linspace(
                min(self.lims[k]) - np.log10(self._axes_buffer) * 3,
                max(self.lims[k]) + np.log10(self._axes_buffer) * 3,
                2000,
            )

        # set up 3D grid for background colormap
        log_x_2d, log_y_2d = np.meshgrid(
            smooth_log_grids[self.x_var], smooth_log_grids[self.y_var]
        )

        # calculate 2D image of the probability of having an atmosphere
        # (we could probably do this with vmap faster...)
        samples_of_log_P_2d = []
        log_slice_samples = np.linspace(*self.plottable["slice"].lim, N_samples)

        # loop over samples
        for s in range(N_samples):

            # for this map, define the smooth 2D grids for the background probability
            inputs_with_exoatlas_names = {
                self.x_var: log_x_2d,
                self.y_var: log_y_2d,
                self.slice_var: log_slice_samples[s],
            }
            inputs_with_model_names = {
                f"log_{k}": inputs_with_exoatlas_names[k] for k in "vfL"
            }

            # set up parameters for this sample iteration
            input_parameters = dict(
                p=sampled_parameters["p"].values[s],
                q=sampled_parameters["q"].values[s],
                log_f_0=sampled_parameters["log_f_0"].values[s],
                ln_w=sampled_parameters["ln_w"].values[s],
            )

            # calculate the 2D probability and store it as part of
            this_P_2d = self.probability_of_atmosphere(
                **inputs_with_model_names, **input_parameters
            )
            samples_of_log_P_2d.append(this_P_2d)

        log_P_2d = np.mean(samples_of_log_P_2d, axis=0)

        self._shoreline_probability_imshow = plt.pcolormesh(
            10**log_x_2d,
            10**log_y_2d,
            log_P_2d,
            cmap=one2another("burlywood", "lightskyblue"),
            alpha=1,
            zorder=-1e9,
            rasterized=True,
        )
        self._shoreline_probability_contours = plt.contour(
            10**log_x_2d,
            10**log_y_2d,
            log_P_2d,
            levels=[0.05, 0.5, 0.95],
            linestyles="--",  # (0, (5,5)),
            alpha=0.25,
            colors=["gray", "black", "gray"],
        )

    def label_flux_limits(
        self, limits={"magma ocean": 1700 * u.K, "$\sf CO_2$ freezes": 194 * u.K}
    ):

        # figure out which axis is instellation

        for k in self.plottable:
            if isinstance(self.plottable[k], RelativeInstellation):
                instellation_axis = k

        try:
            instellation_axis
        except NameError:
            return
            # print(
            #    "This map doesn't seem to have any instellation axis,"
            #    "so we can't plot flux limits on it."
            # )

        # calculate flux limits
        def calculate_S(T, f=1 / 4, A=0):
            # translate temperature into instellation
            unnormalized_flux = con.sigma_sb * T**4 / f / (1 - A)
            earth_insolation = (1 * u.Lsun / 4 / np.pi / u.AU**2).to(u.W / u.m**2)
            return (unnormalized_flux / earth_insolation).decompose()

        for k in limits:
            S = calculate_S(limits[k])
            line_kw = dict(color="dimgray", alpha=0.25, linewidth=1)
            text_kw = dict(color="dimgray", va="bottom", alpha=1, fontsize=6)
            if instellation_axis == "x":
                plt.axvline(S, **line_kw)
                plt.text(
                    S,
                    self.plottable["y"].lim[1],
                    f"  {k}",
                    ha="right",
                    rotation=90,
                    **text_kw,
                )
            elif instellation_axis == "y":
                plt.axhline(S, **line_kw)
                plt.text(
                    self.plottable["x"].lim[0],
                    S,
                    f"  {k}",
                    ha="left",
                    **text_kw,
                )
            else:
                return

    def refine(self, probability=True, limits=False, **kw):
        """
        Refine a standard shoreline Map.

        Parameters
        ----------
        probability : bool
            Should we imshow the probability of an atmosphere?
        label_flux_limits : bool
            Should we draw flux limits on the plot?
        **kw : dict
            Additional keywords will be ignored.

        """
        if probability:
            self.plot_shoreline_probability()

        if limits:
            self.label_flux_limits()

        # enforce axis limits
        plt.xscale("log")
        plt.yscale("log")
        xlim = self.plottable["x"].lim
        ylim = self.plottable["y"].lim
        plt.xlim(
            min(xlim) / self._axes_buffer,
            max(xlim) * self._axes_buffer,
        )
        plt.ylim(
            min(ylim) / self._axes_buffer,
            max(ylim) * self._axes_buffer,
        )


_teq_limits = [194, 1700] * u.K
_flux_limits = 4 * con.sigma_sb * _teq_limits**4
_earth_insolation = (1 * u.Lsun / 4 / np.pi / u.AU**2).to(u.W / u.m**2)
_f_lim = list((_flux_limits / _earth_insolation).value)
