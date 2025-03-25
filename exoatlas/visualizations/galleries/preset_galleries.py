from .Gallery import *


class physical_summary(Gallery):
    """
    This Gallery is designed to give a quick
    look at populations of transiting exoplanets,
    focusing on some of the physical attributes
    that strongly influence how planets work.
    """

    def setup_panels(self):
        """
        Create the figure + axes, and
        link the axes to Panel objects.
        """
        # create the subplots
        self.fi, self.ax = self.create_subplots(
            nrows=2, ncols=3, width_ratios=[2, 4, 1], height_ratios=[1, 3]
        )

        # attach panels to each ax
        self.panels["flux-radius"] = FluxRadius(ax=self.ax[1, 1])
        self.panels["mass-radius"] = MassRadius(ax=self.ax[1, 0])
        self.panels["flux-escape"] = FluxEscape(ax=self.ax[0, 1])
        self.panels["radius-radius"] = StellarRadiusPlanetRadius(ax=self.ax[1, 2])
        self.panels["mass-escape"] = BubblePanel(
            xaxis=Mass, yaxis=EscapeVelocity, ax=self.ax[0, 0]
        )

    def refine_panels(self):
        """
        Make small changes to Panels, after data are plotted.
        """

        #
        plt.sca(self.panels["flux-radius"].ax)
        self.panels["flux-radius"].plot_hz()
        self.panels["flux-radius"].ticks_simplify_exponents("y")
        self.panels["flux-radius"].add_teqaxis()
        self.panels["flux-radius"].remove_ylabel()
        self.panels["flux-radius"].add_legend(frameon=False)

        plt.sca(self.panels["mass-radius"].ax)
        self.panels["mass-radius"].plot_both_seager(zorder=1e9)
        # self.panels['mass-radius'].ticks_enforce_multiple_oom("x")
        # self.panels['mass-radius'].ticks_simplify_exponents("xy")

        plt.sca(self.panels["mass-escape"].ax)
        self.panels["mass-escape"].remove_xlabel()

        plt.sca(self.panels["flux-escape"].ax)
        self.panels["flux-escape"].remove_xlabel()
        self.panels["flux-escape"].remove_ylabel()
        self.panels["flux-escape"].plot_constant_lambda()

        plt.sca(self.panels["radius-radius"].ax)
        self.panels["radius-radius"].remove_ylabel()
