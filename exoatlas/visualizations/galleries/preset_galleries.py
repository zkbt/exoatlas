from .Gallery import *


class physical_summary(Gallery):
    """
    This Gallery is designed to give a quick
    look at populations of transiting exoplanets,
    focusing on some of the physical attributes
    that strongly influence how planets work.
    """

    def setup_maps(self):
        """
        Create the figure + axes, and
        link the axes to Map objects.
        """
        # create the subplots
        self.fi, self.ax = self.create_subplots(
            nrows=2, ncols=3, width_ratios=[2, 4, 1], height_ratios=[1, 3]
        )

        # attach maps to each ax
        self.maps = {}
        self.maps["flux_x_radius"] = Flux_x_Radius(ax=self.ax[1, 1])
        self.maps["mass_x_radius"] = Mass_x_Radius(ax=self.ax[1, 0])
        self.maps["flux_x_escape"] = Flux_x_EscapeVelocity(ax=self.ax[0, 1])
        self.maps["radius_x_radius"] = StellarRadius_x_PlanetRadius(ax=self.ax[1, 2])
        self.maps["mass_x_escape"] = BubbleMap(
            xaxis=Mass, yaxis=EscapeVelocity, ax=self.ax[0, 0]
        )

    def refine_maps(self):
        """
        Make small changes to Maps, after data are plotted.
        """

        #
        plt.sca(self.maps["flux_x_radius"].ax)
        self.maps["flux_x_radius"].plot_hz()
        self.maps["flux_x_radius"].ticks_simplify_exponents("y")
        self.maps["flux_x_radius"].add_teqaxis()
        self.maps["flux_x_radius"].remove_ylabel()
        self.maps["flux_x_radius"].add_legend(frameon=False, fontsize=8)

        plt.sca(self.maps["mass_x_radius"].ax)
        self.maps["mass_x_radius"].plot_both_seager(zorder=1e9)
        # self.maps['mass-radius'].ticks_enforce_multiple_oom("x")
        # self.maps['mass-radius'].ticks_simplify_exponents("xy")

        plt.sca(self.maps["mass_x_escape"].ax)
        self.maps["mass_x_escape"].remove_xlabel()

        plt.sca(self.maps["flux_x_escape"].ax)
        self.maps["flux_x_escape"].remove_xlabel()
        self.maps["flux_x_escape"].remove_ylabel()
        self.maps["flux_x_escape"].plot_constant_lambda()

        plt.sca(self.maps["radius_x_radius"].ax)
        self.maps["radius_x_radius"].remove_ylabel()


'''class observable_summary(BuildablePlot):
    def plot(
        self,
        pops,
        observables=["depth", "emission", "transmission"],
        note="",
        telescope_name="JWST",
        per_transit=False,
        **wavelengthkw
    ):
        gs = self.create_gridspec(
            2,
            len(observables) + 2,
            wspace=0.05,
            hspace=0.05,
            bottom=0.21,
            top=0.98,
            right=0.98,
            left=0.08,
            figsize=(10, 8),
        )

        if note == "trim":
            for k, p in pops.items():
                pops[k] = p[p.radius < 2 * u.Rearth]

        kw = dict()
        if note != "bw":
            kw["color"] = "log_relative_insolation"
            kw["vmin"] = -1
            kw["vmax"] = 5
            for k, x in pops.items():
                if k == "solarsystem":
                    x.alpha = 1
                else:
                    x.alpha = 0.5

        hidealpha = 0.05
        showalpha = 1.0

        """if 'highlight' in note:
            for k in pops:
                if k in note:
                    pops[k].alpha = showalpha
                else:
                    pops[k].alpha = hidealpha

        if 'toi' in note:
            for k in pops:
                pops[k].alpha = hidealpha
            pops['toi'].alpha = showalpha
        else:
            pops['toi'].alpha = 0.

        if 'all-together' in note:
            pops['toi'].alpha=0.5
            pops['tess'].alpha=0.5
            pops['nontess'].alpha=0.5
        """

        column = 0

        if True:
            pr = Distance_x_Radius(**kw, **wavelengthkw)
            pr.build(pops=pops, ax=plt.subplot(gs[1, column]))

            db = Distance_x_Brightness(
                telescope_name=telescope_name, **kw, **wavelengthkw
            )
            db.build(pops=pops, ax=plt.subplot(gs[0, column]))
            db.remove_xlabel()
            column += 1

        if "depth" in observables:
            dr = Depth_x_Radius(
                size=DepthSNR(telescope_name=telescope_name, per_transit=per_transit),
                **kw,
                **wavelengthkw
            )
            dr.build(pops=pops, ax=plt.subplot(gs[1, column]))
            dr.remove_ylabel()

            x = Depth_x_Brightness(
                telescope_name=telescope_name,
                size=DepthSNR(telescope_name=telescope_name, per_transit=per_transit),
                **kw,
                **wavelengthkw
            )
            x.build(pops=pops, ax=plt.subplot(gs[0, column]))
            x.plot_sigma()
            x.remove_xlabel()
            x.remove_ylabel()
            column += 1

        if "emission" in observables:
            er = Emission_x_Radius(
                size=EmissionSNR(
                    telescope_name=telescope_name, per_transit=per_transit
                ),
                **kw,
                **wavelengthkw
            )
            er.build(pops=pops, ax=plt.subplot(gs[1, column]))
            er.remove_ylabel()

            x = Emission_x_Brightness(
                telescope_name=telescope_name,
                size=EmissionSNR(
                    telescope_name=telescope_name, per_transit=per_transit
                ),
                **kw,
                **wavelengthkw
            )
            x.build(pops=pops, ax=plt.subplot(gs[0, column]))
            x.plot_sigma()
            x.remove_xlabel()
            x.remove_ylabel()
            column += 1

        if "reflection" in observables:
            rr = Reflection_x_Radius(
                size=ReflectionSNR(
                    telescope_name=telescope_name, per_transit=per_transit
                ),
                **kw,
                **wavelengthkw
            )
            rr.build(pops=pops, ax=plt.subplot(gs[1, column]))
            rr.remove_ylabel()

            x = Reflection_x_Brightness(
                telescope_name=telescope_name,
                size=ReflectionSNR(
                    telescope_name=telescope_name, per_transit=per_transit
                ),
                **kw,
                **wavelengthkw
            )
            x.build(pops=pops, ax=plt.subplot(gs[0, column]))
            x.plot_sigma()
            x.remove_xlabel()
            x.remove_ylabel()
            column += 1

        if "transmission" in observables:
            tr = Transmission_x_Radius(
                size=TransmissionSNR(
                    telescope_name=telescope_name, per_transit=per_transit
                ),
                **kw,
                **wavelengthkw
            )
            tr.build(pops=pops, ax=plt.subplot(gs[1, column]))
            tr.remove_ylabel()

            x = Transmission_x_Brightness(
                telescope_name=telescope_name,
                size=TransmissionSNR(
                    telescope_name=telescope_name, per_transit=per_transit
                ),
                **kw,
                **wavelengthkw
            )
            x.build(pops=pops, ax=plt.subplot(gs[0, column]))
            x.plot_sigma()
            x.remove_xlabel()
            x.remove_ylabel()
            column += 1

        if True:
            fr = Flux_x_Radius(**kw)
            fr.xlabel = "Bolometric Flux Received\n(relative to Earth)"
            fr.build(pops=pops, ax=plt.subplot(gs[1, column]))
            fr.remove_ylabel()
            column += 1
'''
