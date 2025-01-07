from ..imports import *
from .panels import *


class BuildablePlot(Talker):
    """docstring for BuildablePlot."""

    def __init__(self, pops, subplot=None, stepbystep=True, **kw):
        """

        Parameters
        ----------
        pops : dict
            The populations to plot.
        subplot : matplotlib.gridspec.SubplotSpec
            A gridspec specifier to plot into.
        stepbystep : bool
            Should this plot be saved out in stages
            of each panel being added one by one?
        **kw : dictionary
            All other keywords will be passed to the buildable plot function.
        """

        self.subplot = subplot
        self.stepbystep = stepbystep
        self.name = "exoplanets"

        pops = clean_pops(pops)

        self.plot(pops, **kw)

    def create_gridspec(self, *args, figsize=(10, 5), **kwargs):
        """
        Wrapper for gridspec, so we can make this plot either on its own
        or inside a panel of another gridspec location.
        """
        if self.subplot is None:
            self.figure = plt.figure(figsize=figsize)
            f = plt.matplotlib.gridspec.GridSpec
            # kwargs["figure"] = self.figure
            return f(*args, **kwargs)
        else:
            for k in ["bottom", "top", "left", "right"]:
                if k in kwargs:
                    kwargs.pop(k)
            f = plt.matplotlib.gridspec.GridSpecFromSubplotSpec
            return f(*args, subplot_spec=self.subplot, **kwargs)

    # def pause(self):
    #    if self.stepbystep:
    #
    # def write(self):
    #    label = '-'.join(list(pops.keys()))
    #    plt.savefig(f'{label}.pdf')


class physical_summary(BuildablePlot):
    def plot(self, pops):
        gs = self.create_gridspec(
            2,
            3,
            width_ratios=[1, 4, 1],
            height_ratios=[1, 3],
            wspace=0.04,
            hspace=0.04,
            bottom=0.21,
            top=0.98,
            right=0.98,
            left=0.08,
        )

        fr = FluxRadius()
        fr.build(pops=pops, ax=plt.subplot(gs[1, 1]))
        fr.plot_hz()
        fr.ticks_simplify_exponents("y")
        fr.add_teqaxis()
        if len(pops) > 1:
            fr.add_legend(frameon=False, bbox_to_anchor=(1, 1), loc="lower left")
        fr.remove_ylabel()

        mr = MassRadius()
        mr.build(pops=pops, ax=plt.subplot(gs[1, 0], sharey=fr.ax))
        mr.plot_both_seager(zorder=1e9)
        mr.ticks_enforce_multiple_oom("x")
        mr.ticks_simplify_exponents("xy")

        er = FluxEscape()
        er.build(pops=pops, ax=plt.subplot(gs[0, 1], sharex=fr.ax))
        er.remove_xlabel()
        er.plot_constant_lambda()

        # plot the
        sr = StellarRadiusPlanetRadius()
        sr.build(pops=pops, ax=plt.subplot(gs[1, 2], sharey=fr.ax))
        sr.remove_ylabel()
        mr.ticks_simplify_exponents("y")

        self.panels = {p.name: p for p in [fr, mr, er, sr]}


class observable_summary(BuildablePlot):
    def plot(
        self,
        pops,
        observables=["depth", "emission", "reflection", "transmission"],
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
            pr = DistanceRadius(**kw, **wavelengthkw)
            pr.build(pops=pops, ax=plt.subplot(gs[1, column]))

            db = DistanceBrightness(telescope_name=telescope_name, **kw, **wavelengthkw)
            db.build(pops=pops, ax=plt.subplot(gs[0, column]))
            db.remove_xlabel()
            column += 1

        if "depth" in observables:
            dr = DepthRadius(
                size=DepthSNR(telescope_name=telescope_name, per_transit=per_transit),
                **kw,
                **wavelengthkw
            )
            dr.build(pops=pops, ax=plt.subplot(gs[1, column]))
            dr.remove_ylabel()

            x = DepthBrightness(
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
            er = EmissionRadius(
                size=EmissionSNR(
                    telescope_name=telescope_name, per_transit=per_transit
                ),
                **kw,
                **wavelengthkw
            )
            er.build(pops=pops, ax=plt.subplot(gs[1, column]))
            er.remove_ylabel()

            x = EmissionBrightness(
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
            rr = ReflectionRadius(
                size=ReflectionSNR(
                    telescope_name=telescope_name, per_transit=per_transit
                ),
                **kw,
                **wavelengthkw
            )
            rr.build(pops=pops, ax=plt.subplot(gs[1, column]))
            rr.remove_ylabel()

            x = ReflectionBrightness(
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
            tr = TransmissionRadius(
                size=TransmissionSNR(
                    telescope_name=telescope_name, per_transit=per_transit
                ),
                **kw,
                **wavelengthkw
            )
            tr.build(pops=pops, ax=plt.subplot(gs[1, column]))
            tr.remove_ylabel()

            x = TransmissionBrightness(
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
            fr = FluxRadius(**kw)
            fr.xlabel = "Bolometric Flux Received\n(relative to Earth)"
            fr.build(pops=pops, ax=plt.subplot(gs[1, column]))
            fr.remove_ylabel()
            column += 1
