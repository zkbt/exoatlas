from ..imports import *
from .panels import *


def make_sure_panel_is_iniated(x, **kw):
    if isinstance(x, Panel):
        return x
    else:
        return x(**kw)


class MultiPanelPlot(Talker):
    """
    Make and modify row or a column of connected exoatlas plot panels.
    """

    def __init__(
        self,
        panels=[MassRadius, FluxRadius, StellarRadiusPlanetRadius, DistanceRadius],
        horizontal=True,
        figsize=(8, 6),
        **kw
    ):
        """
        Set up the plotting panels.
        Pause before plotting, in case we want to make any modifications.

        Parameters
        ----------

        panels : list
            A list of Panel class definitions to include.
        horizontal : bool
            True = make a row of panels side-by-side
            False = make a column of panels on top of each other
        figsize : tuple
            What's the (width, height) of the figure to make.
        **kw : dict
            Other keywords will be passed to gridspec to set the
            widths, heights, margins, and spacing of the panels.
            These keywords commonly include:

                left, bottom, right, top,
                wspace, hspace,
                width_ratios, height_ratios
        """

        self.panels = {}

        # store list of panel names
        def get_name(x):
            try:
                return x.name
            except:
                try:
                    return x.__name__
                except AttributeError:
                    return x.__class__.__name__

        # create a dictionary of panel objects
        for p in panels:
            this_panel = make_sure_panel_is_iniated(p, **kw)
            key = get_name(this_panel)
            if key in self.panels:
                key += "+"
            self.panels[key] = this_panel
        self.panel_names = list(self.panels.keys())

        # set up the geometry of the figure
        self.horizontal = horizontal
        if horizontal:
            # set up a row of axes
            nrows, ncols = 1, len(self.panels)
            sharex, sharey = False, True
        else:
            # set up a column of axes
            nrows, ncols = len(self.panels), 1
            sharex, sharey = True, False

        # set some default gridspec options, but allow updating
        gridspec_kw = dict(hspace=0.1, left=0.15, right=0.95, bottom=0.15, wspace=0.05)
        gridspec_kw.update(**kw)

        # set up the plotting figure and axes
        self.figure, axgrid = plt.subplots(
            nrows,
            ncols,
            sharex=sharex,
            sharey=sharey,
            figsize=figsize,
            gridspec_kw=gridspec_kw,
        )

        # store the axes in a dictionary, with panel names as keys
        self.ax = {}
        for i, k in enumerate(self.panel_names):
            if nrows * ncols == 1:
                self.ax[k] = axgrid
            else:
                self.ax[k] = axgrid[i]

    def build(self, pops={}, **kw):
        """
        Actually make the plot by building up each panel.

        Parameters
        ----------
        pops : dict
            A dictionary of populations, with keys as labels values as
            initialized populations. For example,
                pops = {'solarsystem':SolarSystem(),
                        'transiting':TransitingExoplanets()}
        **kw : dict
            Any extra keywords will be passed on to all panels' `build`
        """

        self.pops = clean_pops(pops)

        # plot each population in each panel
        for i, k in enumerate(self.panel_names):
            self.panels[k].build(pops=self.pops, ax=self.ax[k], **kw)

        # clean up unnecessary labels
        if self.horizontal:
            for k in self.panel_names[1:]:
                plt.setp(self.ax[k].get_yticklabels(), visible=False)
                self.ax[k].set_ylabel("")
        else:
            for k in self.panel_names[:-1]:
                plt.setp(self.ax[k].get_xticklabels(), visible=False)
                self.ax[k].set_xlabel("")

        return self


FourPanels = MultiPanelPlot
