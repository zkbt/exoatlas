from ..imports import *
from .panels import *


def make_sure_panel_is_iniated(x, **kw):
    if isinstance(x, Panel):
        return x
    else:
        return x(**kw)


class MultiPanelPlot:
    """
    Make and modify row or a column of connected exoatlas plot panels.
    """

    def __init__(
        self,
        panels=[
            MassRadius(),
            FluxRadius(),
            StellarRadiusPlanetRadius(),
            DistanceRadius(),
        ],
        horizontal=True,
        figsize=(12, 8),
        label=None,
        **kw,
    ):
        """
        Set up the plotting panels.

        The initialization creates the subplot axes,
        and then the `.build()` function is used
        to population those with data.

        Parameters
        ----------

        panels : list
            A list of Panel class definitions to include.
        horizontal : bool
            True = make a row of panels side-by-side
            False = make a column of panels on top of each other
        figsize : tuple
            What's the (width, height) of the figure to make.
        label : str
            What is a label to use for this plot for filenames?
        **kw : dict
            Other keywords will be passed to gridspec to set the
            widths, heights, margins, and spacing of the panels.
            These keywords commonly include:

                left, bottom, right, top,
                wspace, hspace,
                width_ratios, height_ratios
        """

        # set up a dictionary for panels and their matplotlib axes
        self.panels = {}
        self.ax = {}

        self.setup_panels()

        # give this MultiPanel plot a label
        self.label = label or "-".join([p.label for p in self.panels.values()])

    def setup_panels(self):
        # populate a dictionary of panel objects
        for p in panels:
            this_panel = make_sure_panel_is_iniated(p, **kw)
            key = this_panel.label
            if key in self.panels:
                key += "+"
            self.panels[key] = this_panel
        self.panel_keys = list(self.panels.keys())

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
        for i, k in enumerate(self.panel_keys):
            if nrows * ncols == 1:
                self.ax[k] = axgrid
            else:
                self.ax[k] = axgrid[i]

            # connect this ax to this panel
            self.panels[k].ax = self.ax[k]

    def plot(self, pop, **kw):
        """
        Add the points for one population to this panel.

        Parameters
        ----------
        pop : Population, str, int
            The population to plot. This could be either an
            actual Population object, or a key referring to
            an element of the `self.populations` dictionary.
        annotate_kw : dict
            Keywords for labeling the planet names.
        **kw : dict
            Any extra keywords will be passed on to
            each panel's `plot()` method
        """
        for k in self.panel_keys:
            self.panels[k].plot(pop, **kw)

        self.refine_panel_labels()

    def refine_panel_labels(self):
        """
        Clean up the look of the panels.

        This mostly removes unnecessary axis labels.
        """

        # clean up unnecessary labels between panels
        if self.horizontal:
            for k in self.panel_keys[1:]:
                plt.setp(self.ax[k].get_yticklabels(), visible=False)
                self.ax[k].set_ylabel("")
        else:
            for k in self.panel_keys[:-1]:
                plt.setp(self.ax[k].get_xticklabels(), visible=False)
                self.ax[k].set_xlabel("")

    def build(
        self,
        pops={},
        save=False,
        format="png",
        savefig_kw={},
        legend=True,
        legend_kw={},
        **kw,
    ):
        """
        Build up this panel by plotting every population into it.

        Parameters
        ----------
        pops : dict
            A dictionary of populations, with keys as labels values as
            initialized populations. For example,

                pops = {'solarsystem':SolarSystem(),
                        'transiting':TransitingExoplanets()}
        save : bool
            Should we save the figure(s) to files?
        format : str
            What kind of graphics file should be saved?
            Any format allowed by `plt.savefig` is fine;
            we recommend `png` or `pdf`.
        savefig_kw : dict
            Keywords to pass to `plt.savefig` for saving figures.
        **kw : dict
            Any extra keywords will be passed on to all panels' `build()`
        """

        # make sure we're working
        self.populations = clean_pops(pops)

        # plot each population in each panel
        for i, k in enumerate(self.panel_keys):
            self.panels[k].build(
                pops=self.populations, legend=False, **kw
            )  # ax=self.ax[k],

        # clean up unnecessary labels
        self.refine_panel_labels()

        # save the final figure
        if save:
            filename = f"{self.label}+" + "+".join(self.populations.keys())
            plt.savefig(f"{filename}.{format}", **savefig_kw)
