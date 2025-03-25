from ..panels import *


def clean_panels(panels):
    """
    Make sure that a collection of Panels
    is stored as a dictionary.

    Users might sometimes provide a list
    or sometimes a dictionary. This helper
    makes sure it's a dictionary.
    """
    if isinstance(panels, dict):
        # if it's a dictionary, make su
        panels_as_dictionary = panels
    else:
        # otherwise, assume it's a list
        panels_as_dictionary = {p.label: p for p in panels}
    for p in panels_as_dictionary.values():
        assert isinstance(p, Panel)


class Gallery:

    def __init__(self, panels=[]):
        """
        Initialize a Gallery of Maps.
        """
        self.panels = clean_panels(panels)

        # create the figure + axes (but leave them empty)
        self.setup_panels()

    def create_subplots(
        self,
        nrows=1,
        ncols=1,
        sharex="col",
        sharey="row",
        figsize=(12, 8),
        width_ratios=None,
        height_ratios=None,
        wspace=None,
        hspace=None,
    ):
        """
        This helper just sets some reasonable defaults
        for constructing the figure and axes into
        which Panels will be organized.

        Parameters
        ----------
        nrows : int
            How many rows?
        ncols : int
            How many columns?
        sharex :
            How should the x-axis be shared across panels?
            (default 'col') shares along columns
        sharey : bool, str
            How should the y-axis be shared across panels?
                (default 'row') shares along rows
        figsize : tuple
            How big should the figure be?
        width_ratios : tuple
            How wide should the columns be relative to each other?
        height_ratios : tuple
            How tall should the rows be relative to each other?
        wspace : float
            How much space between columns?
        hspace : float
            How much space between rows?
        """
        return plt.subplots(
            nrows=nrows,
            ncols=ncols,
            sharex=sharex,
            sharey=sharey,
            figsize=figsize,
            constrained_layout=True,
            gridspec_kw=dict(
                width_ratios=width_ratios,
                height_ratios=height_ratios,
                wspace=wspace,
                hspace=hspace,
            ),
        )

    def remove_unused_axes(self):
        """
        If any axes aren't connected to panels, remove them.
        """
        for a in self.ax.flatten():
            is_used = np.any([a == p.ax for p in self.panels.values()])
            if is_used == False:
                a.axis("off")

    def setup_panels(self):
        """
        Create the figure + axes, and
        link the axes to Panel objects.

        This should create two internal variables:
            `self.fi` = the figure
            `self.ax` = a 1D or 2D grid of axes
        It should also make sure that every panel
        to plot in `self.panels` is connected to an
        ax in `self.ax`, either by creating the panel
        with a `ax=...` keyword, or by setting the
        `.ax` attribute of an already existing Panel.
        """
        pass

    def refine_panels(self):
        """
        Make small changes to Panels, after data are plotted.

        We might want to plot some extra curves on
        a panel, or fuss with its axis labels, or
        nudge its ticks, or add a legend. This
        function gets called after the data have
        been built up into the plot, so those
        modifications can be made.
        """
        pass

    def build_panels(self, pops):
        """
        Populate all Panels in a Gallery,
        using a set of populations.

        Parameters
        ----------
        pops : dict
            The populations to plot.
        """

        # put data into the axes
        for k, p in self.panels.items():
            p.build(pops=pops, legend=False)

        # remove any empty axes
        self.remove_unused_axes()

        # make small tweaks to the panels
        self.refine_panels()
