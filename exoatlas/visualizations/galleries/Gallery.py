from ..maps import *


def clean_maps(maps):
    """
    Make sure that a collection of Maps
    is stored as a dictionary.

    Users might sometimes provide a list
    or sometimes a dictionary. This helper
    makes sure it's a dictionary.
    """
    if isinstance(maps, dict):
        # if it's a dictionary, make su
        maps_as_dictionary = maps
    else:
        # otherwise, assume it's a list
        maps_as_dictionary = {p.label: p for p in maps}

    # make deepy copies of all maps to avoid connections between galleries
    for k, v in maps_as_dictionary.items():
        assert isinstance(v, Map)
        maps_as_dictionary[k] = v.copy()
    return maps_as_dictionary


class Gallery:
    def __init__(
        self,
        maps=[],
        label=None,
        **kw,
    ):
        """
        Initialize a Gallery of Maps.

        The default will display a set of Maps in one
        row or column. For building more complicated
        galleries of maps you might want to start
        from here and overwrite two methods:
        - `setup_maps()` to connect maps to axes
        - `refine_maps()` to fuss with details
        The docstrings for these two methods briefly
        sketch what each is for.

        Parameters
        ----------
        maps : dict, list
            A dictionary or list of Maps that should be
            included in this Gallery. If Maps are not
            specified here, they must be defined inside of
            `.setup_maps()`.
        label : str, None
            A string to label this Gallery, which will
            mostly appear in filenames for saved plots.
        **kw : dict
            All additional keywords will be passed to
            `.setup_maps()` to decide how the subplots
            will be arranged.
        """
        self.maps = clean_maps(maps)

        # create the figure + axes (but leave them empty)
        self.setup_maps(**kw)

        # make sure this gallery has a label
        self.label = label or "+".join([p.label for p in self.maps.values()])

    def create_subplots(
        self,
        nrows=1,
        ncols=1,
        sharex="col",
        sharey="row",
        figsize=(9, 6),
        width_ratios=None,
        height_ratios=None,
        wspace=None,
        hspace=None,
    ):
        """
        This helper just sets some reasonable defaults
        for constructing the figure and axes into
        which Maps will be organized.

        Parameters
        ----------
        nrows : int
            How many rows?
        ncols : int
            How many columns?
        sharex :
            How should the x-axis be shared across maps?
            (default 'col') shares along columns
        sharey : bool, str
            How should the y-axis be shared across maps?
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
        If any axes aren't connected to maps, remove them.
        """
        for a in np.atleast_1d(self.ax.flatten()):
            is_used = np.any([a == p.ax for p in self.maps.values()])
            if is_used == False:
                a.axis("off")

    def setup_maps(self, horizontal=True, **kw):
        """
        Create the figure + axes, and
        link the axes to Map objects.

        You might want to write over this for a
        more complicated Gallery of Maps!

        This should create two internal variables:
            `self.fi` = the figure
            `self.ax` = a 1D or 2D grid of axes
        It should also make sure that every map
        to plot in `self.maps` is connected to an
        ax in `self.ax`, either by creating the map
        with a `ax=...` keyword, or by setting the
        `.ax` attribute of an already existing Map.

        Parameters
        ----------
        horizontal : bool
            If True, make a row of Maps.
            If False, make a column of Maps.
        **kw : dict
            All other keywords will be passed to
            `.create_subplots()` for setting up the
            figure and axes. See its docstring.
        """

        # decide rows and columns
        N = len(self.maps)
        self.horizontal = horizontal
        if horizontal:
            rows, cols = 1, N
        else:
            rows, cols = N, 1

        # create the subplots
        self.fi, self.ax = self.create_subplots(nrows=rows, ncols=cols, **kw)

        # attach maps to each ax
        for i, k in enumerate(self.maps):
            self.maps[k].ax = np.atleast_1d(self.ax.flatten())[i]

    def refine_maps(self):
        """
        Make small changes to Maps, after data are plotted.

        You might want to write over this for a
        more complicated Gallery of Maps!

        We might want to plot some extra curves on
        a map, or fuss with its axis labels, or
        nudge its ticks, or add a legend. This
        function gets called after the data have
        been built up into the plot, so those
        modifications can be made.

        By default, this simply removes axis labels
        from
        """
        if self.horizontal:
            for a in self.ax[1:]:
                # plt.setp(a.get_yticklabels(), visible=False)
                a.set_ylabel("")
        else:
            for a in self.ax[:-1]:
                # plt.setp(a.get_xticklabels(), visible=False)
                a.set_xlabel("")

    def build(self, pops, save=False, steps=True, format="png", savefig_kw={}):
        """
        Populate all Maps in a Gallery,
        using a set of populations.

        Parameters
        ----------
        pops : dict
            The populations to plot.
        save : bool
            Should we save the figure(s) to files?
        format : str
            What kind of graphics file should be saved?
            Any format allowed by `plt.savefig` is fine;
            we recommend `png` or `pdf`.
        savefig_kw : dict
            Keywords to pass to `plt.savefig` for saving figures.

        Returns
        -------
        self : Gallery
            In case you want to modify the plot after it's
            generated, this returns the entire Gallery.
        """

        # put data into the axes
        for k, p in self.maps.items():
            p.build(pops=pops, legend=False)

        # remove any empty axes
        self.remove_unused_axes()

        # make small tweaks to the maps
        self.refine_maps()

        # save figure
        if save:
            filename = f"{self.label}"
            plt.savefig(f"{filename}.{format}", **savefig_kw)

        return


class TransitGallery(Gallery):
    def __init__(
        self,
        maps=[
            Mass_x_Radius(),
            Flux_x_Radius(),
            StellarRadius_x_PlanetRadius(),
            Distance_x_Radius(),
        ],
        **kw,
    ):
        Gallery.__init__(self, maps=maps, **kw)
