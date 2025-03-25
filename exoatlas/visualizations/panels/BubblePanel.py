from .Panel import *

__all__ = ["BubblePanel"]

default_size = plt.matplotlib.rcParams["lines.markersize"] ** 2


class BubblePanel(Panel):
    """
    BubblePanel is a general wrapper for making scatter plots
    where planets are represented as bubbles that can have
    informative sizes and/or colors. As bubbles can have
    x, y, size, color; a BubblePanel can be a decent way
    to represent up to four different dimensions.
    """

    color = None
    size = None

    def __init__(
        self,
        xaxis=None,
        yaxis=None,
        size=None,
        color=None,
        size_normalization=100,
        size_vmin=None,
        size_vmax=None,
        size_scale=None,
        cmap="plasma",
        color_vmin=None,
        color_vmax=None,
        color_scale=None,
        **kw,
    ):
        """
        Initialize a plotting panel for drawing bubbles.

        Parameters
        ----------
        xaxis : Plottable, str, None
            What should the x-positions encode?
        yaxis : Plottable, str, None
            What should the y-positions encode?
        size : Plottable, str, float, None
            What should the sizes of points be or encode?
        color : Plottable, str, float, None
            What should the colors of points be or encode?
        size_normalization : float
            If sizes depend on quantities, how big should
            the symbols be for a `normalized_value` of 1?
            Sizes are scaled relative to the default
            symbol area in `matplotlib`, so a value of
            `size_normalization=1` will give default areas.
        size_vmin : float, None
            If the sizes depend on quantities, what
            unnormalized data value corresponds to the
            a symbol area of 0? (This will overwrite
            `lim` set for the size plottable.)
        size_vmax : float, None
            If the sizes depend on quantities, what
            unnormalized data value corresponds to the
            default symbol area? (This will overwrite
            `lim` set for the size plottable.)
        size_scale: float, None
            If the sizes depend on quantities, what
            scaling should be applied when normalizing
            values? Options are `linear` and `log`.
            (This will overwrite `scale` set for the
            size plottable.)
        cmap : str, cmap from plt.matplotlib.cm
            If the colors depend on quantities, what cmap
            should be used for them? Colors will be based
            on `normalized_value` from the Plottable
        color_vmin : float, None
            If the colors depend on quantities, what
            unnormalized data value corresponds to the
            bottom of the color map? (This will overwrite
            `lim` set for the color plottable.)
        color_vmax : float, None
            If the colors depend on quantities, what
            unnormalized data value corresponds to the
            top of the color map? (This will overwrite
            `lim` set for the color plottable.)
        color_scale: float, None
            If the colors depend on quantities, what
            scaling should be applied when normalizing
            values? Options are `linear` and `log`.
            (This will overwrite `scale` set for the
            color plottable.)
        **kw : dict
            Other keywords will be passed on to *all*
            Panel/Plottable initializations (which may
            include x, y, size, and color). If you need
            more fine-grained control over which axis
            gets which keyword, consider initializing
            those panels one-by-one.
        """

        # initialize the basics of the panel with the plottable axes
        Panel.__init__(self, xaxis=xaxis, yaxis=yaxis, **kw)

        # set up how we should scale the sizes of points
        size_to_use = clean_plottable(size or self.size, **kw)
        if isinstance(size_to_use, Plottable):
            self.plottable["size"] = size_to_use
            # overwrite plottable defaults, if requested
            if size_vmin is not None:
                self.plottable["size"].lim[0] = size_vmin
            if size_vmax is not None:
                self.plottable["size"].lim[1] = size_vmax
            if size_scale is not None:
                self.plottable["size"].scale = size_scale
        else:
            self.static_size = size_to_use

        # make sure a size normalization has been defined (relative to default)
        self.size_normalization = size_normalization

        # set up how we should set the colors of points
        color_to_use = clean_plottable(color or self.color, **kw)
        if isinstance(color_to_use, Plottable):
            self.plottable["color"] = color_to_use
            # overwrite plottable defaults, if requested
            if color_vmin is not None:
                self.plottable["color"].lim[0] = color_vmin
            if color_vmax is not None:
                self.plottable["color"].lim[1] = color_vmax
            if color_scale is not None:
                self.plottable["color"].scale = color_scale
        else:
            self.static_color = color_to_use

        # if an actual cmap was provided, use it
        if isinstance(cmap, plt.matplotlib.colors.Colormap):
            self.cmap = cmap
        # otherwise, treat the cmap as a string key
        else:
            self.cmap = plt.matplotlib.colormaps[cmap]

        # apply (x,y) axis labels, scales, limits appropriately
        for axis in "xy":
            for attribute in ["label", "scale", "lim"]:
                setattr(
                    self, f"{axis}{attribute}", getattr(self.plottable[axis], attribute)
                )

    def get_sizes(self):
        """
        The sizes of the bubbles.

        Returns
        -------
        s : an input for plt.scatter
            Either a single scalar, or an array with variable
            sizes for each bubble according to some quantity.
        """

        # should we ignore any variable size instructions?
        if self.pop.respond_to_size == False:
            size = self.pop._plotkw.get("s", None)
        # if desired, set variable sizes
        elif "size" in self.plottable:
            # get the raw values for the sizes
            x = self.plottable["size"].normalized_value(self.pop)
            # calculate the normalized size
            size = default_size * self.size_normalization * x
        # otherwise, set a single size
        else:
            # get default, first from pop and then from panel
            size = getattr(self.pop, "s", None) or self.static_size

        # return a valid input to plt.scatter(s=...)
        return size

    def get_colors(self):
        """
        The colors of the bubbles.

        Returns
        -------
        c : an input for plt.scatter
            Either a single color, or an array with variable
            colors for each bubble according to some quantity.
        """

        # should we ignore any variable color instructions?
        if self.pop.respond_to_color == False:
            color = self.pop.color
        # should we use a variable color?
        elif "color" in self.plottable:
            normalized = self.plottable["color"].normalized_value(self.pop)
            color = self.cmap(normalized)
        # finally, should we just use a default color?
        else:
            # get default, first from pop and then from panel
            color = getattr(self.pop, "color", None) or self.static_color

        # return a valid input to any one of the following:
        #   plt.scatter(c=...)
        #   plt.scatter(edgecolors=...)
        #   plt.scatter(facecolors=...)
        return color

    def kw(self, pop=None, **kw):
        """
        Do a little decision-making about the plotting keyword
        arguments, pulling defaults from each population where
        needed.

        Parameter
        ---------
        pop : Population, str, int
            The population to plot. This could be either an
            actual Population object, or a key referring to
            an element of the `self.populations` dictionary.
        **kw : dict
            All other keywords will be directed toward
            overwriting individual population defaults.
        """

        # make sure we're pointing at the right population
        if pop != None:
            self.point_at(pop)

        # define some default keywords, which can be over-written
        default = dict(
            s=self.get_sizes(),
            marker=self.pop.marker,
            linewidth=self.pop.linewidth,
            alpha=self.pop.alpha,
            zorder=self.pop.zorder,
            label=self.pop.label,
        )

        # sort out whether faces and/or edges should get color
        c = self.get_colors()
        if self.pop.filled:
            default["facecolors"] = c
        else:
            default["facecolors"] = "none"
        if self.pop.outlined:
            default["edgecolors"] = c
        else:
            default["edgecolors"] = "none"

        # if any other keywords are provided, overwrite these defaults
        default.update(**kw)

        return default

    '''def plot(self, key, ax=None, annotate_kw={}, **kw):
        """
        Add the points for a particular population to this panel.

        Parameters
        ----------
        key : str
            The population (as an item in the self.populations dictionary) to add.
        ax :
            Into what ax should we place this plot?
            If None, use default.
        annotate_kw : dict
            Keywords for labeling the planet names.
        **kw : dict
            Any extra keywords will be passed on to `scatter`
        """

        # focus attention on that population
        self.point_at(key)

        # make sure we're plotting into the appropriate axes
        try:
            plt.sca(self.ax)
        except AttributeError:
            self.setup_axes(ax=ax)

        # add the scattered points
        self.scattered[key] = self.ax.scatter(self.x, self.y, **self.kw(key, **kw))

        # set the scales, limits, labels
        self.finish_plot(annotate_kw=annotate_kw)
'''
