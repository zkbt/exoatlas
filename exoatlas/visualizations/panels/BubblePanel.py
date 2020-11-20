from .Panel import *

__all__ = ['BubblePanel']

default_size = plt.matplotlib.rcParams['lines.markersize']**2

class BubblePanel(Panel):
    '''
    BubblePanel is a general wrapper for making scatter plots
    where planets are represented as bubbles that can have
    informative sizes and/or colors.
    '''

    def __init__(self,
                 xaxis=None,
                 yaxis=None,
                 size=None, size_normalization=None,
                 color=None, cmap='plasma', vmin=None, vmax=None, color_normalization=None,
                 **kw):
        '''
        Initialize a plotting panel.

        Parameters
        ----------
        size : PlottableAxis, str, float, None
            What should the sizes of points be or encode?
        size_normalization : float
            If sizes depend on quantities,
            how should they be normalized?
        color : PlottableAxis, str, float, None
            What should the colors of points be or encode?
        cmap : str, cmap from plt.matplotlib.cm
            If the colors depend on quantities,
            what cmap should be used for them?
        vmin : float, astropy.units.quantity.Quantity
            If the colors depend on quantities,
            what should the bottom of the cmap be?
        vmax : float, astropy.units.quantity.Quantity
            If the colors depend on quantities,
            what should the top of the cmap be?
        color_normalization : matplotlib.colors.Normalize
            If color depend on quantities, how should
            the values be normalized. If color_normalization
            is defined, any values provided here for
            vmin and vmax will be ignored.
        **kw : dict
            Other keywords will be passed on to *all*
            Panel/Plottable initializations (which may
            include x, y, size, and color). If you need
            more fine-grained control over which axis
            gets which keyword, consider initializing
            those panels one-by-one.
        '''

        # initialize the basics of the panel with the plottable axes
        Panel.__init__(self, xaxis=xaxis, yaxis=yaxis, **kw)

        # set up how we should scale the sizes of points
        size = clean_axis(size)
        try:
            # try to make a variable size axis
            self.plottable['size'] = size(panel=self, **kw)
            default_size_normalization = self.plottable['size'].size_normalization
        except TypeError:
            # otherwise, use a single size for all points
            self.plottable['size'] = size
            default_size_normalization = 1

        # make sure a size normalization has been defined
        self.size_normalization = size_normalization or default_size_normalization

        # set up how we should set the colors of points
        color = clean_axis(color)
        try:
            # try to make a variable color axis
            self.plottable['color'] = color(panel=self, **kw)
            default_lim = self.plottable['color'].lim
        except TypeError:
            # otherwise, use a single color for all points
            self.plottable['color'] = color
            default_lim = [None, None]

        # if an actual cmap was provided, use it
        if isinstance(cmap, plt.matplotlib.colors.Colormap):
            self.cmap = cmap
        # otherwise, treat the cmap as a string key
        else:
            self.cmap = plt.matplotlib.cm.cmap_d[cmap]

        # make sure the color map limits are set
        self.vmin = vmin or default_lim[0]
        self.vmax = vmax or default_lim[1]

        # if a custom normalization is used, reset vmin + vmax
        self.color_normalization = color_normalization
        if isinstance(self.color_normalization,
                      plt.matplotlib.colors.Normalize):
            # pull the normalization's min/max for information
            self.vmin = color_normalization.vmin
            self.vmax = color_normalization.vmax

        # apply (x,y) axis labels, scales, limits appropriately
        for axis in 'xy':
            for attribute in ['label', 'scale', 'lim']:
                setattr(self,
                        f'{axis}{attribute}',
                        getattr(self.plottable[axis],
                                attribute))

    def get_sizes(self):
        '''
        The sizes of the bubbles.

        Returns
        -------
        s : an input for plt.scatter
            Either a single scalar, or an array with variable
            sizes for each bubble according to some quantity.
        '''

        # if desired, set variable sizes
        if isinstance(self.plottable['size'], PlottableAxis):
            # get the raw values for the sizes
            x = self.plottable['size'].value()
            # calculate the normalized size
            size = default_size*x/self.size_normalization
        # otherwise, set a single size
        else:
            # get default, first from pop and then from panel
            size = self.pop.plotkw.get('s', self.plottable['size'])

        # return a valid input to plt.scatter(s=...)
        return size

    def get_colors(self):
        '''
        The colors of the bubbles.

        Returns
        -------
        c : an input for plt.scatter
            Either a single color, or an array with variable
            colors for each bubble according to some quantity.
        '''

        # should we ignore any variable color instructions?
        if self.pop.ink == False:
            color = self.pop.color
        # should we use a variable color?
        elif isinstance(self.plottable['color'], PlottableAxis):
            # get the raw values to go into the color
            x = self.plottable['color'].value()

            # FIXME - make sure to check vmin/vmax are valid
            #if (self.vmin is None) or (self.vmax is None):
            #    raise AtlasError(f'''
            #    It looks like you're trying to use
            #    {self.plottable['color']} to set variable
            #    colors for bubbles. To do so, please make
            #    sure it has finite values defined for its
            #    .vmin and .vmax attributes.
            #    ''')

            # make sure we have *some* normalizer defined
            f = plt.matplotlib.colors.Normalize
            self.color_normalization = (self.color_normalization
                or f(vmin=self.vmin, vmax=self.vmax))

            normalized = self.color_normalization(x)
            color = self.cmap(normalized)
        # finally, should we just use a default color?
        else:
            # get default, first from pop and then from panel
            color = self.pop.color
            if color is None:
                color = self.plottable['color']

        # return a valid input to any one of the following:
        #   plt.scatter(c=...)
        #   plt.scatter(edgecolors=...)
        #   plt.scatter(facecolors=...)
        return color

    def kw(self, key=None, **kwargs):
        '''
        Do a little decision-making about the plotting keyword
        arguments, pulling defaults from each population where
        needed.

        Parameter
        ---------
        key : str
            The population for which we should pull keywords.
            If None, go with the current population.
        **kwargs : dict
            All other keywords will be directed toward
            overwriting individual population defaults.
        '''

        # identify the population we're working with
        if key is None:
            key = self.key
        else:
            self.point_at(key)

        # define some default keywords, which can be over-written
        default = dict(s=self.get_sizes(),
                       marker=self.pop.marker,
                       linewidth=self.pop.linewidth,
                       alpha=self.pop.alpha,
                       zorder=self.pop.zorder,
                       label=self.pop.label)

        # sort out whether faces and/or edges should get color
        c=self.get_colors()
        if self.pop.filled:
            default['facecolors'] = c
        else:
            default['facecolors'] = 'none'
        if self.pop.outlined:
            default['edgecolors'] = c
        else:
            default['edgecolors'] = 'none'

        # if any other keywords are provided, overwrite these defaults
        for k, v in kwargs.items():
            default[k] = v

        return default

    def plot(self, key, ax=None, labelkw={}, **kwargs):
        '''
        Add the points for a particular population to this panel.

        Parameters
        ----------
        key : str
            The population (as an item in the self.pops dictionary) to add.
        ax :
            Into what ax should we place this plot?
            If None, use default.
        labelkw : dict
            Keywords for labeling the planet names.
        **kwargs : dict
            Any extra keywords will be passed on to `scatter`
        '''

        # focus attention on that population
        self.point_at(key)

        # make sure we're plotting into the appropriate axes
        try:
            plt.sca(self.ax)
        except AttributeError:
            self.setup(ax=ax)

        # add the scattered points
        self.scattered[key] = self.ax.scatter(self.x, self.y, **self.kw(key,**kwargs))

        # set the scales, limits, labels
        self.finish_plot(labelkw=labelkw)
