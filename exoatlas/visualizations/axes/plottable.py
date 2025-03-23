from ...imports import *


class Plottable:
    scale = "log"
    lim = [None, None]
    size_normalization = 1

    """
    General class definition for a plottable axis,
    which can appear as either an xaxis or yaxis.

    The standard items to include in the definition
    of a plottable object are:

        source = Generally a string, indicating a
                 retrievable population attribute,
                 which will be used as the value
                 to plot. The values-getting can
                 also be over-written by defining
                 a .value() method, if the desired
                 values are more complicated than
                 a simple attribute/property (e.g.
                 if it depends on an input, like
                 wavelength).

        label  = Human-friendly string describing
                 this axis, which will appear as
                 its label in the plot.

        scale  = String ('linear', 'log', ?) to
                 indicate what scale should be
                 used for this axis.

        lim    = Tuple to indicate (lower, upper)
                 limits for the plot. These might
                 also be used to set the vmin and
                 vmax of a color map, if this
                 plottable is being used to color
                 points in a BubblePanel.

        size_normalization = (optional) What to multiply
                 the values by to convert them
                 into sizes for a scatter plot.

    """

    def __init__(self, panel=None, orientation=None, **kw):
        """
        Initialize a plottable axis, connecting it
        to some parent panel that will handle all
        of the population cycling and building.

        Parameters 
        ----------
        panel : Panel
            To which plotting panel is this plottable connected? 
        orientation : str, None
            'horizontal' or 'vertical' for x or y,
            which might be used to format axis labels 
        kw : dict 
            All other keywords will be passed to .get()
            when retrieving values or uncertainties. 
            For example, for source='teq', we might 
            consider passing 'albedo' or 'f' keywords, 
            to therefore get the values and uncertainties 
            associated with a particular parameters. 
        """
        self.panel = panel
        self.orientation = orientation
        self.kw = kw

    def __call__(self, panel=None, orientation=None, **kw):
        """
        As a backup, in case we mistake a Plottable 
        object that has already been initialized, 
        this will produce a new independent instance. 
        Ideally this shouldn't happen, but since the 
        interface allows a few different ways to 
        specify a plottable, we want to have this catch.
        """
        new_instance = copy.deepcopy(self)
        new_instance.panel = panel
        new_instance.orientation = orientation

        return new_instance

    def __repr__(self):
        '''
        How should this plottable axis be represented as a string? 
        '''
        return f"🧮|{self.label}".replace("\n", " ")

    def value(self):
        """
        Extract the values for this plottable axis.

        This extracts values for population that's
        currently in focus for a plotting panel. 
        Generally, plotting panels might loop 
        through multiple populations to include
        each on the plot. 

        By default, this is done by using the 
        string in `self.source` to pull results 
        from a `Population` method. Any keywords 
        stored in `self.kw` will also be passed 
        to the call for that method.

        Write over this function in order to make
        more complicated function calls, if necessary.

        Returns
        -------
        value : Quantity
            The values for the currently in-focus population,
            as an astropy Quantity array. 
        """
        return self.panel.pop.get(self.source, **self.kw)

    def uncertainty_lowerupper(self):
        """
        Extract the upper and lower uncertainties. 

        This extracts values for population that's
        currently in focus for a plotting panel. 
        Generally, plotting panels might loop 
        through multiple populations to include
        each on the plot. 

        By default, this is done by using the 
        string in `self.source` to pull uncertainties 
        from a `Population` method. Any keywords 
        stored in `self.kw` will also be passed 
        to the call for that method.


        You might want to overwrite this method for 
        quantities where the uncertainty returned  
        by the Population's normal uncertainty 
        propagation is, for whatever reason, not 
        exactly what you want. 

        Returns 
        -------
        lower : np.array
            The magnitude of the lower uncertainties (x_{-lower}^{+upper}), 
            for the currently in-focus population, 
            as an astropy Quantity array. 
        upper : np.array
            The magnitude of the upper uncertainties (x_{-lower}^{+upper}),
            for the currently in-focus population, 
            as an astropy Quantity array. 
        """
        lower, upper = self.panel.pop.get_uncertainty_lowerupper(self.source, **self.kw)
        return lower, upper

        # FIXME!
        # this old code might not be necessary with the new uncertainty framework
        # it can probably be deleted after a bit more testing
        '''
        try:
            ul = self.panel.pop.get_uncertainty_lowerupper(self.source, **self.kw)
            return ul
        except (AtlasError, AttributeError, KeyError):
            sigma = self.panel.pop.get_uncertainty(self.source, **self.kw)
            return sigma, sigma
        '''


def clean_plottable(initial):
    """
    Make sure the axis initializer is a Plottable class.

    We allow fairly flexible calls to plotting panels,
    where it's possible that a Plottable might be 
    requested using a variety of different formats. 
    This wrapper tries to make sure we get to what we 
    need, a class definition, no matter what.

    Parameters
    ----------
    initial : Plottable, Plottable class, str, None
        One of a few ways to flexibly define a Plottable.
        Here's how the different options will be interpreted:

        Plottable = 
            If already an instance, the Plottable object itself 
            will be used for plotting. This would be common for 
            something like `StellarBrightness(wavelength=5*u.micron)`,
            where keywords are necessary for what's being plotted.
        Plottable class = 
            If a class definition, the Plottable will be 
            created using the defaults for that class. 
        str = 
            If a string, a Plottable will be created using 
            the method with that name from the populations. 
            Keywords should (CHECKME?!?!!?) still work, 
            but it will be more difficult to set bespoke 
            axis labels and/or plotting defaults. 
        None = 
            The parent panel probably has a Plottable 
            connected to a particular axis; None just means
            not to overwrite that default. 
    Returns
    -------
    plottable : (various)
        Either a Plottable, or something that can be interpreted 
        to create plottable when assigned data to plot in 
        a visualization Panel.
    """

    if isinstance(initial, Plottable):
        # pass through an actual Plottable object
        return initial
    elif type(initial) is type:
        # pass through an actual Plottable definition
        return initial
    elif type(initial) is str:
        # create a temporary Plottable from this string
        class GenericPlottable(Plottable):
            source = initial
            label = initial
            scale = "log"
            lim = [None, None]
        return GenericPlottable
    elif initial is None:
        # pass through None so panel can use its own axis
        return None
    else:
        # complain otherwise
        raise ValueError(
            f"""
        It' not clear how to turn {initial} into
        a definition for a Plottable.
        """
        )
