from ...imports import *

class PlottableAxis:
    scale = 'log'
    lim = [None, None]
    size_normalization = 1

    '''
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

    '''
    def __init__(self, panel=None,
                       orientation=None,
                       **kw):
        '''
        Initialize a plottable axis, connecting it
        to some parent panel that will handle all
        of the population cycling and building.
        '''
        self.panel = panel
        self.orientation = orientation
        self.kw = kw

    def __call__(self, panel=None, orientation=None, **kw):
        '''
        As a backup, in case we call something that
        looks like initializing this PlottableAxis,
        make sure that we connect to the appropriate
        panel.
        '''
        self.panel = panel
        self.orientation = orientation

        return self

    def __repr__(self):
        return f"<Plottable | {self.label}>".replace('\n', ' ')

    def value(self):
        '''
        Extract the values for this plottable axis.
        By default, this is done by pulling using
        the string in `source` to pull an attribute
        from a population.

        Write over this function in order to make
        more complicated function calls
        '''
        return getattr(self.panel.pop, self.source)

    def value_lowerupper(self):
        '''
        Extract the upper and lower uncertainties
        for this plottable axis. This function
        will likely need to be overwritten for
        any attribute that doesn't directly have
        an uncertainty defined inside the population.
        '''

        try:
            ul = self.panel.pop.uncertainty_lowerupper(self.source)
            return ul
        except AtlasError:
            sigma = self.value()*0.0
            return sigma, sigma

def clean_axis(initial):
    '''
    Make sure the axis initializer is a PlottableAxis
    class definition.

    Parameters
    ----------
    initial : PlottableAxis class definition, string

    Returns
    -------
    axis : PlottableAxis class definition
    '''

    if initial is None:
        # pass through None so panel can use its own axis
        return None
    elif type(initial) is type:
        # pass through an actual PlottableAxis definition
        return initial
    elif isinstance(initial, PlottableAxis):
        # pass through an actual PlottableAxis object
        return initial
    elif type(initial) is str:
        # create a temporary PlottableAxis from this string
        class GenericPlottableAxis(PlottableAxis):
            source = initial
            label = initial
            scale = 'log'
            lim = [None, None]
        return GenericPlottableAxis
    else:
        # complain otherwise
        raise ValueError(f'''
        It' not clear how to turn {initial} into
        a definition for a PlottableAxis.
        ''')
