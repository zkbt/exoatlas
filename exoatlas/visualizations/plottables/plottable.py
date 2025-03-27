from ...imports import *
from ...populations import Population
from matplotlib.colors import Normalize, LogNorm


class Plottable:
    source = None
    label = None
    scale = "log"
    lim = [None, None]
    unit = None
    """
    General class definition for a plottable axis,
    which can appear as either an xaxis or yaxis.

    The standard items to include in the definition
    of a plottable object are:


    """

    def __init__(self, source=None, label=None, scale=None, lim=None, unit=None, **kw):
        """
        Initialize a Plottable quantity.

        This `Plottable` serves as a helper for plotting
        `exoatlas` quantities. It provides functions to
        return values and uncertainties for a `Population`,
        but it also stores important information about
        labels, scales, and limits that are needed for
        good automatic visualizations (but are not needed
        for simply doing calculations).

        Often a `Plottable` will be connected to some kind
        of a visualization `Map`, which will handle rendering
        data on plots and cycling through populations.

        Parameters
        ----------
        source : str
            Generally a string, indicating a
            retrievable population attribute,
            which will be used as the value
            to plot. The values-getting can
            also be over-written by defining
            a .value() method, if the desired
            values are more complicated than
            a simple attribute/property (e.g.
            if it depends on an input, like
            wavelength).
        label : str
            Human-friendly string describing
            this axis, which will appear as
            its label in the plot.
        scale : str
            String ('linear', 'log', ?) to
            indicate what scale should be
            used for this axis.
        lim : tuple
            Tuple to indicate (lower, upper)
            limits for the plot. These might
            also be used to set the vmin and
            vmax of a color map, if this
            plottable is being used to color
            points in a BubbleMap. `lim`
            should have no units attached.
            For compariing different populations,
            it's probably a very good idea to
            explicitly set the limits and not
            leave them as None; otherwise
            individual populations might normalize
            themselves separately.
        unit : astropy.units.Unit
            The astropy unit to be used for plotting.
            It must be something to which the
            quantity can be converted.
        kw : dict
            All other keywords will be passed to .get()
            when retrieving values or uncertainties.
            For example, for source='teq', we might
            consider passing 'albedo' or 'f' keywords,
            to therefore get the values and uncertainties
            associated with a particular parameters.
        """

        # set properties to be defaults, updated by keywords
        self.source = source or self.source
        self.label = label or self.label
        self.scale = scale or self.scale
        self.lim = lim or self.lim
        self.unit = unit or self.unit

        # set keywords to be defaults for source, updated by **kw
        self.kw = self._get_default_keyword_arguments() | kw

        # after .source and .kw have been defined, update label, if necessary
        self._update_label()

    def _update_label(self):
        """
        Update the label.

        Sometimes, we might not know everything
        we need to write an appropriate label
        until after the `Plottable` has been
        created. For example, if a quantity is
        wavelength-dependent, we need to make
        an instance of the `Plottable` with a
        real wavelength, before we know how to
        set the wavelength label. If so,
        this method should be overwritten
        when defining a new `Plottable` class.

        When called, it should change
        `self.label`, and it can count on
        `self.source` and `self.kw` already
        being defined.
        """
        pass

    def _get_default_keyword_arguments(self):
        """
        Figure out default science keyword arguments.

        For the `source` method that returns the quantities
        to be plotted, figure out what the default keyword
        arguments are, so if this `Plottable` gets created
        without custom keyword values, we can still
        make use of the keyword inputs.
        """
        try:
            # get method for calculating quantity
            f = getattr(Population, self.source)
        except AttributeError:
            # Table columns that don't have explicit
            # method definitions probably won't show
            # up in the Population class until after
            # an instance is created. That's OK, because
            # those won't have custom keywords we need
            # care about.
            return {}

        # figure out default keywords
        kw = {
            k: v.default
            for k, v in inspect.signature(f).parameters.items()
            if v.default is not inspect.Parameter.empty
        }
        try:
            kw.pop("distribution")
        except KeyError:
            pass
        return kw

    def __repr__(self):
        """
        How should this plottable axis be represented as a string?
        """
        return f"📏 {self.source or self.label or '?'}".replace("\n", " ")

    def _convert_unit(self, x):
        """
        Convert quantity into desired unit,
        unless units aren't specified.

        Parameter
        ---------
        x : quantity
            The quantity to convert to `self.unit`
        """
        if self.unit is None:
            return x
        else:
            try:
                return x.to(self.unit)
            except AttributeError:
                return x

    def value(self, pop):
        """
        Extract the values for a population.

        Generally, plotting maps might loop
        through multiple populations to include
        each on the plot.

        By default, this is done by using the
        string in `self.source` to pull results
        from a `Population` method. Any keywords
        stored in `self.kw` will also be passed
        to the call for that method.

        Write over this function in order to make
        more complicated function calls, if necessary.

        Parameters
        ----------
        pop : Population
            The `exoatlas` population for which
            values will be retreived.

        Returns
        -------
        value : Quantity
            The value an astropy Quantity array.
        """
        value = pop.get(self.source, **self.kw)
        return self._convert_unit(value)

    def uncertainty_lowerupper(self, pop):
        """
        Extract the upper and lower uncertainties for a population.

        Generally, plotting maps might loop
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

        Parameters
        ----------
        pop : Population
            The `exoatlas` population for which
            values will be retreived.

        Returns
        -------
        lower : np.array
            The magnitude of the lower uncertainties (x_{-lower}^{+upper}),
            as an astropy Quantity array.
        upper : np.array
            The magnitude of the upper uncertainties (x_{-lower}^{+upper}),
            as an astropy Quantity array.
        """
        lower, upper = pop.get_uncertainty_lowerupper(self.source, **self.kw)
        return self._convert_unit(lower), self._convert_unit(upper)

    def uncertainty(self, pop):
        """
        Extract the uncertainties for a population.

        Generally, plotting maps might loop
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

        Parameters
        ----------
        pop : Population
            The `exoatlas` population for which
            values will be retreived.

        Returns
        -------
        sigma : np.array
            The uncertainties an astropy Quantity array.
        """
        sigma = pop.get_uncertainty(self.source, **self.kw)
        return self._convert_unit(sigma)

    def normalized_value(self, pop):
        """
        Extract the values for a population,
        normalized to fall between 0 and 1.

        This normalizes the values for a population
        according to scale

        Write over this function in order to make
        more complicated function calls, if necessary.

        Parameters
        ----------
        pop : Population
            The `exoatlas` population for which
            values will be retreived.

        Returns
        -------
        normalized_value : Quantity
            The value, normalized
        """

        # get unitless array of values
        value = remove_unit(self._convert_unit(pop.get(self.source, **self.kw)))

        # set non-nan limits if none are provided
        vmin, vmax = self.lim
        if vmin == None:
            vmin = np.nanmin(value)
        if vmax == None:
            vmax = np.nanmax(value)

        # call different scaling functions
        reversed = vmin > vmax
        if reversed:
            kw = dict(vmin=vmax, vmax=vmin)
        else:
            kw = dict(vmin=vmin, vmax=vmax)
        if self.scale == "linear":
            normalize = Normalize(**kw)
        elif self.scale == "log":
            normalize = LogNorm(**kw)
        if reversed:
            return 1 - normalize(value)
        else:
            return normalize(value)


def clean_plottable(initial, **kw):
    """
    Make sure something is actually a Plottable object.

    We allow fairly flexible calls to plotting maps,
    where it's possible that a Plottable might be
    requested using a variety of different formats.
    This wrapper tries to make sure we get to what we
    need, a class definition, no matter what.

    (This used to accept a string to initialize a
    basic )

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
        None =
            The parent map probably has a Plottable
            connected to a particular axis; None just means
            not to overwrite that default.
    **kw : dict
        Extra keywords that should be passed in when
        initializing `Plottable` from a class or a string.
    Returns
    -------
    plottable_or_not : (various)
        Either an instance of a Plottable object,
        or something else entirely.
    """

    if isinstance(initial, Plottable):
        # pass through an actual Plottable object
        return initial
    elif type(initial) == type:
        if issubclass(initial, Plottable):
            # create an instance from a Plottable class, using defaults
            return initial(**kw)
    else:
        # pass others through
        return initial
