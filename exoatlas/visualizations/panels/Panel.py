# general class for plotting exoplanet populations
# all other panels derive from this one

from ...imports import *
from ..axes.plottable import *
from ...populations import Population

# set the default aspect ratios
aspect = 768 / 1024.0
figwidth = 7


def clean_pops(initial):
    """
    Make sure the populations are stored as a dictionary,
    with reasonable labels for each as the keys.

    Parameters
    ----------
    initial : dict, list, Population
        A dictionary containing zero or more populations.
        Or, a list of one or more populations.

    Returns
    -------
    pops : dict
        A dictionary containing one or more populations.
    """

    if type(initial) == dict:
        return initial
    elif type(initial) == list:
        return {i: p for i, p in enumerate(initial)}
    else:
        # (otherwise, assume it's already a single population)
        return {initial.label: initial}


class Panel(Talker):
    # define some defaults
    title = None
    xaxis = None 
    yaxis = None 
    
    def __init__(self, xaxis=None, yaxis=None, label=None, **kw):
        """
        Initialize a plotting panel.

        Parameters
        ----------
        """

  

        # the 2-4 quantities being experssed in this plot
        self.plottable = {}

        # the matplotlib objects to access the plotted points themselves
        self.scattered = {}

        # the matplotlib objects to access text labels in the plot
        self.labeled = {}

        # the populations that will be plotted in this panel 
        self.populations = {}

        # set up the x and y axes
        self.plottable["x"] = clean_plottable(xaxis or self.xaxis, **kw)
        self.plottable["y"] = clean_plottable(yaxis or self.yaxis, **kw)

        # store a label for this panel
        self.label = label or f'{self.plottable["x"].source}'

        # apply axis labels, scales, limits appropriately
        for axis in "xy":
            for attribute in ["label", "scale", "lim"]:
                setattr(
                    self, f"{axis}{attribute}", getattr(self.plottable[axis], attribute)
                )

    def __repr__(self):
        """
        How should this plottable axis be represented as a string?
        """
        individual_plottable_strings = [
            f"{k:>6} = {v}" for k, v in self.plottable.items() if v is not None
        ]
        plottable_string = "\n".join(individual_plottable_strings)
        return f"🖼️ {self.__class__.__name__}\n{plottable_string}\n"

    def point_at(self, pop):
        """
        Point this panel at a particular population. This will be called
        when building up multiple populations on the same plot.

        Parameters
        ----------
        pop : Population, str, int
            The population to plot. This could be either an 
            actual Population object, or a key referring to 
            an element of the `self.populations` dictionary.
        """

        if isinstance(pop, Population):
            # if a Population, focus on it and make sure it is in dictionary
            needs_to_be_added_to_dictionary = True 
            for k, v in self.populations.items():
                # check whether populations match (= don't trust keys)
                if pop == v:
                    self.pop_key = k 
                    needs_to_be_added_to_dictionary = needs_to_be_added_to_dictionary 
            if needs_to_be_added_to_dictionary:
                self.pop_key = pop.label
                self.populations[self.pop_key] = pop

            self.pop = self.populations[self.pop_key]

        elif pop in self.populations:
            # if a key, focus on Population from dictionary
            self.pop = self.populations[pop]
            self.pop_key = pop 
        else:
            raise ValueError(f'''
            It's not clear how to interpret {pop}
            as a population at which we might point 
            {self}''')

    @property
    def x(self):
        return self.plottable["x"].value(self.pop)

    @property
    def y(self):
        return self.plottable["y"].value(self.pop)

    @property
    def x_lowerupper(self):
        return self.plottable["x"].uncertainty_lowerupper(self.pop)

    @property
    def y_lowerupper(self):
        return self.plottable["y"].uncertainty_lowerupper(self.pop)

    def setup_axes(self, ax=None):
        """
        Set up the axes for this Panel.

        Parameters
        ----------
        ax :
            The axes into which this panel should be placed.
            If no axes are provided, a new figure + axes will be created
        """

        # make sure we're pointing at the axes for this panel
        if ax is None:
            # if need be, create a new axes for this panel
            try:
                self.fi, self.ax
            except AttributeError:
                self.fi, self.ax = plt.subplots(figsize=(figwidth, figwidth * aspect), dpi=300, constrained_layout=True)
        else:
            # if an ax is provided, point at that
            self.ax = ax
            self.fi = ax.figure
        plt.sca(self.ax)

    def kw(self, pop=None, **kwargs):
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
        **kwargs : dict
            All other keywords will be directed toward
            overwriting individual population defaults.
        """

        # make sure we're pointing at the right population
        if pop != None:        
            self.point_at(pop)

        # define some default keywords, which can be over-written
        default = dict(
            alpha=self.pop.alpha,
            zorder=self.pop.zorder,
            label=self.pop.label,
        )

        # if any other keywords are provided, overwrite these defaults
        default.update(**kwargs)

        return default

    def plot(self, pop, ax=None, labelkw={}, **kwargs):
        """
        Add the points for a particular population to this panel.

        Parameters
        ----------
        pop : Population, str, int
            The population to plot. This could be either an 
            actual Population object, or a key referring to 
            an element of the `self.populations` dictionary.
        ax :
            Into what ax should we place this plot?
            If None, create a new plotting axes.
        labelkw : dict
            Keywords for labeling the planet names.
        **kwargs : dict
            Any extra keywords will be passed on to `scatter`
        """

        # focus attention on that population
        self.point_at(pop)

        # make sure we're plotting into the appropriate axes
        self.setup_axes(ax=ax)

        # add the scattered points
        these_scattered_points = self.ax.scatter(self.x, self.y, **self.kw(**kwargs))
        if self.pop_key in self.scattered:
            warnings.warn(f'''
            Key '{self.pop_key}' already exists. This might be fine for plotting, 
            but if you want to access any of the plotted elements to modify them 
            later, you might not be able to because they may have been overwritten.
            Might we please encourage you to give your populations unique labels 
            via the `population.label = "here's some neat label"`?
            '''
            )
        self.scattered[self.pop_key] = these_scattered_points

        # set the scales, limits, labels
        self.refine_axes()

        # add planet or hostname labels
        self.add_system_annotations(labelkw=labelkw)

    def refine_axes(self):
        """
        Update the basic axes properties after plotting.

        This applies the default scale, limits, and axis 
        labels to this Panel, after data have been plotted 
        in it. Technically, it only needs to be called after 
        all populations have been added to a plot, but 
        we call it at the end of each `.plot()` to make 
        sure it happens at least once.
        """
        # implement the default scales, limiets, labels for the axes
        self.ax.set_xscale(self.xscale)
        self.ax.set_yscale(self.yscale)
        self.ax.set_ylim(*self.ylim)
        self.ax.set_xlim(*self.xlim)
        self.ax.set_xlabel(self.xlabel)
        self.ax.set_ylabel(self.ylabel)

    def add_system_annotations(self, labelkw={}):
        '''
        Add text annotations for planet names.
        FIXME! 
        '''
        # if requested, label the individual planets in the plot
        if self.pop.label_planets:
            kw = dict(**labelkw)
            if getattr(self.pop, "zorder", None) is not None:
                kw["zorder"] = self.pop.zorder
            kw.update(**getattr(self.pop, "labelkw", {}))
            self.label_planets(self.pop, **kw)

    def build(self, pops={}, **kw):
        """
        Build up this panel by plotting every population into it.

        Parameters
        ----------
        pops : dict
            A dictionary of populations, with keys as labels values as
            initialized populations. For example,

                pops = {'solarsystem':SolarSystem(),
                        'transiting':TransitingExoplanets()}
        """

        # start with a dictionary of populations
        self.populations.update(clean_pops(pops))

        # loop over all the populations
        for key in self.populations.keys():

            # plot each population, passing plotting keywords to it
            self.plot(key, **kw)

    def label_planets(self, pop, before="\n", after="", restrictlimits=False, **kwargs):
        """
        Label the planets in whatever population we're pointed at.

        Parameters
        ----------
        pop : Population, str, int
            The population to plot. This could be either an 
            actual Population object, or a key referring to 
            an element of the `self.populations` dictionary.
        before : str
            Add the string at the start of each name.
        after : str
            Add the string at the end of each name.
        restrictlimits : bool
            Should we plot names only for those planets that fall within
            the panel's x and y limits?
        **kwargs : dict
            Any additional keywords will be passed to the `text` command.
        """

        # make sure we're pointing at this population
        self.point_at(pop)

        # make sure we're set to the current axes
        plt.sca(self.ax)

        # loop over the elements in the population
        for i in range(len(self.x)):
            # pull out the positions and the name
            x, y, name = self.x[i], self.y[i], self.pop.name()[i]

            # skip over the planets that aren't within limits
            if restrictlimits:
                try:
                    with np.errstate(invalid="ignore"):
                        if x < np.min(self.xlim) or x > np.max(self.xlim):
                            if y < np.min(self.ylim) or y > np.max(self.ylim):
                                continue
                except TypeError:  # (are there other errors we need to add?)
                    pass

            # define some defaults for the text
            textkw = dict(
                ha="center",
                va="top",
                fontsize=6,
                color=self.pop.color,
                alpha=self.pop.alpha,
                clip_on=True,
            )

            # think this is just as Python 3 thing
            textkw.update(**kwargs)

            # store the text plot, so it can be modified
            try:
                assert np.isfinite(x * y)
                self.labeled[name] = plt.text(x, y, before + name + after, **textkw)
            except AssertionError:  # (do we need to add other errors?)
                pass

    def label_hosts(
        self, pop, before="\n", after="", restrictlimits=False, once=False, **kwargs
    ):
        """
        Label the planet hosts in whatever population we're pointed at.

        Parameters
        ----------
        pop : Population, str, int
            The population to plot. This could be either an 
            actual Population object, or a key referring to 
            an element of the `self.populations` dictionary.
        before : str
            Add the string at the start of each name.
        after : str
            Add the string at the end of each name.
        restrictlimits : bool
            Should we plot names only for those planets that fall within
            the panel's x and y limits?
        **kwargs : dict
            Any additional keywords will be passed to the `text` command.
        """

        # make sure we're pointing at this population
        self.point_at(pop)

        # make sure we're set to the current axes
        plt.sca(self.ax)

        # loop over the elements in the population
        for i in range(len(self.x)):
            # pull out the positions and the name
            x, y, name = self.x[i], self.y[i], self.pop.hostname()[i]

            if once:
                if name in self.labeled:
                    continue

            # skip over the planets that aren't within limits
            if restrictlimits:
                try:
                    with np.errstate(invalid="ignore"):
                        if x < np.min(self.xlim) or x > np.max(self.xlim):
                            if y < np.min(self.ylim) or y > np.max(self.ylim):
                                continue
                except TypeError:  # (are there other errors we need to add?)
                    pass

            # define some defaults for the text
            textkw = dict(
                ha="center",
                va="top",
                fontsize=6,
                color=self.pop.color,
                alpha=self.pop.alpha,
                clip_on=True,
            )

            # think this is just as Python 3 thing
            textkw.update(**kwargs)

            # store the text plot, so it can be modified
            try:
                assert np.isfinite(x * y)
                self.labeled[name] = plt.text(x, y, before + name + after, **textkw)

            except AssertionError:  # (do we need to add other errors?)
                pass

    def connect_planets(self, pop, **kwargs):
        """
        Identify all the multiplanet systems,
        and draw lines connecting the planets
        within each individual system.

        Parameters
        ----------
        pop : Population, str, int
            The population to plot. This could be either an 
            actual Population object, or a key referring to 
            an element of the `self.populations` dictionary.
        **kwargs : dict
            Any keywords will overwrite the defaults
            going into the `.plot()` function call.
        """

        # make sure we're pointing at this population
        self.point_at(pop)

        # make sure we're set to the current axes
        plt.sca(self.ax)

        # define some defaults for the lines
        linekw = dict(color=self.pop.color, alpha=self.pop.alpha * 0.5, zorder=-100)

        # overwrite those defaults if kwargs are provided
        linekw.update(**kwargs)

        original_pop = self.pop
        for hostname in np.unique(self.pop.hostname()):
            friends = self.pop[hostname]
            self.pop = friends
            x, y = self.x, self.y
            i = np.argsort(x)
            plt.plot(x[i], y[i], **linekw)
            self.pop = original_pop


    def remove_xlabel(self):
        """
        Remove the xlabel and xticklabels from this panel.
        """
        plt.setp(self.ax.get_xticklabels(), visible=False)
        self.ax.set_xlabel("")

    def remove_ylabel(self):
        """
        Remove the ylabel and xticklabels from this panel.
        """
        plt.setp(self.ax.get_yticklabels(), visible=False)
        self.ax.set_ylabel("")

    def ticks_simplify_exponents(self, which="xy"):
        if "x" in which:
            self.ax.xaxis.set_major_formatter(FormatStrFormatter("%d"))
        if "y" in which:
            self.ax.yaxis.set_major_formatter(FormatStrFormatter("%d"))

    def ticks_enforce_multiple_oom(self, which="xy"):
        if "x" in which:
            loc = LogLocator(base=10, numticks=3)
            self.ax.xaxis.set_major_locator(loc)
        if "y" in which:
            loc = LogLocator(base=10, numticks=3)
            self.ax.yaxis.set_major_locator(loc)

    def add_legend(self, outside=True, frameon=False, **kw):
        legendkw = dict()
        legendkw["frameon"] = frameon
        if outside:
            legendkw["bbox_to_anchor"] = (1, 1)
            legendkw["loc"] = "upper left"
        legendkw.update(**kw)
        plt.legend(**legendkw)
