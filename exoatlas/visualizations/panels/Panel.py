# general class for plotting exoplanet populations
# all other panels derive from this one

from ...imports import *
from ..axes.plottable import *

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
    title = ""

    def __init__(self, xaxis=None, yaxis=None, name="?", **kw):
        """
        Initialize a plotting panel.

        Parameters
        ----------
        """

        # store the name of this panel
        self.name = name

        # the 2-4 quantities being experssed in this plot
        self.plottable = {}

        # the matplotlib objects to access the plotted points themselves
        self.scattered = {}

        # the matplotlib objects to access text labels in the plot
        self.labeled = {}

        # set up the x and y axes
        self.plottable["x"] = clean_plottable(xaxis or self.xaxis, **kw)
        self.plottable["y"] = clean_plottable(yaxis or self.yaxis)

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

    def setup(self, ax=None):
        """
        Set up this panel.

        Parameters
        ----------
        ax :
            The axes into which this panel should be placed.
            If no axes are provided, a new figure will be created
        """

        # what's the aspect ratio of this plot
        self.aspect = aspect

        # make sure we're pointing at the axes for this panel
        if ax is None:
            # if need be, create a new axes for this panel
            fi, ax = plt.figure(self.title, figsize=(figwidth, figwidth * self.aspect), dpi=300)
            gs = plt.matplotlib.gridspec.GridSpec(1, 1, wspace=0, hspace=0)
            self.ax = plt.subplot(gs[0])
        else:
            # if an ax is provided, point at that
            self.ax = ax
            plt.sca(self.ax)

    def point_at(self, key):
        """
        Point this panel at a particular population. This will be called
        when building up multiple populations on the same plot.

        Parameters
        ----------
        key : str
            The key of a particular population.
        """

        # set the current population
        self.pop = self.pops[key]

        # record the name of the current population
        self.key = key

        return self

    def fileprefix(self, key):
        """
        Define a fileprefix to use for saving this plot.
        """
        return self.title.replace(" ", "") + "_" + key.title()

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
        self.pops = clean_pops(pops)

        # loop over all the populations
        for key in self.pops.keys():

            # plot each population, passing plotting keywords to it
            self.plot(key, **kw)

        return self

    def label_planets(self, before="\n", after="", restrictlimits=False, **kwargs):
        """
        Label the planets in whatever population we're pointed at.

        Parameters
        ----------
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
        self, before="\n", after="", restrictlimits=False, once=False, **kwargs
    ):
        """
        Label the planet hosts in whatever population we're pointed at.

        Parameters
        ----------
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

        # make sure we're set to the current axes
        plt.sca(self.ax)

        # loop over the elements in the population
        for i in range(len(self.x)):
            # pull out the positions and the name
            x, y, name = self.x[i], self.y[i], self.pop.hostname[i]

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

    def connect_planets(self, **kwargs):
        """
        Identify all the multiplanet systems,
        and draw lines connecting the planets
        within each individual system.

        Parameters
        ----------
        **kwargs : dict
            Any keywords will overwrite the defaults
            going into the `.plot()` function call.
        """

        # make sure we're set to the current axes
        plt.sca(self.ax)

        # define some defaults for the lines
        linekw = dict(color=self.pop.color, alpha=self.pop.alpha * 0.5, zorder=-100)

        # overwrite those defaults if kwargs are provided
        linekw.update(**kwargs)

        original_pop = self.pop
        for hostname in np.unique(self.pop.hostname):
            friends = self.pop[hostname]
            self.pop = friends
            x, y = self.x, self.y
            i = np.argsort(x)
            plt.plot(x[i], y[i], **linekw)
            self.pop = original_pop

    def finish_plot(self, labelkw={}):
        """
        Do some general things to finish up all panel plots.
        """
        # implement the default scales, limiets, labels for the axes
        self.ax.set_xscale(self.xscale)
        self.ax.set_yscale(self.yscale)
        self.ax.set_ylim(*self.ylim)
        self.ax.set_xlim(*self.xlim)
        self.ax.set_xlabel(self.xlabel)
        self.ax.set_ylabel(self.ylabel)

        # if requested, label the individual planets in the plot
        # for p in self.pops:
        #    self.point_at(p)
        if self.pop.label_planets:
            kw = dict(**labelkw)
            if getattr(self.pop, "zorder", None) is not None:
                kw["zorder"] = self.pop.zorder
            kw.update(**getattr(self.pop, "labelkw", {}))
            self.label_planets(**kw)

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
