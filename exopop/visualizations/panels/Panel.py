# general class for plotting exoplanet populations
# all other panels derive from this one

from exopop.imports import *

# set the aspect ratios
aspect = 768/1024.0
figwidth = 7

class Panel(Talker):
    # define some defaults
    title = ''
    xlabel = '?'
    ylabel = '?'
    xscale = 'linear'
    yscale = 'linear'
    xlim = [None, None]
    ylim = [None, None]

    def __init__(self, pops={}):
        '''
        Initialize a plotting panel.

        Parameters
        ----------
        pops : dict
            A dictionary of populations, with keys as labels values as
            initialized populations. For example,

                pops = {'solarsystem':SolarSystem(),
                        'transiting':TransitingExoplanets()}
        '''

        # store the populations inside this panel
        self.pops = pops

        # dictionaries to store access points to items that appear on the plot
        self.scattered = {}
        self.labeled = {}

    def __repr__(self):
        '''
        How to represent this panel as a string?
        '''
        return f'<{self.__class__.__name__} panel | {len(self.pops)} populations>'

    @property
    def x(self):
        raise NotImplementedError(f"""
        The 'x' quantity hasn't been defined for
        {self}.
        """)

    @property
    def y(self):
        raise NotImplementedError(f"""
        The 'y' quantity hasn't been defined for
        {self}.
        """)

    def setup(self, ax=None):
        '''
        Set up this panel.

        Parameters
        ----------
        ax :
            The axes into which this panel should be placed.
        '''

        # what's the aspect ratio of this plot
        self.aspect = aspect

        # make sure we're pointing at the axes for this panel
        if ax is None:
            # if need be, create a new axes for this panel
            plt.figure(self.title, figsize=(figwidth, figwidth*self.aspect), dpi=300)
            gs = plt.matplotlib.gridspec.GridSpec(1,1,wspace=0,hspace=0)
            self.ax = plt.subplot(gs[0])
        else:
            # if an ax is provided, point at that
            self.ax = ax
            plt.sca(self.ax)

    def point_at(self, key):
        '''
        Point this panel at a particular population. This will be called
        when building up multiple populations on the same plot.

        Parameters
        ----------
        key : str
            The key of a particular population.
        '''

        # set the current population
        self.pop = self.pops[key]

        # record the name of the current population
        self.key = key


    def fileprefix(self, key):
        '''
        Define a fileprefix to use for saving this plot.
        '''
        return self.title.replace(' ','') + '_' + key.title()

    def build(self, **kw):
        '''
        Build up this panel by plotting every population into it.
        '''

        # loop over all the populations
        for key in self.pops.keys():

            # plot each population, passing plotting keywords to it
            self.plot(key, **kw)

    def label_planets(self, before='\n', after='', restrictlimits=False, **kwargs):
        '''
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
        '''

        # make sure we're set to the current axes
        plt.sca(self.ax)

        # loop over the elements in the population
        for i in range(len(self.x)):

            # pull out the positions and the name
            x, y, name = self.x[i], self.y[i], self.pop.name[i]

            # skip over the planets that aren't within limits
            if restrictlimits:
                try:
                    if x < np.min(self.xlim) or x > np.max(self.xlim):
                        if y < np.min(self.ylim) or y > np.max(self.ylim):
                            continue
                except TypeError: # (are there other errors we need to add?)
                    pass

            # define some defaults for the text
            textkw = dict(ha='center', va='top',
                          fontsize=6,
                          color=self.pop.color,
                          alpha=self.pop.alpha)

            # think this is just as Python 3 thing
            textkw.update(**kwargs)

            # store the text plot, so it can be modified
            try:
                assert(np.isfinite(x*y))
                self.labeled[name] = plt.text(x, y,
                                              before + name + after,
                                              **textkw)
            except AssertionError: # (do we need to add other errors?)
                pass

    def finish_plot(self, labelkw={}):
        '''
        Do some general things to finish up all panel plots.
        '''
        # implement the default scales, limiets, labels for the axes
        self.ax.set_xscale(self.xscale)
        self.ax.set_yscale(self.yscale)
        self.ax.set_ylim(*self.ylim)
        self.ax.set_xlim(*self.xlim)
        self.ax.set_xlabel(self.xlabel)
        self.ax.set_ylabel(self.ylabel)

        # if requested, label the individual planets in the plot
        if self.pop.label_planets:
            self.label_planets(**labelkw)
