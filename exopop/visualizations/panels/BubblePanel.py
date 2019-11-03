from .Panel import *

__all__ = ['BubblePanel']

class BubblePanel(Panel):
    '''
    BubblePanel is a general wrapper for making scatter plots
    where planets are represented as bubbles that can have
    informative sizes and/or colors.
    '''

    def size(self):
        '''
        The sizes of the bubbles.
        '''
        return None

    def color(self):
        '''
        The colors of the bubbles.
        '''
        return None

    def kw(self, key=None, **kwargs):
        '''
        Do a little decision-making about the plotting keyword arguments,
        pulling in the defaults from each population.

        Parameter
        ---------
        key : str
            The population for which we should pull keywords.
            If None, go with the current population.
        **kwargs : dict
            Any other keywords will go into overwriting defaults.
        '''

        # identify the population we're working with
        if key is None:
            key = self.key
        else:
            self.point_at(key)

        pop = self.pop

        # define some default keywords, which can be over-written
        default = dict(s=self.size(),
                       c=self.color(),
                       cmap='plasma',
                       marker='o',
                       linewidth=1,
                       alpha=pop.alpha,
                       zorder=pop.zorder,
                       label=pop.label)

        # if not using the cmap, just draw point edges
        if default['c'] is None:
            default['edgecolors']=pop.color,
            default['facecolors']='none'

        # if any other keywords are provided, overwrite these defaults
        for k, v in kwargs.items():
            default[k] = v

        # FIXME - we should make it a bit more intuitive for playing with
        # the colors of the points and cmaps

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
