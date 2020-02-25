from .Panel import *
from ..ink_errorbar import *

__all__ = ['ErrorPanel']

def remove_unit(x):
    if type(x) == u.quantity.Quantity:
        return x.value
    else:
        return x

class ErrorPanel(Panel):
    '''
    Error is a general wrapper for making scatter plots
    where planets are represented with 2D error bars, with
    their intensity scaled to some overall visual weight.
    '''

    #@property
    #def x_lowerupper(self):
    #    raise NotImplementedError(f"""
    #    The 'x_lowerupper' quantity hasn't been defined for
    #    {self}.
    #    """)

    #@property
    #def y_lowerupper(self):
    #    raise NotImplementedError(f"""
    #    The 'x_lowerupper' quantity hasn't been defined for
    #    {self}.
    #    """)

    @property
    def x_unc(self):
        l, u = self.x_lowerupper
        #return np.sqrt(l*u)
        return 0.5*(l + u)

    @property
    def y_unc(self):
        l, u = self.y_lowerupper
        #return np.sqrt(l*u)
        return 0.5*(l + u)

    def intensity(self):
        '''
        What visual intensity should each datapoint have?

        By default, this will be the product of the fractional uncertainties
        on both x and y.

        Returns
        -------
        intensity : np.array
            The intensity of the points for the current population,
            as a numeric value between 0 and 1.
        '''

        dlnx = self.x_unc/self.x
        dlny = self.y_unc/self.y

        # things with bigger errors should have lower weight
        weight = (1 - np.sqrt(dlnx**2  + dlny**2)*1.5)

        # clip the weights above 1 (shouldn't exist?) or below zero
        # clipped = np.minimum(np.maximum(weight, 0), 1)
        return remove_unit(weight)

    def plot(self, key, ax=None, labelkw={}, **kw):
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
        **kw : dict
            Any extra keywords will be passed on to `errorbar`
        '''

        # focus attention on that population
        self.point_at(key)

        # make sure we're plotting into the appropriate axes
        try:
            plt.sca(self.ax)
        except AttributeError:
            self.setup(ax=ax)

        # define the data we're trying to plot
        x = remove_unit(self.x)
        y = remove_unit(self.y)

        # if the entire population is exact (Solar System), don't include errors
        exact = self.pop.exact #(self.x_unc == 0).all() and (self.y_unc == 0).all()
        if exact:
            plotkw = dict(color=self.pop.color,
                          edgecolor=self.pop.color,
                          **kw)
            plotkw['alpha'] = 1
            plotkw['zorder'] = 1e9
            self.scattered[key] = plt.scatter(x, y, **plotkw)
        else:

            # define the error bars to be plotting
            xl, xu = self.x_lowerupper
            x_unc = remove_unit(np.vstack([xl, xu])) # the vstack strips the units?

            yl, yu = self.y_lowerupper
            y_unc = remove_unit(np.vstack([yl, yu])) # the vstack strips the units?

            #y_unc = np.empty((2, len(xl)))*yl.unit
            #y_unc[0, :] = yl
            #y_unc[1, :] = yu

            width=1
            kw = dict(#marker='o',
                      linewidth=0,
                      elinewidth=width,
                      alpha=1.0,
                      #capthick=width,
                      #capsize=2,
                      #markersize=3)
                      #color=self.pop.color,
                      #markeredgecolor=self.pop.color,
                      )

            # pull out the default color of this population
            color = self.pop.color


            # define an Nx4 array of RGBA colors for the N points
            weights = self.intensity()

            # remove everything with infinite errorbars
            ok = (np.isfinite(xl) &
                  np.isfinite(xu) &
                  np.isfinite(x)  &
                  np.isfinite(yl) &
                  np.isfinite(yu) &
                  np.isfinite(y)  &
                  np.isfinite(weights))

            n_nouncertainty = sum(ok == False)
            self.speak(f'skipping {n_nouncertainty} planets that are missing data or uncertainties')

            # kludge to remove those that cross-zero

            with np.errstate(invalid='ignore'):
                ok *= (self.x - xl) > 0
                ok *= (self.y - yl) > 0

            n_consistentwithzero = sum(ok == False) - n_nouncertainty
            self.speak(f'skipping {n_consistentwithzero} planets that are consistent with zero')

            if (len(x) > 1) & (self.pop.plotkw.get('ink', True)):

                self.speak('plotting inked errorbars, this may take a while')
                # FIXME -- make the "invisible" color more flexible than white,
                # in case we're plotting on a dark background
                self.scattered[key] = ink_errorbar(x[ok], y[ok],
                                                   yerr=y_unc[:, ok], xerr=x_unc[:, ok],
                                                   c=weights[ok],
                                                   cmap=one2another(bottom='white',
                                                                    top=color,
                                                                    alphabottom=1.0,
                                                                    alphatop=1.0),
                                                   **kw)
            else:
                self.scattered[key] = self.ax.errorbar(x[ok], y[ok],
                                                       yerr=y_unc[:, ok], xerr=x_unc[:, ok],
                                                       **kw)

        # set the scales, limits, labels
        self.finish_plot(labelkw=labelkw)
