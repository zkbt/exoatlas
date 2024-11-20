from .Panel import *
from .BubblePanel import BubblePanel
from ..ink_errorbar import *

__all__ = ["ErrorPanel"]


def remove_unit(x):
    """
    A simple wrapper to help make sure that we're dealing
    just with numerical values, instead of astropy Quantitys
    with units.
    """

    if type(x) == u.quantity.Quantity:
        return x.value
    else:
        return x


class ErrorPanel(BubblePanel):
    """
    Error is a general wrapper for making scatter plots
    where planets are represented with 2D error bars, with
    their intensity scaled to some overall visual weight.
    """

    @property
    def x_unc(self):
        l, u = self.x_lowerupper
        # return np.sqrt(l*u)
        return 0.5 * (l + u)

    @property
    def y_unc(self):
        l, u = self.y_lowerupper
        # return np.sqrt(l*u)
        return 0.5 * (l + u)

    def intensity(self, invisible_fraction=0.75, x_power=2, y_power=2):
        """
        What visual intensity should each datapoint have?

        By default, this will be set by the product of the
        fractional uncertainties on both x and y.

        Parameters
        ----------
        invisible_fraction : float
            The aproximate 2D fractional uncertainty at which
            a point will fade to invisible. For example, if
            set to 0.7, then points will disappear when
              (sigma_x/x)**2 + (sigma_y/y)**2 > 0.7**2

        x_power : float
            The power to which x is raised in the quantity
            used to define the visual weight. For example,
            if we want visual weight to scale with the
            fractional uncertainty on density (mass/radius**3)
            for x showing mass and y showing radius, we'd
            choose x_power=1 and y_power=3.

        y_power : float
            (see `x_power` above)

        Returns
        -------
        intensity : np.array
            The intensity of the points for the current
            population, as a numeric value between 0 and 1.
            By default,
        """

        dlnx = self.x_unc / self.x * np.abs(x_power)
        dlny = self.y_unc / self.y * np.abs(y_power)

        # things with bigger errors should have lower weight
        weight = 1 - np.sqrt(dlnx**2 + dlny**2) / invisible_fraction

        # clip the weights above 1 (shouldn't exist?) or below zero
        # clipped = np.minimum(np.maximum(weight, 0), 1)

        # return the visual weight
        return remove_unit(weight)

    def plot(self, key, ax=None, labelkw={}, **kw):
        """
        Add the points for a particular population to this panel.

        Parameters
        ----------
        key : str
            The item in the self.pops dictionary to add.
        ax :
            Into what ax should we place this plot?
            If None, use default.
        labelkw : dict
            Keywords for labeling the planet names.
        **kw : dict
            Any extra keywords will be passed on to `errorbar`
        """

        # focus attention on that population
        self.point_at(key)

        # do a bubble instead, if requested
        if getattr(self.pop, "bubble_anyway", False):
            BubblePanel.plot(self, key=key, ax=ax, labelkw=labelkw, **kw)
            return

        # make sure we're plotting into the appropriate axes
        try:
            plt.sca(self.ax)
        except AttributeError:
            self.setup(ax=ax)

        # define the data we're trying to plot
        x = remove_unit(self.x)
        y = remove_unit(self.y)

        # set the base color to use throughout
        default_color = plt.scatter([], []).get_facecolor()[0]
        color = self.pop.color or default_color
        marker = self.pop.marker or "o"

        # if the entire population is exact (e.g., Solar System),
        # then don't include any errors when plotting
        if self.pop.exact:
            # define plotting keywords without errorbars
            plotkw = dict(color=color, edgecolor=color, marker=marker, **kw)
            plotkw["alpha"] = 1
            plotkw["zorder"] = 1e9
            self.scattered[key] = plt.scatter(x, y, **plotkw)
            # FIXME, 5/25/20: I think BubblePanel is doing
            # something a little more clever with being able
            # to manage face and edge colors separately.
            # Perhaps we should set things up so that we might
            # inherit some of this skills here in ErrorPanel
        else:
            # define the error bars to be plotting
            xl, xu = self.x_lowerupper
            x_unc = remove_unit(np.vstack([xl, xu]))

            yl, yu = self.y_lowerupper
            y_unc = remove_unit(np.vstack([yl, yu]))

            width = 1
            kw = dict(  # marker='o',
                linewidth=0,
                elinewidth=width,
                alpha=1.0,
                # capthick=width,
                # capsize=2,
                # markersize=3)
                # color=self.pop.color,
                # markeredgecolor=self.pop.color,
            )

            # define an Nx4 array of RGBA colors for the N points
            weights = self.intensity()

            # remove everything with infinite errorbars
            ok = (
                np.isfinite(xl)
                & np.isfinite(xu)
                & np.isfinite(x)
                & np.isfinite(yl)
                & np.isfinite(yu)
                & np.isfinite(y)
                & np.isfinite(weights)
            )

            n_nouncertainty = sum(ok == False)
            # self._speak(
            #    f"skipping {n_nouncertainty} planets that are missing data or uncertainties"
            # )

            # kludge to remove those that cross zero
            with np.errstate(invalid="ignore"):
                ok *= (self.x - xl) > 0
                ok *= (self.y - yl) > 0
            # FIXME, 5/25/2020: This kludge always useful on
            # logarithmic axes (which are the only that have
            # been defined so far), but at some point we
            # might want to use linear axes too, where we might
            # not want to throw out values that might go
            # negative.

            n_consistentwithzero = sum(ok == False) - n_nouncertainty
            # self._speak(
            #    f"skipping {n_consistentwithzero} planets that are consistent with zero"
            # )

            if (len(x) > 1) & (self.pop._plotkw.get("ink", True)):
                # self.speak("plotting inked errorbars, this may take a while")
                # FIXME, 5/25/2020: We should make the
                # "invisible" color be something more flexible
                # than white, in case we're plotting on a dark
                # background. Remember, things look messy if
                # we use alpha to do the visual weighting for
                # these errorbars, because it introduces many
                # more intersecting lines.
                self.scattered[key] = ink_errorbar(
                    x[ok],
                    y[ok],
                    yerr=y_unc[:, ok],
                    xerr=x_unc[:, ok],
                    c=weights[ok],
                    cmap=one2another(
                        bottom="white", top=color, alphabottom=1.0, alphatop=1.0
                    ),
                    **kw,
                )
            else:
                self.scattered[key] = self.ax.errorbar(
                    x[ok],
                    y[ok],
                    yerr=y_unc[:, ok],
                    xerr=x_unc[:, ok],
                    color=self.pop.color,
                    **kw,
                )

        # set the scales, limits, labels
        self.finish_plot(labelkw=labelkw)
