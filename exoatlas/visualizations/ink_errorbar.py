from ..imports import *
from matplotlib.colors import Normalize

__all__ = ["ink_errorbar"]


def ink_errorbar(
    x,
    y,
    yerr=None,
    xerr=None,
    c=None,
    cmap=one2another("white", "black"),
    vmin=0,
    vmax=1,
    alpha=1.0,
    **kw
):
    """
    Draw errorbars, allowing each to have its own color.

    This is meant to be a drop-in replacement for plt.errorbar that
    can more easily give different colors to each point, using the c=
    and cmap= keywords. It's a little slow because it needs to draw the
    individual errorbars one-by-one, but for datasets smaller than a
    few thousand elements, it shouldn't be too dreadful.

    Parameters
    ----------

    x : array
        The x values
    y : array
        The y values
    yerr : array
        The y uncertainty value
    xerr : array
        The x uncertainty value
    c : array
        The array of numerical values to convert to colors.
        By default, these should fall between 0 and 1.
    cmap : colormap
        A colormap to convert c= values to colors.
    vmin, vmax : float, float
        The minimumu and maximum values to normalize the cmap.
    alpha : float or array
        The transparency by which *all* points will be multiplied.
    **kw : dict
        All other keywords will be passed to plt.errorbar.

    """

    i_sorted = np.argsort(c)

    # calculate Nx4 array of colors
    norm = Normalize(vmin=vmin, vmax=vmax)
    colors = cmap(norm(c))

    colors[:, -1] *= alpha

    # loop over data points
    for i in i_sorted:  # tqdm()
        # pick this data point's color
        # color = np.array(colors[i]).astype(np.float)
        # assert(len(color) == 4)

        # make more transparent, if alpha is <1
        # try:
        #    color[-1] *= alpha[i]
        # except TypeError:
        #    color[-1] *= alpha

        # rgba = (color[0], color[1], color[2], color[3])

        # set the color of every part of the point to be plotted
        kw["color"] = colors[i, :]
        kw["ecolor"] = colors[i, :]
        # kw['markeredgecolor'] = colors[i,:]
        # kw['markerfacecolor'] = colors[i,:]

        # try:
        # assert(len(zorder) > 1)
        # assert(len(zorder) == len(x))
        #    kw['zorder'] = zorder[i]
        # except (TypeError, KeyError):
        #    kw['zorder'] = zorder

        # reshape the errors as needed (for asymmetric errors)
        if yerr is not None:
            if len(yerr.shape) == 1:
                yerrtoplot = yerr[i]
            if len(yerr.shape) == 2:
                yerrtoplot = yerr[:, i].reshape(2, 1)
        else:
            yerrtoplot = None

        if xerr is not None:
            if len(xerr.shape) == 1:
                xerrtoplot = xerr[i]
            if len(xerr.shape) == 2:
                xerrtoplot = xerr[:, i].reshape(2, 1)
        else:
            xerrtoplot = None

        # actually draw the error bar for this one
        # print(kw['zorder'])
        assert np.isfinite(x[i])
        assert np.isfinite(y[i])
        assert np.isfinite(xerrtoplot).all()
        assert np.isfinite(yerrtoplot).all()

        plt.errorbar(x[i], y[i], xerr=xerrtoplot, yerr=yerrtoplot, **kw)
