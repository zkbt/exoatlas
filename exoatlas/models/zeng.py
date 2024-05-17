from ..imports import *

__all__ = ["plot_zeng", "plot_three_zeng"]


# plot the mass radius models
def plot_zeng(which="75pSi_25pFe", **kw):
    """
    Plot a single mass-radius model from Zeng & Sasselov (2013).

    Parameters
    ----------
    which : str
        A string to indicate which model to plot.

    """

    file = os.path.join(code_directory, f"models/data/zengandsasselov/MR_{which}.csv")

    # load the values
    m, r = np.transpose(np.loadtxt(file, delimiter=","))

    # plot the values
    return plt.plot(m, r, **kw)


def plot_three_zeng(**kw):
    """
    Plot both the rock and the cold H/He models from Seager et al. (2007)
    """
    zkw = dict(linewidth=4, alpha=0.25, zorder=1.0)
    zkw.update(**kw)
    plot_zeng("75pSi_25pFe", label="75% Si, 25% Fe", color="sienna", **zkw)
    plot_zeng("50pSi_50pFe", label="50% Si, 25% 50", color="brown", **zkw)
    plot_zeng("25pSi_75pFe", label="25% Si, 75% Fe", color="maroon", **zkw)
