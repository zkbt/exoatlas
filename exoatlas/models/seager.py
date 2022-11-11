from ..imports import *

__all__ = ["plot_seager", "plot_both_seager"]

# plot the mass radius models
def plot_seager(which="rock", **kw):
    """
    Plot a single model from Seager et al. (2007). These were
    originally point-clicked out of the GJ1214b discovery paper.

    Parameters
    ----------
    which : str
        Options are 'rock' or 'HHe'
    """

    assert (which.lower() == "rock") or (which.lower() == "hhe")

    file = resource_filename(__name__, f"data/seager/gj1214mr_{which}.txt")

    # load the values
    m, r = np.transpose(np.loadtxt(file))

    # fit a simply polynomial to smooth out the messiness
    fit = np.poly1d(np.polyfit(np.log(m), np.log(r), 5))

    # plot over a grid of masses
    mm = np.logspace(np.log10(min(m)), np.log10(max(m)), 1000)
    return plt.plot(mm, np.exp(fit(np.log(mm))), **kw)


def plot_both_seager(**kw):
    """
    Plot both the rock and the cold H/He models from Seager et al. (2007)
    """
    zkw = dict(linewidth=4, alpha=0.25, zorder=1.0)
    zkw.update(**kw)
    plot_seager("hhe", label="cold H/He", color="darkorange", **zkw)
    plot_seager("rock", label="rock", color="black", **zkw)
