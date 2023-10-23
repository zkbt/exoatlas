"""
Tools to summarize a particular population.
"""

from ..imports import *
from .Population import attribute_columns, method_columns

__all__ = ["plot_histograms", "summarize_planet"]


def split_cols(pop):
    # find the columns that aren't strings
    quant, qual = [], []
    for x in attribute_columns:
        try:
            pop.standard[x][0] + "a"
            qual.append(x)
        except TypeError:
            quant.append(x)
    return quant, qual


def plot_histograms(pop):
    """
    Make histograms of all the necessary columns in the population.
    These can help see what's missing and what's there.
    """

    # pull out the quantitative columns
    quant_cols, string_cols = split_cols(pop)

    cols = []
    for x in string_cols:
        if len(np.unique(pop.standard[x])) < 50:
            cols.append(x)
    cols += quant_cols

    # create a grid of histograms
    ncols = 3
    nrows = np.ceil(len(cols) / ncols).astype(np.int)
    scale = 2
    fi, ax = plt.subplots(nrows, ncols, figsize=(ncols * 2 * scale, nrows * scale))

    for i, x in enumerate(cols):
        if x in quant_cols:
            good = np.isfinite(pop.standard[x])
        else:
            good = np.ones(len(pop.standard)).astype(bool)
        bad = good == False
        badfraction = sum(bad) / len(bad)

        plt.sca(ax.flatten()[i])
        plt.hist(pop.standard[x][good], color="black")
        plt.axvspan(*plt.xlim(), 0, badfraction, color="red", alpha=0.5, zorder=-1)
        plt.xlabel(x)
        plt.title(f"{x} lacks {sum(bad)}/{len(bad)} ({badfraction:.0%})")

    plt.tight_layout()


def summarize_planet(planets):
    """
    Quick tool to print out everything we need about a planet.
    """
    for p in planets:
        for k in attribute_columns:
            print(f"{k:>20} = {getattr(p, k)}")

        for k in method_columns:
            print(f"{k:>20} = {getattr(p, k)()}")
        print()
