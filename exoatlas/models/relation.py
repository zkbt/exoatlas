# imports from general utilities
import os
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate
import astropy.io.ascii
from importlib.resources import files
import string
import warnings


def remove_punctuation(s):
    """
    This helper function cleans a string of all punctuation.
    """
    translator = s.maketrans("", "", " " + string.punctuation)
    return s.translate(translator)


class Relation:
    """
    Base class for astrophysical relations,
    defining tools to read tables, define methods, etc...
    """

    def __init__(self, filename, **kwargs):
        """
        Initialize a Relation object.
        """

        # make sure the basic settings are included
        self.name = "Some General Relation"
        self.bibcode = "(no ADS bibcode defined)"

        # store the data filename
        self.filename = filename

        # load the data table
        self.load(**kwargs)

    def load(self, **kwargs):
        """
        Load the data table, and store it internally to this object.
        """

        # figure out the path to the data file, relative to this package
        # path = os.path.dirname(__file__) + '/'+ self.filename
        path = files(__name__) / self.filename

        # print("loading data from {0}".format(path))

        # load as an astropy table
        self.table = astropy.io.ascii.read(
            path, fill_values=[("...", np.nan)], **kwargs
        )

        # give an update
        # print("   ...success!")

    @property
    def possible(self):
        """
        What are the possible columns in this relation?
        """
        try:
            # if custom columns are defined, use those
            return self.columns
        except AttributeError:
            # otherwise, use all the columns in the table
            return self.table.keys()

    def limits(self, key):
        """
        Return the range of values spanned by
        a particular aspect of the relation.
        """
        return np.nanmin(self.table[key]), np.nanmax(self.table[key])

    def describe(self, key):
        """Describe one column in the table."""
        # print("{} = {}".format(key, self.descriptions[key]))

    def summarize(self):
        """Summarize all the possible columns in the table."""
        # print("")
        # print("The columns in {} are:".format(self.name))
        # for key in self.possible:
        # print("{:>20} = {}".format(key, self.descriptions.get(key, "???")))

    def tofrom(self, outkey, verbose=True):
        """
        Return a function that returns a function
        that will interpolate the values for given output (`outkey`),
        for some specified input (`inkey`).

        For example

            interpolator = Mamajek.tofrom('BCv')('Teff')

        will return a function called interpolator
        that takes values for Teff as input,
        and returns values for BCv.
        """

        # create a function that takes one input as a keyword arg
        def function(inkey):
            return self.interpolator(inkey=inkey, outkey=outkey)

        return function

    def interpolator(self, inkey=None, outkey=None):
        """
        Create an interpolator function, going from inkey to outkey.
        """
        # print("creating interpolator to convert {0} to {1}".format(inkey, outkey))
        try:
            x = self.table[inkey]
        except:
            warnings.warn(
                "it seems like the attempted input key {0} isn't valid".format(inkey)
            )
            return None
        try:
            y = self.table[outkey]
        except:
            warnings.warn(
                "it seems like the attempted output key {0} isn't valid".format(outkey)
            )
            return None

        # make sure to include only finite x-values
        ok = np.isfinite(x)
        return scipy.interpolate.interp1d(
            x[ok], y[ok], bounds_error=False, fill_value=np.nan
        )

    def plotone(self, inkey, outkey):
        """
        Plot one pair of columns.
        """
        try:
            self.ax.cla()
        except AttributeError:
            plt.figure("Relations Possible for " + self.__class__.__name__)
            self.ax = plt.subplot()

        if inkey == outkey:
            self.ax.text(
                0.5,
                0.5,
                inkey,
                fontsize=30,
                ha="center",
                va="center",
                transform=self.ax.transAxes,
            )
            self.ax.get_xaxis().set_visible(False)
            self.ax.get_yaxis().set_visible(False)

        else:
            interpolator = self.interpolator(inkey=inkey, outkey=outkey)
            x = np.linspace(np.nanmin(interpolator.x), np.nanmax(interpolator.x), 100)
            self.ax.plot(
                x, interpolator(x), alpha=0.3, color="sienna", linewidth=5, zorder=-1
            )

            ok = np.isfinite(self.table[inkey])
            realx, realy = np.array(self.table[inkey][ok]), np.array(
                self.table[outkey][ok]
            )
            i = np.argsort(realx)
            self.ax.plot(realx[i], realy[i], marker="o", alpha=0.5, color="black")
            self.ax.set_xlabel(inkey)
            self.ax.set_ylabel(outkey)

    def plot(self, keys=None):
        """
        Plot all possible pairs of columns.
        """
        if keys == None:
            keys = self.possible

        N = len(keys)
        size = 5

        fig = plt.figure(
            "Relations Possible for " + self.__class__.__name__,
            figsize=(N * size, N * size),
            dpi=100,
        )
        gs = plt.matplotlib.gridspec.GridSpec(N, N, figure=fig)

        for i in range(N):
            outkey = keys[i]
            for j in range(N):
                inkey = keys[j]
                self.ax = plt.subplot(gs[i, j])
                try:
                    self.plotone(inkey, outkey)
                    # print("plotted {0} to {1}".format(inkey, outkey))
                except (ValueError, TypeError):
                    pass
                    # print("failed to plot {0} to {1}".format(inkey, outkey))
                self.ax.set_title("{} -> {}".format(inkey, outkey))
                self.ax.set_xlabel("x = {}".format(inkey))
                self.ax.set_ylabel("f(x) = {}".format(outkey))
            # print("{0}/{1}".format(i + 1, len(self.possible)))

        outfile = "{}_{}.pdf".format(
            self.name, "+".join([remove_punctuation(k) for k in keys])
        )
        # ax = fig.add_subplot(111)
        # ax.set_xlabel('Input [x]')
        # ax.set_ylabel('Output [f(x)]')
        # print("saving summary figure to {}".format(outfile))
        plt.savefig(outfile)
        plt.draw()
