from .imports import *
from .population import Population


class Custom(PredefinedPopulation):
    def __init__(self, listofdictionaries=None, **kwargs):
        """Initialize a population of KOIs, from the downloaded CSV."""
        Population.__init__(
            self, label="New", listofdictionaries=listofdictionaries, **kwargs
        )
        self.color = "darkorange"
        self.zorder = 10

    def loadFromScratch(self, **kwargs):
        pass

    def trim_raw(self, **kwargs):
        pass

    def load_raw(self, **kwargs):
        pass

    def create_standard(self, listofdictionaries=None, remake=False):
        """Load a standardized population table, attempting...
        ...first from an .npy file (fast)
        ...then from a text file."""

        d = listofdictionaries

        self.standard = Table(d)
        #    # and resave it as a numpy table (for faster loading next time)
        #    np.save(standard_numpy, self.standard)
