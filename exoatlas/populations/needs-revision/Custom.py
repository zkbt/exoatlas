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

    def download_raw_data(self, **kwargs):
        pass

    def create_standardardized(self, listofdictionaries=None, remake=False):
        """Load a standardized population table, attempting...
        ...first from an .npy file (fast)
        ...then from a text file."""

        d = listofdictionaries

        self.standard = QTable(d)
        #    # and resave it as a numpy table (for faster loading next time)
        #    np.save(standard_numpy, self.standard)
