from imports import *
from Population import Population

class Custom(Population):
    def __init__(self, listofdictionaries=None, **kwargs):
        '''Initialize a population of KOIs, from the downloaded CSV.'''
        Population.__init__(self, label='New',  listofdictionaries=listofdictionaries, **kwargs)
        self.color='orangered'
        self.zorder = 10

    def loadFromScratch(self, **kwargs):
        pass

    def trimRaw(self, **kwargs):
        pass

    def loadRaw(self, **kwargs):
        pass

    def createStandard(self,  listofdictionaries=None, remake=False):
        '''Load a standardized population table, attempting...
            ...first from an .npy file (fast)
            ...then from a text file.'''

        d = listofdictionaries

        self.standard = astropy.table.Table(d)
        #    # and resave it as a numpy table (for faster loading next time)
        #    np.save(standard_numpy, self.standard)

        # add all of the table columns as attributes of the object
        self.propagate()
