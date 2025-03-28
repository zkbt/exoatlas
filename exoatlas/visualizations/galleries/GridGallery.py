from .Gallery import *
from ..plottables import *


class GridGallery(Gallery):
    def __init__(
        self,
        rows=[Declination, Radius],
        cols=[RightAscension, Flux],
        mapsize=(5, 3),
        **kw
    ):
        self.mapsize = mapsize
        self.setup_maps(rows, cols, **kw)

    def setup_maps(self, rows=[], cols=[], map_type=BubbleMap, **kw):

        # create the subplots
        nrows, ncols = len(rows), len(cols)
        figsize = (self.mapsize[0] * len(cols), self.mapsize[1] * len(rows))
        self.fi, self.ax = self.create_subplots(
            nrows=nrows, ncols=ncols, figsize=figsize, **kw
        )

        self.ax = np.atleast_2d(self.ax).reshape((nrows, ncols))

        # attach maps to each ax
        self.maps = {}
        for i, r in enumerate(rows):
            for j, c in enumerate(cols):
                m = map_type(xaxis=c, yaxis=r)
                m.ax = self.ax[i, j]
                self.maps[m.label] = m

    def refine_maps(self):
        """
        Make small changes to Maps, after data are plotted.
        """
        for a in np.atleast_1d(self.ax[:, 1:]).flatten():
            a.set_ylabel("")
        for a in np.atleast_1d(self.ax[:-1, :]).flatten():
            a.set_xlabel("")
