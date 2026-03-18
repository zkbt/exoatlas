from .Gallery import *
from ..plottables import *


class GridGallery(Gallery):
    def __init__(
        self,
        rows=[Declination, Radius],
        cols=[RightAscension, Flux],
        mapsize=(5, 3),
        **kw,
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

    def refine_maps(self, **kw):
        """
        Make small changes to Maps, after data are plotted.
        """
        for a in np.atleast_1d(self.ax[:, 1:]).flatten():
            a.set_ylabel("")
        for a in np.atleast_1d(self.ax[:-1, :]).flatten():
            a.set_xlabel("")


class SliceGridGallery(GridGallery):
    """
    A grid gallery to show multiple slices side by side.
    """

    def __init__(self, map_to_slice, N=4, **kw):
        """
        Initialize a grid showing different slices of the same 3D space.

        This starts from a Map that with a defined 'slice' axis, and
        creates N slices along that axis to display in a grid of panels.

        Parameters
        ----------
        map_to_slice : Map
            One Map with a 'slice' axis defined, with finite limits.
        N : int
            The number of slices to produce,
            evenly splitting the limits
            of the 'map_to_slice'.
        **kw : dict
            All other keywords will be passed to 'setup_maps'.
        """
        # self.map_to_slice = map_to_slice
        # self.N = N
        self.setup_maps(map_to_slice=map_to_slice, N=N, **kw)

    def setup_maps(self, map_to_slice, N=4, figsize=(10, 3), mapsize=None, **kw):

        # create the subplots
        nrows, ncols = 1, N
        if mapsize is not None:
            figsize = (mapsize[0] * ncols, mapsize[1] * nrows)
        self.fi, self.ax = self.create_subplots(
            nrows=nrows, ncols=ncols, figsize=figsize, **kw
        )

        self.ax = np.atleast_2d(self.ax).reshape((nrows, ncols))

        overall_limits = map_to_slice.plottable["slice"].lim
        edges = np.linspace(*overall_limits, N + 1)

        # attach maps to each ax
        self.maps = {}
        for i, a in enumerate(self.ax.flatten()):
            m = copy.deepcopy(map_to_slice)
            m.plottable["slice"].lim = edges[i : i + 2]
            m.ax = a
            self.maps[i] = m

    # something clever to iterate functions over all the maps
    def __getattr__(self, x):
        """
        A sneaky wrapper to be able to call functions
        to apply to each slice of the gallery.

        For example, if `g` is a Gallery, then
        calling `g.refine()` will loop through
        all maps `m` within the gallery and run
        `m.refine()`. It's a handy way to make
        plotting refinements to all maps at once.
        """

        def f(*args, **kwargs):
            for i, m in self.maps.items():
                plt.sca(m.ax)
                getattr(m, x)(*args, **kwargs)

        return f
