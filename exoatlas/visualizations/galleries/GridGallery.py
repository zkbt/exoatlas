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

    def refine_maps(self):
        """
        Make small changes to Maps, after data are plotted.
        """
        for a in np.atleast_1d(self.ax[:, 1:]).flatten():
            a.set_ylabel("")
        for a in np.atleast_1d(self.ax[:-1, :]).flatten():
            a.set_xlabel("")

    def add_panel_letters(self, preset="inside",  **kw):
        """
        Add (a), (b), (c) labels to a group of axes.

        Parameters
        ----------
        ax : list, AxesSubplot
            The axes into which the labels should be drawn.
        preset : str
            A few presets for where to put the labels relative to
            upper left corner of each panel. Options are ['inside', 'above']
        kw : dict
            All addition keywords will be passed to `plt.text`,
            and they will overwrite defaults.
        """

        axes = [m.ax for m in self.maps.values()]


        textkw = dict(x=0.02, y=0.98, va="top", ha="left")
        if preset == "inside":
            textkw.update(x=0.02, y=0.98, va="top")
        elif preset == "outside":
            textkw.update(x=0, y=1.02, va="bottom")
        textkw.update(**kw)

        letters = "abcdefghijklmnopqrstuvwxyz"
        for i, a in enumerate(axes):
            textkw["s"] = f"({letters[i]})"
            textkw["transform"] = a.transAxes
            a.text(**textkw)
