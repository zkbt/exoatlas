from .GridGallery import *
from ..plottables import *
from ..animation import *


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


class SliceAnimatedGallery(SliceGridGallery):
    """
    A gallery to flip through multiple slices as an animation.
    """

    def __init__(self, map_to_slice, N=4, **kw):
        """
        Initialize an animation showing different slices of the same 3D space.

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
        self.setup_maps(map_to_slice=map_to_slice, N=N, **kw)

    def setup_maps(self, map_to_slice, N=4, figsize=(6, 4), **kw):

        # create the subplots (just one, in this case)
        nrows, ncols = 1, 1
        self.fi, self.ax = self.create_subplots(
            nrows=nrows, ncols=ncols, figsize=figsize, **kw
        )
        self.ax = np.atleast_2d(self.ax).reshape((nrows, ncols))

        overall_limits = map_to_slice.plottable["slice"].lim
        edges = np.linspace(*overall_limits, N + 1)

        # create the maps and attach all to the same axes
        self.maps = {}
        for i in range(N):
            m = copy.deepcopy(map_to_slice)
            m.plottable["slice"].lim = edges[i : i + 2]
            m.ax = self.ax[0, 0]
            self.maps[i] = m

    def animate(
        self,
        *args,
        filename="exoatlas-movie.mp4",
        fps=2,
        metadata={},
        codec=None,
        bitrate=None,
        bounce=True,
        pause=True,
        **kwargs,
    ):

        # decide animation writer based on the filename
        writer = get_animation_writer(
            filename, fps=fps, codec=codec, bitrate=bitrate, metadata=metadata
        )

        # set up the animation writer
        with writer.saving(self.fi, filename, self.fi.get_dpi()):

            keys = list(self.maps.keys())
            if bounce:
                if pause:
                    what_to_loop_over = keys + keys[::-1]
                else:
                    what_to_loop_over = keys + keys[1:-1][::-1]
            else:
                what_to_loop_over = keys
            print(what_to_loop_over)
            # loop over slices
            for i in tqdm(what_to_loop_over):

                m = self.maps[i]

                # entirely clear the current axes
                m.ax.cla()
                self.ax[0, 0].cla()
                plt.gca().cla()

                # build the plot, including populations
                m.build(*args, **kwargs)

                # refine the plot
                m.refine(**kwargs)

                # save this snapshot to a movie frame
                writer.grab_frame()
