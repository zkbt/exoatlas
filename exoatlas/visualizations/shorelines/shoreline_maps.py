from ...calculations.shoreline import *
from ...visualizations.maps import ErrorMap, BubbleMap


class ShorelineErrorMap(ErrorMap, Shoreline):
    """
    Define a cosmic shoreline map, for looking at
    some slice of the 3D cosmic shoreline.

    It contains the powers both of an ErrorMap
    and of a Shoreline with model posteriors.
    """

    def __init__(self, **kw):
        super().__init__(**kw)
        Shoreline.__init__(self, **kw)

    def add_colorbar(self):
        """
        Add a colorbar showing the shoreline probability.

        (This can be called only if self.plot_shoreline_probability() has been too.)
        """

        # put colors on the colorbar
        colorbar = plt.colorbar(
            self._shoreline_probability_imshow,
            label="P(atmosphere)",
            norm=plt.matplotlib.colors.Normalize(vmin=0, vmax=1),
        )

        # set tick labels
        colorbar.set_ticks([0, 0.5, 1])

        # add lines from the contours, with dashed linestyle
        colorbar.add_lines(self._shoreline_probability_contours)
        for l in colorbar.lines:
            l.set_linestyles("--")


class ShorelineBubbleMap(BubbleMap, Shoreline):
    """
    Define a cosmic shoreline map, for looking at
    some slice of the 3D cosmic shoreline.

    It contains the powers both of an BubbleMap
    and of a Shoreline with model posteriors.
    """

    def __init__(self, **kw):
        super().__init__(**kw)
        Shoreline.__init__(self, **kw)
