from ...populations.calculations.shoreline import *
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
