"""
Define a subset of all exoplanets for only those known to transit. 
"""

from ...imports import *
from .exoplanets import *

__all__ = ["TransitingExoplanets"]


class TransitingExoplanets(Exoplanets):
    def __init__(self, **kw):
        Exoplanets.__init__(self, **kw)

        # set the label
        if "label" not in kw:
            self.label = "Transiting Exoplanets"

        # tidy up this population
        self._remove_nontransiting()

    def _remove_nontransiting(self):
        """
        Remove non-transiting planets from the population.

        Returns
        -------
        (modifies the input population in place, but returns nothing)
        """
        ok = self.detected_in_transit == 1
        self.standard = self.standard[ok]
