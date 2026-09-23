"""
Define some commonly useful subsets of all transiting exoplanets,
giving them unique names and colors to simplify plotting
"""

from ...imports import *
from .transiting_exoplanets import *

__all__ = [
    "TransitingExoplanetsSubset",
    "Kepler",
    "NonKepler",
    "TESS",
    "NonTESS",
    "NonKeplerNonTESS",
    "Space",
    "Ground",
    "GoodMass",
    "BadMass",
]


class TransitingExoplanetsSubset(TransitingExoplanets):
    def __init__(self, label="Subset", **kw):
        TransitingExoplanets.__init__(self, **kw)

        # set the label
        self.label = label

        # trim to just the data we want
        self.table = self.table[self.to_include()]

        self._plotkw["color"] = kw.get("color", None)
        self._plotkw["c"] = kw.get("color", None)

    def to_include(self):
        raise NotImplementedError(
            "Please define `.to_include()` for this ExoplanetSubset!"
        )


class Kepler(TransitingExoplanetsSubset):
    def __init__(self, **kw):
        kw = dict( label="Kepler", color="royalblue", zorder=0) | kw

        TransitingExoplanetsSubset.__init__(
            self, **kw
        )

    def to_include(self):
        foundbykepler = (self.discovery_facility() == "Kepler") | (
            self.discovery_facility() == "K2"
        )
        return foundbykepler


class NonKepler(TransitingExoplanetsSubset):
    def __init__(self, **kw):
        kw = dict( label="Non-Kepler", color="black", zorder=0) | kw

        TransitingExoplanetsSubset.__init__(
            self, **kw
        )

    def to_include(self):
        foundbykepler = (self.discovery_facility() == "Kepler") | (
            self.discovery_facility() == "K2"
        )
        return foundbykepler == False


class TESS(TransitingExoplanetsSubset):
    def __init__(self, **kw):
        kw = dict( label="TESS", color="orange", zorder=0) | kw

        TransitingExoplanetsSubset.__init__(
            self, **kw
        )

    def to_include(self):
        foundbytess = (
            self.discovery_facility() == "Transiting Exoplanet Survey Satellite (TESS)"
        )
        return foundbytess == True


class NonTESS(TransitingExoplanetsSubset):
    def __init__(self, **kw):
        kw = dict(label="NonTESS", color="black", zorder=0) | kw
        TransitingExoplanetsSubset.__init__(
            self,  **kw
        )

    def to_include(self):
        foundbytess = (
            self.discovery_facility() == "Transiting Exoplanet Survey Satellite (TESS)"
        )
        return foundbytess == False

class NonKeplerNonTESS(TransitingExoplanetsSubset):
    def __init__(self, **kw):
        kw = dict(label="No Kepler, No TESS", color="black", zorder=0) | kw
        TransitingExoplanetsSubset.__init__(
            self,  **kw
        )

    def to_include(self):
        foundbytess = (
            self.discovery_facility() == "Transiting Exoplanet Survey Satellite (TESS)"
        )
        foundbykepler = (self.discovery_facility() == "Kepler") | (
            self.discovery_facility() == "K2"
        )
        return (foundbytess == False)*(foundbykepler == False)



space_telescopes = [
    "Transiting Exoplanet Survey Satellite (TESS)",
    "K2",
    "Kepler",
    "CoRoT",
    "Hubble Space Telescope",
]


class Space(TransitingExoplanetsSubset):
    def __init__(self, **kw):
        TransitingExoplanetsSubset.__init__(
            self, label="Space-based", color="orchid", zorder=0, **kw
        )

    def to_include(self):
        foundfromspace = np.zeros(len(self)).astype(bool)
        for x in space_telescopes:
            foundfromspace = foundfromspace | (self.discovery_facility() == x)
        return foundfromspace


class Ground(TransitingExoplanetsSubset):
    def __init__(self, **kw):
        kw = dict(label="Ground-based", color="black", zorder=0) | kw
        TransitingExoplanetsSubset.__init__(
            self,  **kw
        )

    def to_include(self):
        foundfromspace = np.zeros(len(self)).astype(bool)
        for x in space_telescopes:
            foundfromspace = foundfromspace | (self.discovery_facility() == x)
        return foundfromspace == False


sigma = 3


def mass_is_good(pop):
    with np.errstate(invalid="ignore"):
        # the uncertainty must be greater than 0
        exists = pop.get_uncertainty("mass") > 0

        # the uncertainty must be less than a maximum
        fractional = pop.get_uncertainty("mass") / pop.get("mass")
        small = fractional < pop.maximum_uncertainty

        return small & exists


class GoodMass(TransitingExoplanetsSubset):
    def __init__(self, sigma=sigma, **kw):
        self.maximum_uncertainty = 1 / sigma
        kw = dict(label="Good Mass") | kw
        TransitingExoplanetsSubset.__init__(self,  **kw)

    def to_include(self):
        return mass_is_good(self)


class BadMass(TransitingExoplanetsSubset):
    def __init__(self, sigma=sigma, **kw):
        self.maximum_uncertainty = 1 / sigma
        kw = dict(label="Bad Mass", color="lightblue") | kw

        TransitingExoplanetsSubset.__init__(
            self,  **kw
        )

    def to_include(self):
        return mass_is_good(self) == False
