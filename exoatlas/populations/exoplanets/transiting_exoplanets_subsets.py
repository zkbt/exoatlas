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
        self.standard = self.standard[self.to_include()]

    def to_include(self):
        raise NotImplementedError(
            "Please define `.to_include()` for this ExoplanetSubset!"
        )


class Kepler(TransitingExoplanetsSubset):
    def __init__(self, **kw):
        TransitingExoplanetsSubset.__init__(
            self, label="Kepler", color="royalblue", zorder=0, **kw
        )

    def to_include(self):
        foundbykepler = (self.discovery_facility == "Kepler") | (
            self.discovery_facility == "K2"
        )
        return foundbykepler


class NonKepler(TransitingExoplanetsSubset):
    def __init__(self, **kw):
        TransitingExoplanetsSubset.__init__(
            self, label="Non-Kepler", color="black", zorder=0, **kw
        )

    def to_include(self):
        foundbykepler = (self.discovery_facility == "Kepler") | (
            self.discovery_facility == "K2"
        )
        return foundbykepler == False


class TESS(TransitingExoplanetsSubset):
    def __init__(self, **kw):
        TransitingExoplanetsSubset.__init__(
            self, label="TESS", color="orangered", zorder=0, **kw
        )

    def to_include(self):
        foundbytess = (
            self.discovery_facility == "Transiting Exoplanet Survey Satellite (TESS)"
        )
        return foundbytess == True


class NonTESS(TransitingExoplanetsSubset):
    def __init__(self, **kw):
        TransitingExoplanetsSubset.__init__(
            self, label="TESS", color="black", zorder=0, **kw
        )

    def to_include(self):
        foundbytess = (
            self.discovery_facility == "Transiting Exoplanet Survey Satellite (TESS)"
        )
        return foundbytess == False


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
        foundfromspace = np.zeros(self.n).astype(bool)
        for x in space_telescopes:
            foundfromspace = foundfromspace | (self.discovery_facility == x)
        return foundfromspace


class Ground(TransitingExoplanetsSubset):
    def __init__(self, **kw):
        TransitingExoplanetsSubset.__init__(
            self, label="Ground-based", color="black", zorder=0, **kw
        )

    def to_include(self):
        foundfromspace = np.zeros(self.n).astype(bool)
        for x in space_telescopes:
            foundfromspace = foundfromspace | (self.discovery_facility == x)
        return foundfromspace == False


sigma = 2.5


def mass_is_good(pop):
    with np.errstate(invalid="ignore"):
        # the uncertainty must be greater than 0
        exists = pop.get_uncertainty("mass") > 0

        # the uncertainty must be less than a maximum
        fractional = pop.get_uncertainty("mass") / pop.mass
        small = fractional < pop.maximum_uncertainty

        return small & exists


class GoodMass(TransitingExoplanetsSubset):
    def __init__(self, sigma=sigma, **kw):
        self.maximum_uncertainty = 1 / sigma
        TransitingExoplanetsSubset.__init__(self, label="Good Mass", **kw)

    def to_include(self):
        return mass_is_good(self)


class BadMass(TransitingExoplanetsSubset):
    def __init__(self, sigma=sigma, **kw):
        self.maximum_uncertainty = 1 / sigma
        TransitingExoplanetsSubset.__init__(self, label="Bad Mass", color="gray", **kw)

    def to_include(self):
        return mass_is_good(self) == False


"""
class lateM(TransitingExoplanetsSubset):
    def __init__(self, threshold=threshold):
        TransitingExoplanetsSubset.__init__(self, label="T$_{eff}$<3400K", color='darkred', zorder=-100)

    def toRemove(self):
        return self.stellar_teff > 3400

class earlyM(TransitingExoplanetsSubset):
    def __init__(self, threshold=threshold):
        TransitingExoplanetsSubset.__init__(self, label="3400K<T$_{eff}$<3800K", color='darkred', zorder=-100)

    def toRemove(self):
        return (self.stellar_teff > 3800) | (self.stellar_teff < 3400)

class M(TransitingExoplanetsSubset):
    def __init__(self, threshold=threshold):
        TransitingExoplanetsSubset.__init__(self, label="M", color='darkred', zorder=-100)

    def toRemove(self):
        return (self.stellar_teff > 3800)

class K(TransitingExoplanetsSubset):
    def __init__(self, threshold=threshold):
        TransitingExoplanetsSubset.__init__(self, label="K", color='darkred', zorder=-100)

    def toRemove(self):
        return (self.stellar_teff > 5300) | (self.stellar_teff < 3800)

class G(TransitingExoplanetsSubset):
    def __init__(self, threshold=threshold):
        TransitingExoplanetsSubset.__init__(self, label="G", color='darkred', zorder=-100)

    def toRemove(self):
        return (self.stellar_teff > 6000) | (self.stellar_teff < 5300)

class F(TransitingExoplanetsSubset):
    def __init__(self, threshold=threshold):
        TransitingExoplanetsSubset.__init__(self, label="F", color='darkred', zorder=-100)

    def toRemove(self):
        return (self.stellar_teff > 7200) | (self.stellar_teff < 6000)
"""
