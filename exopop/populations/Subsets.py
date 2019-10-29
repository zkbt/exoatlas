from .Exoplanets import *

__all__ = ['Subset',
           'Kepler', 'NonKepler',
           'TESS', 'NonTESS',
           'Space', 'Ground',
           'GoodMass', 'BadMass']

class Subset(TransitingExoplanets):
    def __init__(self, label, **kw):

        TransitingExoplanets.__init__(self, **kw)

        # set the label
        self.label = label

        # trim to just the data we want
        self.standard = self.standard[self.to_include()]

    def to_include(self):
        raise NotImplementedError('!')

class Kepler(Subset):
    def __init__(self, **kw):
        Subset.__init__(self, label="Kepler", color='royalblue', zorder=0, **kw)

    def to_include(self):
        foundbykepler = (self.discoverer == 'Kepler') | (self.discoverer == 'K2')
        return foundbykepler

class NonKepler(Subset):
    def __init__(self, **kw):
        Subset.__init__(self, label="Non-Kepler", color='black', zorder=0, **kw)

    def to_include(self):
        foundbykepler = (self.discoverer == 'Kepler') | (self.discoverer == 'K2')
        return  foundbykepler == False

class TESS(Subset):
    def __init__(self, **kw):
        Subset.__init__(self, label="TESS", color='orangered', zorder=0, **kw)

    def to_include(self):
        foundbytess = (self.discoverer == 'Transiting Exoplanet Survey Satellite (TESS)')
        return foundbytess == True


class NonTESS(Subset):
    def __init__(self, **kw):
        Subset.__init__(self, label="TESS", color='orchid', zorder=0, **kw)

    def to_include(self):
        foundbytess = (self.discoverer == 'Transiting Exoplanet Survey Satellite (TESS)')
        return foundbytess == False


space_telescopes = ['Transiting Exoplanet Survey Satellite (TESS)',
                    'K2', 'Kepler',
                    'CoRoT',
                    'Hubble Space Telescope']
class Space(Subset):
    def __init__(self, **kw):
        Subset.__init__(self, label="Space-based", color='orangered', zorder=0, **kw)

    def to_include(self):
        foundfromspace = np.zeros(self.n).astype(np.bool)
        for x in space_telescopes:
            foundfromspace = foundfromspace | (self.discoverer == x)
        return foundfromspace

class Ground(Subset):
    def __init__(self, **kw):
        Subset.__init__(self, label="Ground-based", color='orangered', zorder=0, **kw)

    def to_include(self):
        foundfromspace = np.zeros(self.n).astype(np.bool)
        for x in space_telescopes:
            foundfromspace = foundfromspace | (self.discoverer == x)
        return foundfromspace == False

sigma = 2.5
def mass_is_good(pop):
    # the uncertainty must be greater than 0
    exists = pop.uncertainty('planet_mass') > 0

    # the uncertainty must be less than a maximum
    fractional = (pop.uncertainty('planet_mass')/pop.planet_mass)
    small = fractional < pop.maximum_uncertainty

    return small & exists

class GoodMass(Subset):
    def __init__(self, sigma=sigma, **kw):
        self.maximum_uncertainty = 1/sigma
        Subset.__init__(self, label="Good Mass", **kw)
    def to_include(self):
        return mass_is_good(self)

class BadMass(Subset):
    def __init__(self, sigma=sigma, **kw):
        self.maximum_uncertainty = 1/sigma
        Subset.__init__(self, label="Bad Mass", **kw)
    def to_include(self):
        return mass_is_good(self) == False

'''
class lateM(Subset):
    def __init__(self, threshold=threshold):
        Subset.__init__(self, label="T$_{eff}$<3400K", color='darkred', zorder=-100)

    def toRemove(self):
        return self.stellar_teff > 3400

class earlyM(Subset):
    def __init__(self, threshold=threshold):
        Subset.__init__(self, label="3400K<T$_{eff}$<3800K", color='darkred', zorder=-100)

    def toRemove(self):
        return (self.stellar_teff > 3800) | (self.stellar_teff < 3400)

class M(Subset):
    def __init__(self, threshold=threshold):
        Subset.__init__(self, label="M", color='darkred', zorder=-100)

    def toRemove(self):
        return (self.stellar_teff > 3800)

class K(Subset):
    def __init__(self, threshold=threshold):
        Subset.__init__(self, label="K", color='darkred', zorder=-100)

    def toRemove(self):
        return (self.stellar_teff > 5300) | (self.stellar_teff < 3800)

class G(Subset):
    def __init__(self, threshold=threshold):
        Subset.__init__(self, label="G", color='darkred', zorder=-100)

    def toRemove(self):
        return (self.stellar_teff > 6000) | (self.stellar_teff < 5300)

class F(Subset):
    def __init__(self, threshold=threshold):
        Subset.__init__(self, label="F", color='darkred', zorder=-100)

    def toRemove(self):
        return (self.stellar_teff > 7200) | (self.stellar_teff < 6000)
'''
