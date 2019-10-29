from .Exoplanets import *

class Subset(TransitingExoplanets):
    def __init__(self, label, **plotkw):

        TransitingExoplanets.__init__(self, **plotkw)

        # set the label
        self.label = label

        # trim to just the data we want
        self.standard = self.standard[self.to_include()]

    def to_include(self):
        raise NotImplementedError('!')


class Kepler(Subset):
    def __init__(self):
        Subset.__init__(self, label="Kepler", color='royalblue', zorder=0)

    def to_include(self):
        foundbykepler = (self.discoverer == 'Kepler') | (self.discoverer == 'K2')
        return foundbykepler

class NonKepler(Subset):
    def __init__(self):
        Subset.__init__(self, label="Non-Kepler", color='black', zorder=0)

    def to_include(self):
        foundbykepler = (self.discoverer == 'Kepler') | (self.discoverer == 'K2')
        return  foundbykepler == False

class TESS(Subset):
    def __init__(self):
        Subset.__init__(self, label="TESS", color='orangered', zorder=0)

    def to_include(self):
        foundbytess = (self.discoverer == 'Transiting Exoplanet Survey Satellite (TESS)')
        return foundbytess == False

space_telescopes = ['Transiting Exoplanet Survey Satellite (TESS)',
                    'K2', 'Kepler',
                    'CoRoT',
                    'Hubble Space Telescope']
class Space(Subset):
    def __init__(self):
        Subset.__init__(self, label="Space-based", color='orangered', zorder=0)

    def to_include(self):
        foundfromspace = np.zeros(self.n).astype(np.bool)
        for x in space_telescopes:
            foundfromspace = foundfromspace | (self.discoverer == x)
        return foundfromspace

class Ground(Subset):
    def __init__(self):
        Subset.__init__(self, label="Ground-based", color='orangered', zorder=0)

    def to_include(self):
        foundfromspace = np.zeros(self.n).astype(np.bool)
        for x in space_telescopes:
            foundfromspace = foundfromspace | (self.discoverer == x)
        return foundfromspace == False


'''


#
# a pair of subsamples, for those with and without good mass measurements
#
threshold = 2.5
def hasMass(pop):
    mass_uncertainty = (np.abs(pop.planet_mass_uncertainty_lower) + np.abs(pop.planet_mass_uncertainty_upper))/2.0
    smallEnough = mass_uncertainty/pop.planet_mass < pop.maximum_uncertainty
    exists = (mass_uncertainty > 0)*np.isfinite(mass_uncertainty)
    return smallEnough*exists

class GoodMass(Subset):
    def __init__(self, threshold=threshold):
        self.maximum_uncertainty = 1.0/threshold
        Subset.__init__(self, label="GoodMass", color='black', zorder=0)
        self.ink=True

    def toRemove(self):
        return hasMass(self) == False

class BadMass(Subset):
    def __init__(self, threshold=threshold):
        self.maximum_uncertainty = 1.0/threshold
        Subset.__init__(self, label="BadMass", color='blue', zorder=-100)
        self.ink=True

    def toRemove(self):
        return hasMass(self) == True

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

class Highlight(TransitingExoplanets):
    def __init__(self, name):
        TransitingExoplanets.__init__(self)
        self.removeRows(np.array([name in n.lower().replace(' ', '') for n in self.name]) == False)
'''
