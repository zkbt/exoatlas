# exoplanet population of all "confirmed" exoplanets from exoplanet archive
from imports import *
from Population import Population

from curation.Confirmed import correct


url = 'http://exoplanetarchive.ipac.caltech.edu/cgi-bin/nstedAPI/nph-nstedAPI?table=exoplanets&format=bar-delimited&select=*'

initial_filename = directories['data'] + 'exoplanetArchiveConfirmedPlanets.psv'
def downloadLatest():
    print 'downloading the latest list of confirmed exoplanets from the Exoplanet Archive'
    urllib.urlretrieve(url, initial_filename)

class Confirmed(Population):
    def __init__(self, label='Confirmed', **kwargs):
        '''Initialize a population of Known planets, from the NASA Exoplanet archive.'''

        # set up the population
        Population.__init__(self, label=label, **kwargs)
        correct(self)
        self.saveStandard()
        # defing some plotting parameters
        self.color = 'black'
        self.zorder = -1

    def loadFromScratch(self):

        # load from a NASA Exoplanet Archive csv file
        try:
            self.table = astropy.io.ascii.read(initial_filename)
        except IOError:
            downloadLatest()
            self.table = astropy.io.ascii.read(initial_filename)

        self.speak('loaded Exoplanet Archive planets from {0}'.format(initial_filename))

        # report original size
        self.speak('original table contains {0} elements'.format(len(self.table)))

        # trim non-transiting
        good = self.table['pl_tranflag'] == 1
        self.table = self.table[good]
        self.speak('trimmed to {0} transiting planets'.format(np.sum(good)))

    def trimRaw(self):
        ok = (self.table['st_j'] > 1.0)*(self.table['st_rad'] > 0)#*(self.table['pl_rade'] < 3.5)
        self.trimmed = self.table[ok]
        self.speak('trimmed from {0} to {1}'.format(len(self.table), len(self.trimmed)))
        print "  having removed"
        print self.table[ok == False]

    def createStandard(self):
        t = self.trimmed
        s = astropy.table.Table()
        s['name'] = [t['pl_hostname'][i] + t['pl_letter'][i] for i in range(len(t))]

        s['period'] = t['pl_orbper']

        s['teff'] = t['st_teff']
        s['stellar_radius'] = t['st_rad']
        s['J'] = t['st_j']

        # planet radius
        s['planet_radius'] = t['pl_rade']
        s['planet_radius_upper'] = t['pl_radeerr1']
        s['planet_radius_lower'] = t['pl_radeerr2']


        #KLUDGE?
        s['a_over_r'] = t['pl_ratdor']


        #KLUDGE?
        s['rv_semiamplitude'] =  t['pl_rvamp'] #t.MaskedColumn(t['K'], mask=t['K']==0.0)

        s['planet_mass'] = t['pl_masse']
        s['planet_mass_upper'] = t['pl_masseerr1']
        s['planet_mass_lower'] = t['pl_masseerr2']


        s['radius_ratio'] = t['pl_ratror']

        badpos = (t['ra'] ==0.0)*(t['dec'] == 0.0)
        s['ra'] = t.MaskedColumn(t['ra'], mask=badpos)
        s['dec'] = t.MaskedColumn(t['dec'], mask=badpos)

        s['b'] = t['pl_imppar']


        s['stellar_distance'] = t['st_dist']
        s['stellar_distance_upper'] = t['st_disterr1']
        s['stellar_distance_lower'] = t['st_disterr2']

        # a little kludge
        #s['teff'][s['name'] == 'GJ 436b'] = 3400.0
        #s['teff'][s['name'] == 'Qatar-1b'] = 4860.0
        #s['stellar_radius'][s['name'] == 'WASP-100b'] = 1.5#???
        s.sort('name')
        self.standard = s

class Subset(Confirmed):
    def __init__(self, label, color='black', zorder=0):

        # set the label
        self.label=label
        self.color=color
        self.zorder=zorder
        try:
            # first try to load this population
            Talker.__init__(self)
            self.loadStandard()
        except IOError:
            # if that fails, recreate it from the confirmed population
            Confirmed.__init__(self)
            self.label=label
            self.selectSubsample()

    def selectSubsample(self):
        tr = self.toRemove()
        self.speak('removing {0} rows'.format(np.sum(tr)))
        self.removeRows(tr)
        self.speak('leaving {0} rows'.format(self.n))
        self.saveStandard()

#
# a pair of subsamples, for those discovered by Kepler or not
#

stringsIndicatingKepler = ['Kepler', 'K2', 'KIC', 'KOI', 'PH', '116454']
def discoveredByKepler(pop):
    d = np.zeros(len(pop.standard)).astype(np.bool)
    for s in stringsIndicatingKepler:
        match = np.array([s in name for name in pop.standard['name']])
        d = d + match
    return d

class Kepler(Subset):
    def __init__(self):
        Subset.__init__(self, label="Kepler", color='black', zorder=0)

    def toRemove(self):
        return discoveredByKepler(self) == False

class NonKepler(Subset):
    def __init__(self):
        Subset.__init__(self, label="NonKepler", color='blue', zorder=0)

    def toRemove(self):
        return discoveredByKepler(self) == True

#
# a pair of subsamples, for those with and without good mass measurements
#
threshold = 2.5
def hasMass(pop):
    mass_uncertainty = (np.abs(pop.planet_mass_lower) + np.abs(pop.planet_mass_upper))/2.0
    smallEnough = mass_uncertainty/pop.planet_mass < pop.maximum_uncertainty
    exists = (mass_uncertainty > 0)*np.isfinite(mass_uncertainty)
    return smallEnough*exists

class GoodMass(Subset):
    def __init__(self, threshold=threshold):
        self.maximum_uncertainty = 1.0/threshold
        Subset.__init__(self, label="GoodMass", color='black', zorder=0)

    def toRemove(self):
        return hasMass(self) == False

class BadMass(Subset):
    def __init__(self, threshold=threshold):
        self.maximum_uncertainty = 1.0/threshold
        Subset.__init__(self, label="BadMass", color='blue', zorder=-100)

    def toRemove(self):
        return hasMass(self) == True
