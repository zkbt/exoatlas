# exoplanet population of all "confirmed" exoplanets from exoplanet archive
from .imports import *
from .Population import Population
from .curation.KOI import correct


url ='http://exoplanetarchive.ipac.caltech.edu/cgi-bin/nstedAPI/nph-nstedAPI?table=cumulative&format=bar-delimited&select=*'
initial_filename = directories['data'] + 'exoplanetArchiveCumulativeCandidates.psv'

def downloadLatest():
    print('downloading the latest list of confirmed exoplanets from the Exoplanet Archive')
    request.urlretrieve(url, initial_filename)

class KOI(PredefinedPopulation):
    def __init__(self, label='KOI', **kwargs):
        '''Initialize a population of KOI's, from Exoplanet archive.'''

        # set up the population
        Population.__init__(self, label=label, **kwargs)
        correct(self)
        self.save_standard()
        # defing some plotting parameters
        self.color = 'gray'
        self.zorder = -1
        self.ink=True

    def loadFromScratch(self):

        # load from a NASA Exoplanet Archive csv file
        try:
            self.table = ascii.read(initial_filename)
        except IOError:
            downloadLatest()
            self.table = ascii.read(initial_filename)

        self.speak('loaded Exoplanet Archive KOIs from {0}'.format(initial_filename))

        # report original size
        self.speak('original table contains {0} elements'.format(len(self.table)))

    def trim_raw(self):
        ok = (self.table['koi_jmag'] > 1.0)* \
                (self.table['koi_srad'] > 0)* \
                (self.table['koi_disposition'] != 'FALSE POSITIVE')* \
                (self.table['koi_prad'] < 20)* \
                (self.table['koi_prad_err1'] != 0.0)*\
                (self.table['koi_impact'] < 1.0)#*(self.table['pl_rade'] < 3.5)
        self.trimmed = self.table[ok]
        self.speak('trimmed from {0} to {1}'.format(len(self.table), len(self.trimmed)))
        print("  having removed")
        print(self.table[ok == False])



    def create_standard(self):
        t = self.trimmed
        n = len(t)

        s = Table()
        s['name'] = t['kepler_name']#[t['kepler_name'][i] + t['pl_letter'][i] for i in range(len(t))]
        for i in range(n):
            if t['kepler_name'][i] != '0':
                s['name'][i] = t['kepler_name'][i].replace(' ','')
            else:
                s['name'][i] = t['kepoi_name'][i]
        s['period'] = t['koi_period']

        s['transit_duration'] = t['koi_duration']/24.0

        s['stellar_teff'] = t['koi_sstellar_teff']
        s['stellar_radius'] = t['koi_srad']
        s['J'] = t['koi_jmag']

        # planet radius
        s['planet_radius'] = t['koi_prad']
        s['planet_radius_uncertainty_upper'] = t['koi_prad_err1']
        s['planet_radius_uncertainty_lower'] = t['koi_prad_err2']


        #KLUDGE?
        s['a_over_r'] = t['koi_dor']

        #KLUDGE?
        #s['rv_semiamplitude'] =  t['pl_rvamp'] #t.MaskedColumn(t['K'], mask=t['K']==0.0)

        s['planet_mass'] = np.zeros(n) + np.nan
        s['planet_mass_uncertainty_upper'] = np.zeros(n) + np.nan
        s['planet_mass_uncertainty_lower'] = np.zeros(n) + np.nan


        s['radius_ratio'] = t['koi_ror']

        badpos = (t['ra'] ==0.0)*(t['dec'] == 0.0)
        s['ra'] = t.MaskedColumn(t['ra'], mask=badpos)
        s['dec'] = t.MaskedColumn(t['dec'], mask=badpos)

        s['b'] = t['koi_impact']


        s['stellar_distance'] = np.zeros(n) + np.nan
        s['stellar_distance_uncertainty_upper'] = np.zeros(n) + np.nan
        s['stellar_distance_uncertainty_lower'] = np.zeros(n) + np.nan

        s['disposition'] = t['koi_disposition']

        # a little kludge
        #s['stellar_teff'][s['name'] == 'GJ 436b'] = 3400.0
        #s['stellar_teff'][s['name'] == 'Qatar-1b'] = 4860.0
        #s['stellar_radius'][s['name'] == 'WASP-100b'] = 1.5#???
        s.sort('name')
        self.standard = s

class Subset(KOI):
    def __init__(self, label, color='black', zorder=0):

        # set the label
        self.label=label
        self.color=color
        self.zorder=zorder
        try:
            # first try to load this population
            Talker.__init__(self)
            self.load_standard()
        except IOError:
            # if that fails, recreate it from the confirmed population
            KOI.__init__(self)
            self.label=label
            self.selectSubsample()
        self.ink=True

    def selectSubsample(self):
        tr = self.toRemove()
        self.speak('removing {0} rows'.format(np.sum(tr)))
        self.removeRows(tr)
        self.speak('leaving {0} rows'.format(self.n))
        self.save_standard()

class UnconfirmedKepler(Subset):
    def __init__(self):
        Subset.__init__(self, label="Kepler (candidates)", color='gray', zorder=-1e6)
        self.ink=True

    def toRemove(self):
        isconfirmed = self.standard['disposition'] == 'CONFIRMED'
        isjunk = self.distance == 10.0
        return isconfirmed | isjunk
