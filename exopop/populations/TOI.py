

# exoplanet population of all "confirmed" exoplanets from exoplanet archive
from .imports import *
from .Population import Population
import pandas as pd
from astropy.table import Table, join
from astroquery.mast import Catalogs



#from .curation.KOI import correct


url = 'https://exofop.ipac.caltech.edu/tess/download_ctoi.php?sort=ctoi&output=pipe'
initial_filename = directories['data'] + 'TOI-from-exofop.psv'

def downloadLatest():
    print('downloading the latest list of TOI candidates from ExoFOP')
    request.urlretrieve(url, initial_filename)

class TOI(Population):
    def __init__(self, label='TOI', **kwargs):
        '''Initialize a population of KOI's, from Exoplanet archive.'''

        # set up the population
        Population.__init__(self, label=label, **kwargs)
        #correct(self)
        self.saveStandard()
        # defing some plotting parameters
        self.color = 'gray'
        self.zorder = -1
        self.ink=True

    def loadFromScratch(self):

        # load from a NASA Exoplanet Archive csv file
        try:
            t = pd.read_csv(initial_filename, sep='|')
            self.table = Table.from_pandas(t)
        except IOError:
            downloadLatest()
            t = pd.read_csv(initial_filename, sep='|')
            self.table = Table.from_pandas(t)
        self.speak('loaded TOIs from {0}'.format(initial_filename))

        # report original size
        self.speak('original table contains {0} elements'.format(len(self.table)))

    def trimRaw(self):
        self.trimmed = self.table


    def createStandard(self):
        t = self.trimmed
        n = len(t)

        s = Table()

        # pull out the name as the CTOI
        s['name'] = ['CTOI{:.2f}'.format(c) for c in t['CTOI']]
        s['period'] = t['Period (days)']
        s['transit_duration'] = t['Duration (hrs)']/24

        s['stellar_teff'] = t['Stellar Eff Temp (K)']
        s['stellar_radius'] = t['Stellar Radius (R_Sun)']

        #t['TIC'].name = 'TIC ID'

        # search for the table from MAST
        print('downloading the TIC. this may take a minute? or restart?')
        tic_table = Catalogs.query_criteria(catalog="Tic", ID=np.unique(t['TIC ID']))

        # preface all the columns with TIC so we can keep them straight
        for k in tic_table.colnames:
            tic_table[k].name = f'TIC {k}'

        # make sure the indices are integers
        tic_table['TIC ID'] = np.array(tic_table['TIC ID']).astype(np.int)

        withtic = join(t, tic_table, 'TIC ID', join_type='left')

        # pull out some magnitudes
        s['J'] = withtic['TIC Jmag']# KLUDGE!!!!!t['koi_jmag']
        s['V'] = withtic['TIC Vmag']
        s['G'] = withtic['TIC GAIAmag']
        s['T'] = withtic['TIC Tmag']

        # planet radius
        s['planet_radius'] = withtic['Radius (R_Earth)']
        s['planet_radius_upper'] = withtic['Radius (R_Earth) Error']
        s['planet_radius_lower'] = -withtic['Radius (R_Earth) Error']


        #KLUDGE?
        s['a_over_r'] = withtic['a/Rad_s']

        #KLUDGE?
        #s['rv_semiamplitude'] =  t['pl_rvamp'] #t.MaskedColumn(t['K'], mask=t['K']==0.0)

        s['planet_mass'] = withtic['Mass (M_Earth)']
        s['planet_mass_upper'] = withtic['Mass (M_Earth) Error']
        s['planet_mass_lower'] = -withtic['Mass (M_Earth) Error']



        s['radius_ratio'] = withtic['Radius (R_Earth)']*u.Rearth/withtic['Stellar Radius (R_Sun)']/u.Rsun

        s['ra'] = withtic['RA']
        s['dec'] = withtic['Dec']

        s['b'] = withtic['Impact Param']

        # is this usually Gaia??
        s['stellar_distance'] = withtic['Stellar Distance (pc)']
        s['stellar_distance_upper'] = np.zeros(n) + np.nan
        s['stellar_distance_lower'] = np.zeros(n) + np.nan

        s['disposition'] = withtic['User Disposition']

        # a little kludge
        #s['stellar_teff'][s['name'] == 'GJ 436b'] = 3400.0
        #s['stellar_teff'][s['name'] == 'Qatar-1b'] = 4860.0
        #s['stellar_radius'][s['name'] == 'WASP-100b'] = 1.5#???
        s.sort('name')
        self.standard = s

class Subset(TOI):
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
        self.saveStandard()


"""class UnconfirmedKepler(Subset):
    def __init__(self):
        Subset.__init__(self, label="Kepler (candidates)", color='gray', zorder=-1e6)
        self.ink=True

    def toRemove(self):
        isconfirmed = self.standard['disposition'] == 'CONFIRMED'
        isjunk = self.distance == 10.0
        return isconfirmed | isjunk"""
