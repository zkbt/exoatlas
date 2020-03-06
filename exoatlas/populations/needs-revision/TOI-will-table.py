# exoplanets from Will's summary table

from ..imports import *
from .Population import PredefinedPopulation

import pandas as pd
from astroquery.mast import Catalogs


from .downloaders import toi_exofop


#url = 'https://exofop.ipac.caltech.edu/tess/download_ctoi.php?sort=ctoi&output=pipe'
#initial_filename = directories['data'] + 'TOI-from-exofop.psv'
#
#def downloadLatest():
#    print('downloading the latest list of TOI candidates from ExoFOP')
#    request.urlretrieve(url, initial_filename)

class TOI(PredefinedPopulation):

    def __init__(self, label='TOI', remake=False, **kw):
        '''
        Initialize a population of TESS Objects of Interest
        from a table downloaded from the NASA Exoplanet Archive.
        '''
        # set up the population
        PredefinedPopulation.__init__(self, label=label, remake=remake, **kw)


    def load_raw(self, remake=False):
        '''
        Load the raw table of TOI data from the NASA Exoplanet Archive.
        '''

        # load (or download) the table of composite exoplanet properties
        raw = toi_exofop.get(remake=remake)

        # for debugging, hang on to the raw table as a hidden attribute
        self._raw = raw

        return raw

    def create_standard(self, trimmed):
        '''
        Create a standardized table, pulling at least the necessary columns
        from the raw table and potentially including others too.
        '''


        t = trimmed
        n = len(t)
        s = Table()

        # PICK UP FROM HERE!!!!!!!!!!!!!

        # pull out the name as the CTOI
        s['name'] = ['TOI{:.2f}'.format(c) for c in t['CTOI']]
        s['hostname'] = ['TOI{:.0f}'.format(c) for c in t['CTOI']]
        
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
        s['Jmag'] = withtic['TIC Jmag']# KLUDGE!!!!!t['koi_jmag']
        s['Vmag'] = withtic['TIC Vmag']
        s['G'] = withtic['TIC GAIAmag']
        s['T'] = withtic['TIC Tmag']

        # planet radius
        s['radius'] = withtic['Radius (R_Earth)']
        s['radius_uncertainty_upper'] = withtic['Radius (R_Earth) Error']
        s['radius_uncertainty_lower'] = -withtic['Radius (R_Earth) Error']


        #KLUDGE?
        s['transit_ar'] = withtic['a/Rad_s']

        #KLUDGE?
        #s['rv_semiamplitude'] =  t['pl_rvamp'] #t.MaskedColumn(t['K'], mask=t['K']==0.0)

        s['mass'] = withtic['Mass (M_Earth)']
        s['mass_uncertainty_upper'] = withtic['Mass (M_Earth) Error']
        s['mass_uncertainty_lower'] = -withtic['Mass (M_Earth) Error']



        s['radius_ratio'] = withtic['Radius (R_Earth)']*u.Rearth/withtic['Stellar Radius (R_Sun)']/u.Rsun

        s['ra'] = withtic['RA']
        s['dec'] = withtic['Dec']

        s['transit_b'] = withtic['Impact Param']

        # is this usually Gaia??
        s['distance'] = withtic['Stellar Distance (pc)']
        s['distance_uncertainty_upper'] = np.zeros(n) + np.nan
        s['distance_uncertainty_lower'] = np.zeros(n) + np.nan

        s['disposition'] = withtic['User Disposition']

        # a little kludge
        #s['stellar_teff'][s['name'] == 'GJ 436b'] = 3400.0
        #s['stellar_teff'][s['name'] == 'Qatar-1b'] = 4860.0
        #s['stellar_radius'][s['name'] == 'WASP-100b'] = 1.5#???
        s.sort('name')
        self.standard = s

class ExoplanetSubsets(TOI):
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


"""class UnconfirmedKepler(ExoplanetSubsets):
    def __init__(self):
        ExoplanetSubsets.__init__(self, label="Kepler (candidates)", color='gray', zorder=-1e6)
        self.ink=True

    def toRemove(self):
        isconfirmed = self.standard['disposition'] == 'CONFIRMED'
        isjunk = self.distance == 10.0
        return isconfirmed | isjunk"""
