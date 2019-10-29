# exoplanet population of all "confirmed" exoplanets from exoplanet archive
from ..imports import *
from .Population import PredefinedPopulation

#from .curation.TransitingExoplanets import correct
from .downloaders import merged_exoplanets

class Exoplanets(PredefinedPopulation):
    def __init__(self, label='All Exoplanets', remake=False, **plotkw):
        '''
        Initialize a population of all known exoplanets,
        from a table downloaded from the NASA Exoplanet Archive.
        '''
        # set up the population
        PredefinedPopulation.__init__(self, label=label, remake=remake, **plotkw)

    def load_raw(self, remake=False):
        '''
        Load the raw table of data from the NASA Exoplanet Archive.
        '''

        # load (or download) the table of composite exoplanet properties
        raw = merged_exoplanets.get(remake=remake)

        # for debugging, hang on to the raw table as a hidden attribute
        self._raw = raw

        return raw

    def trim_raw(self, raw):
        '''
        Trim the raw table down to ones with reasonable values.
        '''

        # trim the raw table
        N = len(raw)
        ok = np.ones(N).astype(np.bool)
        trimmed = raw[ok]

        # for debugging, hang onto the trimmed table as a hidden attribute
        self._trimmed = trimmed

        # report how many points got trimmed away
        ntotal = len(raw)
        ntrimmed = len(trimmed)
        self.speak(f'trimmed down to {ntrimmed}/{ntotal} rows')

        return trimmed

    def create_standard(self, trimmed):
        '''
        Create a standardized table, pulling at least the necessary columns
        from the raw table and potentially including others too.
        '''

        # define the table from which we're deriving everying
        t = trimmed

        s = Table()

        # what's the name of the planet?
        s['name'] = [x.replace(' ', '') for x in t['pl_name']]

        #badpos = (t['ra'] ==0.0)*(t['dec'] == 0.0)
        s['ra'] = t['ra']*u.deg#t.MaskedColumn(t['ra'], mask=badpos)
        s['dec'] = t['dec']*u.deg#t.MaskedColumn(t['dec'], mask=badpos)


        # what's the period of the planet?
        s['period'] = t['pl_orbper']*u.day

        # what are the observed transit properties?
        s['transit_epoch'] = t['pl_tranmid']*u.day
        s['transit_duration'] = t['pl_trandur']*u.day
        s['transit_depth'] = t['pl_trandep']

        # what are the basic stellar properties?
        s['stellar_teff'] = t['st_teff']*u.K
        s['stellar_radius'] = t['st_rad']*u.Rsun
        s['stellar_mass'] = t['st_mass']*u.Msun

        # what are the stellar magnitudes?
        s['J'] = t['st_j']

        # what is the planet radius?
        s['planet_radius'] = t['pl_rade']*u.Rearth

        # pull out the radius uncertainties
        upper = t['pl_radeerr1'].data
        lower = t['pl_radeerr2'].data
        bad = upper.mask | lower.mask
        upper[bad] = np.inf
        lower[bad] = np.inf
        s['planet_radius_uncertainty_upper'] = upper*u.Rearth
        s['planet_radius_uncertainty_lower'] = lower*u.Rearth

        # what are the (often) transit-derived properties?
        s['a_over_r'] = t['pl_ratdor']
        s['b'] = t['pl_imppar']

        #KLUDGE?
        #rsovera = (3*np.pi/u.G/period**2/stellar_density/(1.0+mass_ratio))**(1.0/3.0)


        '''
        bad = (np.isfinite(s['a_over_r']) == False)+( s['a_over_r'].data == 0.0)
        period = s['period']
        stellar_radius = s['stellar_radius']
        stellar_mass = s['stellar_mass']
        otherestimate = (craftroom.units.G*(period*craftroom.units.day)**2*(stellar_mass*craftroom.units.Msun)/4/np.pi**2/(stellar_radius*craftroom.units.Rsun)**3)**(1./3.)

        s['a_over_r'][bad] = otherestimate[bad]
        '''

        #KLUDGE?
        s['rv_semiamplitude'] =  t['pl_rvamp'] #t.MaskedColumn(t['K'], mask=t['K']==0.0)

        s['planet_mass'] = t['pl_masse']*u.Mearth

        # pull out the mass uncertainties
        upper = t['pl_masseerr1'].data
        lower = t['pl_masseerr2'].data
        bad = upper.mask | lower.mask
        upper[bad] = np.inf
        lower[bad] = np.inf
        s['planet_mass_uncertainty_upper'] = upper*u.Mearth
        s['planet_mass_uncertainty_lower'] = lower*u.Mearth


        # how far away is the star?
        s['stellar_distance'] = t['st_dist']*u.pc
        s['stellar_distance_uncertainty_upper'] = t['st_disterr1']*u.pc
        s['stellar_distance_uncertainty_lower'] = t['st_disterr2']*u.pc

        # what facility discovered the planet?
        s['discoverer'] = t['pl_facility']

        # sort these planets by their names
        s.sort('name')

        # set the fill_value for any numeric columns to nans
        for k in s.colnames:
            if not isinstance(s[k][0], str):
                s[k].fill_value = np.nan

        # fill in all the masked elements to make an unmasked array with nans
        standard = s.filled()

        # return that standardized table
        return standard

class TransitingExoplanets(Exoplanets):
    def __init__(self, label='Transiting Exoplanets', remake=False, **plotkw):
        '''
        Initialize a population of all known transiting exoplanets,
        from a table downloaded from the NASA Exoplanet Archive.
        '''
        # set up the population
        Exoplanets.__init__(self, label=label, remake=remake, **plotkw)

    def trim_raw(self, raw):
        '''
        Trim the raw table down to ones with reasonable values.
        '''

        masks = {}

        # does this planet transit?
        masks['transits'] = raw['pl_tranflag'] == 1

        # does this planet have a J magnitude?
        masks['has_J'] = raw['st_j'] > 1.0

        # does this planet have a stellar radius?
        masks['has_radius'] = raw['st_rad'] > 0.0


        ok = np.ones(len(raw)).astype(np.bool)
        for k in masks:
            ok *= masks[k]
            N = sum(ok == True)
            self.speak(f'{N} planets pass the `{k}` filter')

        # trim down the table to just those that are OK
        trimmed = raw[ok]

        # for debugging, hang onto the trimmed table as a hidden attribute
        self._trimmed = trimmed

        ntotal = len(raw)
        ntrimmed = len(trimmed)
        self.speak(f'trimmed down to {ntrimmed}/{ntotal} rows')

        return trimmed
