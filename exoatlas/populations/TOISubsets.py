from ..imports import *
from .TOI import *
from .ExoplanetsSubset import *
from astropy.coordinates import SkyCoord
from astropy import units as u

__all__ = ['TOISubset', 'PreviouslyKnownTOI', 'BrandNewTOI']

class TOISubset(TOI):
    def __init__(self, label, **kw):

        TOI.__init__(self, **kw)

        # set the label
        self.label = label

        # trim to just the data we want
        self.standard = self.standard[self.to_include()]

    def to_include(self):
        raise NotImplementedError('!')

class PreviouslyKnownTOI(TOISubset):
    def __init__(self, label="Previously Known TOI", remake=False, color='black', **kw):
        # create a reference population of transiting exoplanets
        self.all_known = NonTESS(remake=remake)

        TOISubset.__init__(self, label=label, remake=remake, color=color, **kw)

    def to_include(self):
        return self.previouslyknown()

    def previouslyknown(self):
        '''
        Identify the planets that are already known from surveys other than TESS.

        Returns
        -------
        known : array of bools
            An array that is True for previously known planets.
        '''


        # create astropy coordinates for both populations
        toi_coords = SkyCoord(ra=self.ra, dec=self.dec)
        exo_coords = SkyCoord(ra=self.all_known.ra, dec=self.all_known.dec)

        # do a spatial cross match on the sky
        #  (idx gives the index into exo_coords, corresponding to every entry in toi_coords)
        idx, d2d, d3d = toi_coords.match_to_catalog_sky(exo_coords)

        # identify which systems are actually close on the sky
        match = d2d < (1*u.arcmin)

        # create new populations that are linked by spatial position
        i_match = match.nonzero()[0]
        matched_toi = self[i_match]
        matched_known = self.all_known[idx[i_match]]

        if False:
            # make some plots for debugging
            for x in [matched_known, matched_toi]:
                plt.scatter(x.ra, x.dec, alpha=0.3)
            plt.xlabel('RA ($^\circ$)'); plt.ylabel('Dec ($^\circ$)')
            plt.show()

            for x in [matched_known, matched_toi]:
                plt.scatter(x.period, x.radius, alpha=0.3, label=x.label)
            plt.yscale('log');
            plt.legend()
            plt.xlabel('Period (days)'); plt.ylabel('Radius (Earth radii)')
            plt.show()

        # find only those planets that have periods close to each other
        dp = (matched_known.period - matched_toi.period)/matched_known.period
        closeperiod = np.abs(dp) < np.inf# FIXME -- kludge!0.1
        matched_toi.standard['dp']  = dp
        matched_toi.standard['close']  = closeperiod


        # populate a mask with Trues where planets are close to another in space and period
        wasknown = np.zeros(self.n).astype(np.bool)
        wasknown[i_match[closeperiod]] = True


        #weird = (matched_toi.period > 13*u.day)*closeperiod
        #print(matched_toi[weird].standard['name', 'period', 'radius', 'dp', 'close'])
        #print(matched_known[weird].standard['name', 'period', 'radius'])

        # return that array
        return wasknown | self.is_knownplanet

class BrandNewTOI(PreviouslyKnownTOI):
    def __init__(self, label="Brand New TOI", **kw):
        PreviouslyKnownTOI.__init__(self, label=label, color='crimson', **kw)

    def to_include(self):
        return self.previouslyknown() == False
