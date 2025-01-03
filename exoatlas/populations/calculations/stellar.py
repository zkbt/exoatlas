from ...imports import *


def stellar_luminosity_from_table(self, distribution=False, **kw):
    """
    Stellar Luminosity (L*, Lsun)

    Calculate the stellar luminosity from the quoted table value,
    which for the NASA Exoplanet Archive is stored as log10(L/Lsun).

    Parameters
    ----------
    distribution : bool
        If False, return a simple array of values.
        If True, return an astropy.uncertainty.Distribution,
        which can be used for error propagation.
    """
    logL = self.get("stellar_logluminosity", distribution=distribution)
    L = 10**logL * u.Lsun
    return L


def stellar_luminosity_from_radius_and_teff(self, distribution=False, **kw):
    """
    Stellar Luminosity (L*, Lsun)

    Calculate "L*" from the radius and effective temperature.
    This might be used if no stellar luminosity is provided
    in the standardized table.

    Parameters
    ----------
    distribution : bool
        If False, return a simple array of values.
        If True, return an astropy.uncertainty.Distribution,
        which can be used for error propagation.
    """
    T = self.get("stellar_teff", distribution=distribution)
    R = self.get("stellar_radius", distribution=distribution)
    sigma = con.sigma_sb
    L = (4 * np.pi * R**2 * sigma * T**4).to(u.Lsun)
    return L


def stellar_luminosity(self, distribution=False, **kw):
    """
    Stellar Luminosity (L*, Lsun)

    Retrieve "L*" first from the standardized table,
    then from radius + effective temperature.

    Parameters
    ----------
    distribution : bool
        If False, return a simple array of values.
        If True, return an astropy.uncertainty.Distribution,
        which can be used for error propagation.
    """
    L = self._choose_calculation(
        methods=[
            "stellar_luminosity_from_table",
            "stellar_luminosity_from_radius_and_teff",
        ],
        distribution=distribution,
        **kw,
    )
    return L


def distance_modulus(self, distribution=False, **kw):
    """
    Distance Modulus ($\mu$, magnitudes)

    Parameters
    ----------
    distribution : bool
        If False, return a simple array of values.
        If True, return an astropy.uncertainty.Distribution,
        which can be used for error propagation.
    """
    mu = 5 * np.log10(self.distance(distribution=distribution) / (10 * u.pc)) * u.mag
    return mu
