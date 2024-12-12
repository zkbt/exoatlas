from ...imports import *


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
    """
    mu = 5 * np.log10(self.distance() / (10 * u.pc)) * u.mag
    return mu
