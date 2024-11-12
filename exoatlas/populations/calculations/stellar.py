from ...imports import *


@property
def stellar_luminosity(self):
    T = self.stellar_teff
    R = self.stellar_radius
    sigma = con.sigma_sb
    return (4 * np.pi * R**2 * sigma * T**4).to(u.Lsun)


@property
def stellar_luminosity_uncertainty(self):
    T = self.stellar_teff
    R = self.stellar_radius
    T_uncertainty = self.get_uncertainty("stellar_teff")
    R_uncertainty = self.get_uncertainty("stellar_radius")

    sigma = con.sigma_sb

    dLdT = 4 * 4 * np.pi * R**2 * sigma * T**3
    dLdR = 2 * 4 * np.pi * R**1 * sigma * T**4
    return ((dLdT**2 * T_uncertainty**2 + dLdR**2 * R_uncertainty**2) ** 0.5).to("Lsun")


@property
def distance_modulus(self):
    """
    The distance modulus to the system, in magnitudes.
    """
    mu = 5 * np.log10(self.distance / (10 * u.pc))
    return mu
