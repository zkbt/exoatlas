"""
A simple cartoon of the Chen & Kipping (2017) forecaster model,
using the values quoted in their Table 2.
"""

__all__ = ["estimate_radius", "estimate_mass", "plot_chen"]

import numpy as np
import matplotlib.pyplot as plt
from astropy.units import Rearth, Mearth, Rjupiter, Mjupiter, Rsun, Msun

# using log polynomials from https://exoplanetarchive.ipac.caltech.edu/docs/pscp_calc.html


def estimate_radius(M):
    """
    Estimate the radii, given an array of masses.

    Parameters
    ----------
    M : astropy.units.quantity.quantity
        The masses of the planets.
    """

    ### RADIUS FROM MASS!
    M_min, M_max, C, S = {}, {}, {}, {}

    # define the terran range
    M_min["terran"] = 0 * Mearth
    M_max["terran"] = 2.04 * Mearth
    S["terran"] = 0.2790
    C["terran"] = 0.00346

    # define the neptunian range
    M_min["neptunian"] = M_max["terran"]
    M_max["neptunian"] = 132 * Mearth
    S["neptunian"] = 0.589
    C["neptunian"] = -0.0925

    # define the jovian range
    M_min["jovian"] = M_max["neptunian"]
    M_max["jovian"] = 26600 * Mearth
    S["jovian"] = -0.044
    C["jovian"] = 1.25

    # define the jovian range
    M_min["stellar"] = M_max["jovian"]
    M_max["stellar"] = 1 * Msun
    S["stellar"] = 0.881
    C["stellar"] = -2.85

    # create an empty array of radii
    original_shape = np.shape(M)
    logM = np.atleast_1d(np.log10(M.to_value(Mearth)))
    logR = np.zeros(np.shape(logM)) * np.nan

    # work our way down in mass, updating radii as we go
    for this in ["stellar", "jovian", "neptunian", "terran"]:
        # which planets are below this maximum?
        x = M < M_max[this]

        # calculate the radius
        logR[x] = C[this] + logM[x] * S[this]

    # return the same shape
    return 10 ** logR.reshape(original_shape) * Rearth


def estimate_mass(R):
    """
    Estimate the mass, given an array of raii.

    Parameters
    ----------
    R : astropy.units.quantity.quantity
        The radii of the planets.
    """

    ### MASS FROM RADIUS!
    R_min, R_max, C, S = {}, {}, {}, {}

    # define the terran range
    R_min["terran"] = 0 * Rearth
    R_max["terran"] = 1.23 * Rearth
    S["terran"] = 0.2790
    C["terran"] = 0.00346

    # define the neptunian range
    R_min["neptunian"] = R_max["terran"]
    R_max["neptunian"] = 11.1 * Rearth
    S["neptunian"] = 0.589
    C["neptunian"] = -0.0925

    # create an empty array of radii
    original_shape = np.shape(R)
    logR = np.atleast_1d(np.log10(R.to_value(Rearth)))
    logM = np.zeros(np.shape(logR)) * np.nan

    # work our way down in mass, updating radii as we go
    for this in ["neptunian", "terran"]:
        # which planets are below this maximum?
        x = R < R_max[this]

        # calculate the radius
        logM[x] = (logR[x] - C[this]) / S[this]

    # return the same shape
    return 10 ** logM.reshape(original_shape) * Mearth


def plot_chen(independent="mass", **kw):
    """
    Plot the MAP mass-radius curve from Chen and Kipping (2017).
    This is mostly a wrapper to make testing easier.
    """

    if independent[0].lower() == "m":
        M = np.logspace(-1, 5, 1000) * Mearth
        R = estimate_radius(M)
        plt.loglog(M, R, **kw)
    elif independent[0].lower() == "r":
        R = np.logspace(-1, 5, 1000) * Rearth
        M = estimate_mass(R)
        plt.loglog(M, R, **kw)


"""def update_masses_with_chen_and_kipping(population, planets):
    # KLUDGE!?!?!
    # MODIFIES POPULATION IN PLACE!
    population.standard.loc
    mass_is_bad = np.isfinite(population.mass) == False
    population.mass[mass_is_bad] = estimate_mass(population.radius[mass_is_bad])
"""
