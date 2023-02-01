"""
A simple cartoon of the Chen & Kipping (2017) forecaster model,
using the values quoted in their Table 2.
"""

__all__ = ["estimate_radius", "plot_chen"]

import numpy as np
import matplotlib.pyplot as plt
from astropy.units import Rearth, Mearth, Rjupiter, Mjupiter, Rsun, Msun

T_min, T_max, C, S = {}, {}, {}, {}

# define the terran range
T_min["terran"] = 1.0 * Mearth
T_max["terran"] = 2.04 * Mearth
S["terran"] = 0.2790
C["terran"] = 1.008 * Rearth

# define the neptunian range
T_min["neptunian"] = T_max["terran"]
T_max["neptunian"] = 0.414 * Mjupiter
S["neptunian"] = 0.589
C["neptunian"] = (
    C["terran"] * (T_min["neptunian"] / T_min["terran"]).decompose() ** S["terran"]
)

# define the jovian range
T_min["jovian"] = T_max["neptunian"]
T_max["jovian"] = 0.0800 * Msun
S["jovian"] = -0.044
C["jovian"] = (
    C["neptunian"]
    * (T_min["jovian"] / T_min["neptunian"]).decompose() ** S["neptunian"]
)

# define the jovian range
T_min["stellar"] = T_max["jovian"]
T_max["stellar"] = 1 * Msun
S["stellar"] = 0.881
C["stellar"] = (
    C["jovian"] * (T_min["stellar"] / T_min["jovian"]).decompose() ** S["jovian"]
)


def estimate_radius(M):
    """
    Estimate the radii, given an array of masses.

    Parameters
    ----------
    M : astropy.units.quantity.quantity
        The masses of the planets.
    """

    # create an empty array of radii
    original_shape = np.shape(M)
    M = np.atleast_1d(M)
    R = np.zeros(np.shape(M)) * Rearth

    # work our way down in mass, updating radii as we go
    for this in ["stellar", "jovian", "neptunian", "terran"]:
        # which planets are below this maximum?
        x = M < T_max[this]

        # calculate the radius
        R[x] = C[this] * (M[x] / T_min[this]).decompose() ** S[this]

    # return the same shape
    return R.reshape(original_shape)


def plot_chen(**kw):
    """
    Plot the MAP mass-radius curve from Chen and Kipping (2017).
    This is mostly a wrapper to make testing easier.
    """

    M = np.logspace(-1, 5, 1000) * Mearth
    R = estimate_radius(M)
    plt.loglog(M, R, **kw)
