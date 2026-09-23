from ..imports import *

import pandas as pd
from scipy.interpolate import RegularGridInterpolator


class Tang:
    bibcode = "2025ApJ...989...28T"
    age_strings = {"10Myr": 1e7, "100Myr": 1e8, "1Gyr": 1e9, "10Gyr": 1e10}
    gas_strings = ["0.10%", "0.20%", "0.50%", "1%", "2%", "5%", "10%", "20%"]

    def __repr__(self):
        return f"<interpolator for Tang et al. ({self.bibcode}) planet radius evolution models>"

    def __init__(self, pressure="20mbar", metallicity="50x"):
        """
        Initialize an intepolator for the
        Tang et al. (2025) radius models.

        Parameters
        ----------
        pressure : str
            The pressure level at which radius is measured.
            Options include "1nbar", "20mbar", "RCB".
        metallicity : str
            The metallicity of the gas envelope.
            Options include "1x", "50x".
        """
        self.pressure = pressure
        self.metallicity = metallicity
        grid = {}

        for age_string, age in self.age_strings.items():
            filename = f"data/tang/{self.pressure}_{age_string}_{self.metallicity}.csv"
            path = files(__name__) / filename
            grid[age] = pd.read_csv(path, comment="#")

        N_ages = len(grid)
        one_age = grid[1e7]
        self.ages = np.unique(list(grid.keys()))
        self.fluxes = np.unique(one_age["flux (F_e)"])
        self.core_masses = np.unique(one_age["core mass (M_e)"])
        self.gas_fractions = np.array(
            [float(f.strip("%")) / 100.0 for f in self.gas_strings]
        )

        self.radii = np.empty(
            (
                len(self.ages),
                len(self.fluxes),
                len(self.core_masses),
                len(self.gas_fractions),
            )
        )
        for i_age, age in enumerate(self.ages):
            one_age = grid[age]
            for i_flux, flux in enumerate(self.fluxes):
                is_flux = one_age["flux (F_e)"] == flux
                for i_core, core_mass in enumerate(self.core_masses):
                    is_core = one_age["core mass (M_e)"] == core_mass
                    for i_gas, gas_fraction in enumerate(self.gas_strings):
                        this = one_age[gas_fraction][np.nonzero(is_flux * is_core)[0]]
                        is_this = np.nonzero(is_flux * is_core)[0]
                        assert len(this) == 1
                        this = one_age[gas_fraction].iat[is_this[0]]
                        self.radii[i_age, i_flux, i_core, i_gas] = this

        self.log_interpolator = RegularGridInterpolator(
            [
                np.log10(self.ages),
                np.log10(self.fluxes),
                np.log10(self.core_masses),
                np.log10(self.gas_fractions),
            ],
            np.log10(self.radii),
            bounds_error=False,
            fill_value=np.nan,
        )

    def __call__(self, age=1*u.Gyr,
                 flux=20,
                 core=10*u.M_earth,
                 gas=0.1):
        """
        Interpolate the model radius.

        Inputs can be either scalars or 1D arrays,
        but all arrays must have the same length.

        Parameters
        ----------
        age : float, np.array, u.Quantity
            The age, in year.
        flux : float, np.array
            The bolometric flux received by planet, relative to Earth.
        core : float, np.array, u.Quantity
            The core mass, in Earth masses.
        gas : float, np.array
            The H/He fraction, as fraction of the planet (or core?!?) mass.

        Returns
        -------
        radius : float, array
            The radius of the planet, given input parameters.
        """

        # make sure all inputs are 1D arrays
        a = np.atleast_1d(age.to_value('year'))
        f = np.atleast_1d(flux)
        c = np.atleast_1d(core.to_value('M_earth'))
        g = np.atleast_1d(gas)

        # construct array of inputs onto which we will interpolate
        N = np.max([len(a), len(f), len(c), len(g)])
        inputs = np.empty((N, 4))
        inputs[:, 0] = a
        inputs[:, 1] = f
        inputs[:, 2] = c
        inputs[:, 3] = g

        # do the interpolation (in log space!)
        radius = 10**self.log_interpolator(np.log10(inputs))*u.R_earth

        # return the interpolated result
        return radius

    def plot(self):
        ages = [1e7, 1e8, 1e9, 1e10]*u.year
        fi, ax = plt.subplots(1, len(ages), sharey=True, figsize=(8, 3))
        mass = np.logspace(-1, 2)*u.M_earth
        for i, age in enumerate(ages):
            plt.sca(ax[i])
            for g in [0.001, 0.01, 0.1]:
                plt.plot(mass, self(age=age, core=mass, gas=g), label=f"{g:.1%} H/He")
            plt.xlabel("Mass (Earth masses)")
            plt.title(f"$10^{{{np.log10(age.to_value('year')):.0f}}}$ years")
        plt.legend(frameon=False)
        plt.sca(ax[0])
        plt.ylabel("Radius (Earth radii)")
