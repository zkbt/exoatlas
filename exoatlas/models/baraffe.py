from ..imports import *

from scipy.interpolate import LinearNDInterpolator


class Baraffe:
    bibcode = "2015A&A...577A..42B"

    def __repr__(self):
        return f"<interpolator for Baraffe et al. ({self.bibcode}) stellar evolution models>"

    def __init__(self, quantities=["luminosity", "teff", "radius"]):
        """
        Initialize an intepolator for the
        Baraffe et al. (2015) stellar models.

        Parameters
        ----------

        """
        filename = "data/baraffe/BHAC15_tracks.txt"
        path = files(__name__) / filename
        table = ascii.read(
            path,
            comment="!",
            data_start=5,
            names=[
                "mass",
                "logage",
                "teff",
                "logluminosity",
                "logg",
                "radius",
                "lithium",
            ],
        )
        table["logmass"] = np.log10(table["mass"])
        table["logteff"] = np.log10(table["teff"])
        table["logradius"] = np.log10(table["radius"])
        self.units = dict(
            mass=u.Msun, age=u.year, luminosity=u.Lsun, teff=u.K, radius=u.Rsun
        )
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            table["loglithium"] = np.log10(table["lithium"])

        self.log_interpolators = {}
        for k in quantities:
            self.log_interpolators[k] = LinearNDInterpolator(
                points=np.array([table["logmass"], table["logage"]]).T,
                values=table[f"log{k}"],
                fill_value=np.nan,
            )

    def __call__(self, mass=1*u.Msun, age=4.5*u.Gyr):
        """
        Interpolate stellar models to given mass and age.

        Inputs can be either scalars or 1D arrays,
        but all arrays must have the same length.

        Parameters
        ----------
        mass : float, np.array
            The mass, in solar masses.
        age : float, np.array
            The age in yr.

        Returns
        -------
        quantities : tuple
            An array for each of the quantities specified
            during initialization. By default, this will
            be ['luminosity', 'teff', 'radius'], meaning the
            function will return three arrays, one for each,
            corresponding to the input mass(es) and age(s).
        """

        # make sure all inputs are 1D arrays
        m = np.atleast_1d(mass.to_value('Msun'))
        a = np.atleast_1d(age.to_value('year'))

        # construct array of inputs onto which we will interpolate
        N = np.max([len(m), len(a)])
        inputs = np.empty((N, 2))
        inputs[:, 0] = m
        inputs[:, 1] = a

        # do the interpolation (in log space!)
        return [
            10 ** interpolator(np.log10(inputs)) * self.units[k]
            for k, interpolator in self.log_interpolators.items()
        ]

    def imshow(self):
        N = len(self.log_interpolators)
        fi, ax = plt.subplots(N, 1, figsize=(6, N * 4), sharex=True, sharey=True)
        logmass = np.linspace(-1.1, 0.1, 1000)
        logage = np.linspace(6, 10, 1000)
        logmass_2d, logage_2d = np.meshgrid(logmass, logage)
        outputs = [
            z.reshape(logmass_2d.shape)
            for z in self(
                mass=10 ** logmass_2d.flatten()*u.Msun, age=10 ** logage_2d.flatten()*u.year
            )
        ]

        for i, k in enumerate(self.log_interpolators):
            plt.sca(ax[i])
            plt.pcolormesh(logmass, logage, np.log10(outputs[i].value), shading="auto")
            plt.xlabel("log(mass)")
            plt.ylabel("log(age)")
            plt.title(f"{k}")
            plt.colorbar(label=f"log({k}/{outputs[i].unit})")

    def plot(self):
        N = len(self.log_interpolators)
        fi, ax = plt.subplots(
            N,
            2,
            figsize=(6, N * 2),
            constrained_layout=True,
            sharex="col",
            sharey="row",
        )
        for mass in [0.1, 0.4, 0.7, 1.0]*u.Msun:
            ages = np.logspace(6, 10, 1000)*u.year
            outputs = self(mass=mass, age=ages)
            for i, k in enumerate(self.log_interpolators):
                plt.sca(ax[i, 0])
                plt.loglog(ages, outputs[i], label=mass)
                plt.ylabel(k)

            plt.xlabel("Age (yr)")
            plt.legend()
        for age in np.logspace(6, 10, 5)*u.year:
            masses = np.logspace(-1, 0, 1000)*u.Msun
            outputs = self(mass=masses, age=age)
            for i, k in enumerate(self.log_interpolators):
                plt.sca(ax[i, 1])
                plt.loglog(masses, outputs[i], label=f"{age:.3g}")
                plt.ylabel(k)
            plt.xlabel("Mass (solar masses)")
            plt.legend()
