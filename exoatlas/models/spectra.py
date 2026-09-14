# copied and pasted from https://github.com/zkbt/rainbow-connection/
# (trying to avoid a complicated python-version dependency for colour-tools)

from ..imports import *
from .resampling import *

__all__ = ['Spectrum', 'Thermal', 'SumOfThermal']

# define some useful units
spectral_luminosity_unit = u.W / u.micron
spectral_flux_unit = u.W / u.micron / u.m**2
spectral_intensity_unit = u.W / u.micron / u.m**2 / u.sr
luminosity_unit = u.W
flux_unit = u.W / u.m**2
intensity_unit = u.W / u.m**2 / u.sr

unit2name = {
    spectral_flux_unit: "Flux",
    spectral_luminosity_unit: "Luminosity",
    flux_unit: "Flux",
    luminosity_unit: "Luminosity",
    intensity_unit: "Intensity",
    spectral_intensity_unit: "Intensity",
}

def determine_quantity(unit):
    """
    Determine the name of a particular unit.

    Parameters
    ----------
    unit : astropy.units.core.Unit, astropy.units.core.CompositeUnit

    Returns
    -------
    """

    # try all the listed units
    for k in unit2name:
        if unit.is_equivalent(k):
            return unit2name[k]

    # if none match, return "?"
    return "?"


def check_wavelength_unit(w):
    w.to("micron")


bgkw = dict(color="gray", alpha=0.5, zorder=-100)

class Spectrum:
    """
    The Spectrum class is a generic representation of the light
    from some emitting object, particularly for spherically
    symmetric emission.

    By default, the `.spectrum(wavelength)` method will return
    the spectral luminosity of the object. This is a quantity
    with units like W/nm, and it can be integrated over
    wavelength to provide the total luminosity of the
    light-emitting object, in W.

    The `.at(distance)` method creates an object representing
    the spectral flux from the object if seen from some
    particular distance. This is a quantity with units like
    W/nm/m**2, and it can be integrated over wavelength to
    provide the total flux of the light, in W/m**2.

    Classes that inherit from this will likely modify (at least)
    the surface_flux and surface_area methods.
    """

    # the default grid of wavelengths
    default_wavelengths = np.arange(200, 1000) * u.nm

    def __init__(self, wavelength, flux, radius=1 / 4 / np.pi * u.m):
        """
        Initialize a spectrum by providing arrays of
        wavelength and flux. (Normally, some other
        wrapper will be used to create a new Spectrum
        object.)

        Parameters
        ----------
        wavelength : astropy.units.quantity.Quantity
            The wavelength values of the spectrum.
            These must be convertible to units of
            length (nm, m, cm, micron, Angstrom, ...)

        flux : astropy.units.quantity.Quantity
            The flux values of the spectrum. Units
            should be fairly flexible, but W/m**2/nm
            would be pretty reasonable defaults.
        """

        # make sure the wavelengths have some
        check_wavelength_unit(wavelength)

        # assign the hiddgen wavelength and flux values
        self._wavelength = wavelength
        self._flux = flux * u.Unit("")
        self.radius = radius * u.Unit("")

        # set the default wavelengths to be the actual values
        self.wavelength = self._wavelength

    def surface_flux(self, wavelength=None):

        # make sure at least some grid of wavelengths is defined
        w = self.get_wavelength(wavelength)

        original_unit = self._flux.unit
        unitless_flux = self._flux.value

        # bin this spectrum to the particular wavelength grid
        binned = bintogrid(
            x=self._wavelength.to("nm").value,
            y=unitless_flux,
            newx=w.to("nm").value,
            drop_nans=False,
        )
        unitless_neww = binned['x']
        unitless_newf = binned['y']

        # make sure the wavelengths match up
        assert np.all(unitless_neww == w.to("nm").value)

        # make sure the flux units match up
        newf = unitless_newf * original_unit
        assert newf.unit.is_equivalent(self._flux.unit)

        return newf

    def surface_area(self):
        """
        The surface area of the light source,
        in units like m**2.

        Returns
        -------
        surface_area : astropy.units.quantity.Quantity
            The emitting area of the surface, usually in m**2.
        """
        return 4 * np.pi * self.radius**2

    def get_wavelength(self, wavelength=None):
        """
        A wrapper to ensure at least some grid of wavelengths
        gets defined. A default grid will be assumed, unless
        any wavelength array at all is passed.
        """

        # make sure at least some wavelengths are defined
        if wavelength is None:
            wavelength = self.wavelength
        return wavelength.to("micron")

    def spectrum(self, wavelength=None):
        """
        The spectrum of the light source, as spectral luminosity (W/nm)
        or if a distance is defined as spectral flux (W/nm/m**2).

        Parameters
        ----------
        wavelength : astropy.units.quantity.Quantity
            The wavelengths on which we want the spectrum.

        Returns
        -------
        spectrum : astropy.units.quantity.Quantity
            The luminosity (W/nm) or flux (W/nm/m**2).
        """

        # simplify the factor as best we can
        factor = (self.surface_area() / self.normalization()).decompose()

        # return the surface flux with appropriate normalization
        return factor * self.surface_flux(wavelength)

    def normalization(self):
        """
        The normalization by which this Spectrum
        will be divided. It's flexible, to allow
        either spectral luminosity or spectral flux.

        Returns
        -------
        norm : astropy.units.quantity.Quantity
            A normalization
        """
        try:
            # if there's a distance, return a flux
            assert self.distance is not None
            return 4 * np.pi * self.distance**2
        except (AttributeError, AssertionError):
            # by default, simply return a luminosity
            return 1.0

    def angular_size(self):
        """
        The angular size of the light source (if viewed from a distance).
        """

        try:
            # if there's a distance, return an angular size
            assert self.distance is not None
            return np.arctan(self.radius / self.distance).to("deg")
        except (AttributeError, AssertionError):
            # complain if no distance is defined
            raise ValueError(
                """
            This Spectrum has no .distance attribute.
            Please consider using `.at(distance)` to
            create a new light source as viewed from
            a distance.
            """
            )

    def at(self, distance=1 * u.au):
        """
        Create a new Spectrum representing the current
        light source viewed from some distance
        (assuming spherical symmetry).

        Parameters
        ----------
        distance : astropy.units.quantity.Quantity
            The distance at which we're viewing this source.

        Returns
        -------
        flux : Spectrum
            A new Spectrum, with the distance attached.
        """

        # create a copy of the current spectrum
        new = copy.deepcopy(self)

        # update this copy's distance and return
        new.distance = distance
        return new

    def filter(self, f):
        """
        Create a new Spectrum representing the current
        light source viewed through some filter.

        Parameters
        ----------
        f : function
            The filter transmission function. This function
            must take wavelength (with units) as an input.

        Returns
        -------
        spectrum : Spectrum
            A new Spectrum, with the distance attached.
        """

        # create a copy of the current spectrum
        new = copy.deepcopy(self)

        # update this copy's distance and return
        new.filter = f
        try:
            new._flux = self._flux*f(self._wavelength)
        except AttributeError:
            def surface_flux(self, wavelength=None):
                w = self.get_wavelength(wavelength)
                unfiltered_flux = self.surface_flux(wavelength=w)
                filtered_flux = unfiltered_flux*f(w)
                return filtered_flux
            new.surface_flux = surface_flux
        return new #BLERG! DOESN'T WORK!



    # FIXME -- for analytic functions, it'd help to define some kind of
    # a bounding box in wavelength space, so this integral could be done
    # analytically or with scipy.integrate.quad
    def integrate(self, lower=None, upper=None):
        """
        Integrate the spectrum over wavelength.

        It gives a number with units identical to the results of
        `.spectrum()` but without the wavelength (W or W/m**2).

        Parameters
        ----------
        lower : astropy.units.quantity.Quantity
            The lower wavelength limit.

        upper : astropy.units.quantity.Quantity
            The lower wavelength limit.

        Returns
        -------
        integral : astropy.units.quantity.Quantity
            The integral over wavelength
        """

        w = self.wavelength
        f = self.spectrum(w)

        ok = np.ones(np.shape(w)).astype(bool)
        if lower is not None:
            ok *= w >= lower

        if upper is not None:
            ok *= w <= upper
            # raise NotImplementedError('Wavelength limits not yet OK.')

        return np.trapezoid(f[ok], w[ok])

        # np.trapezoid(f.value, w.value)*f.unit*w.unit
        # wlower = lower or self.wavelength[0]
        # wupper = upper or self.wavelength[-1]
        # return quad(self.spectrum, wlower, wupper)

    def __repr__(self):
        """
        How should this object appear as a string?

        Returns
        -------
        s : str
            A simple string representation.
        """
        try:
            assert self.distance is not None
            return f"{self.__class__.__name__} at {self.distance}"
        except (AssertionError, AttributeError):
            return f"{self.__class__.__name__}"

    def set_power(self, power=100 * u.W):
        """
        Change the radius of the object so that it will
        emit a specified power, with the same spectral shape.

        Parameters
        ----------
        power : astropy.units.quantity.Quantity
            The total power we want the object to emit.
        """

        # calculate the total luminosity of this object
        total = self.integrate()

        # make sure we're dealing with an actual luminosity
        assert total.unit.is_equivalent("W")

        # calculation a new normalization
        normalization = (power / total).decompose()

        # change the radius of this object
        self.radius *= np.sqrt(normalization)
        self.power = power

    def mean_intensity(self, wavelength=None):
        """
        Calculate the mean intensity field created by the source.

        The mean intensity field represents the intensity
        of the source, smeared over a full 4pi steradians.

        (This makes sense only for spectra viewed from a distance.)
        """

        # what is the intensity of the disk
        F_disk = self.spectrum(wavelength)

        # make sure we're dealing with a flux
        assert F_disk.unit.is_equivalent(u.W / u.m**2 / u.nm)

        # calculate the mean intensity
        solid_angle = np.pi * self.angular_size() ** 2

        # calculate the mean intensity field
        J = F_disk / 4 / np.pi / u.sr

        return J

    def disk_intensity(self, wavelength=None):
        """
        Calculate the intensity of the disk of the source.

        The disk intensity represents the intensity of staring
        directly at the disk. The flux of the source is its
        intensity integrated over solid angle, so two sources
        with the same intensities can have very different
        fluxes, based on their apparent angular sizes.
        For example, the Sun seen `.at` different distances
        will have different fluxes but the same disk intensity.

        (This makes sense only for spectra viewed from a distance.)
        """

        # what is the intensity of the disk
        F_disk = self.spectrum(wavelength)

        # make sure we're dealing with a flux
        assert F_disk.unit.is_equivalent(u.W / u.m**2 / u.nm)

        # calculate the mean intensity
        solid_angle = np.pi * self.angular_size() ** 2
        I_disk = F_disk / solid_angle

        return I_disk


class Thermal(Spectrum):
    def __init__(self, teff=5800 * u.K, radius=1 * u.Rsun):

        self.teff = teff
        self.radius = radius
        self.wavelength = np.logspace(2, 3, 1000) * u.nm

    def intensity(self, wavelength):
        """
        This function calculates the thermal emission intensity spectrum of a surface.

            Inputs:
                wavelength = numpy array of wavelengths (with astropy units)

            Outputs:
                Returns an array of thermal emission intensities,
                in astropy units of W/(m^2*micron*sr). This is a flux, which has
                already been integrated over solid angle.
        """

        temperature = self.teff

        # define variables as shortcut to the constants we need
        h = con.h
        k = con.k_B
        c = con.c

        # this is the thing that goes into the exponent (it's units better cancel!)
        up = h * c / (wavelength * k * temperature)

        # calculate the intensity from the Planck function
        intensity = (2 * h * c ** 2 / wavelength ** 5 / (np.exp(up) - 1)) / u.steradian

        # return the intensity
        return intensity.to("W/(m**2*micron*sr)")

    def surface_flux(self, wavelength):
        """
        This function calculates the thermal emission flux spectrum of a surface.

            Inputs:
                wavelength = numpy array of wavelengths (with astropy units)
                temperature = a single number, the temperature (with astropy units)

            Outputs:
                Returns an array of thermal emission fluxes,
                in astropy units of W/(m^2*micron). This is a flux, which has
                already been integrated over solid angle.
        """

        # calculate the flux, knowing the angle integral will be pi steradians (for isotropic emission)
        flux = self.intensity(wavelength) * np.pi * u.steradian

        # return the flux, in convenient units
        return flux.to("W/(micron * m**2)")

    def __repr__(self):
        """
        How should this object appear as a string?

        Returns
        -------
        s : str
            A simple string representation.
        """

        basic = f"{self.__class__.__name__} ({self.teff:.0f}, {self.radius})"
        try:
            assert self.distance is not None
            return basic + f" at {self.distance}"
        except (AssertionError, AttributeError):
            return basic

    def integrate(self, lower=None, upper=None):
        """
        Integrate the spectrum over wavelength.

        It gives a number with units identical to the results of
        `.spectrum()` but without the wavelength (W or W/m**2).

        Parameters
        ----------
        lower : astropy.units.quantity.Quantity
            The lower wavelength limit.

        upper : astropy.units.quantity.Quantity
            The lower wavelength limit.

        Returns
        -------
        integral : astropy.units.quantity.Quantity
            The integral over wavelength.
        """

        # if wavelength limits are used, revert back to the numerical integral
        if (lower is not None) or (upper is not None):
            return super().integrate(lower=lower, upper=upper)

        # if there are infinite wavelength limits, do the integral analytically
        surface_flux = con.sigma_sb * self.teff ** 4
        factor = (self.surface_area() / self.normalization()).decompose()
        return factor * surface_flux


class SumOfThermal(Thermal):
    def intensity(self, wavelength, n_grid=90, plots=False):
        # THANKS PAT!!

        temperature = self.teff
        T_dmax = temperature

        # define variables as shortcut to the constants we need
        h = con.h
        k = con.k_B
        c = con.c

        # define a 2D grid of latitude and longitude
        dlat = np.pi / n_grid
        dlon = np.pi / n_grid

        lats = np.arange(-np.pi / 2., np.pi / 2., dlat)
        lons = np.arange(-np.pi, np.pi, dlon)
        lats_mesh, lons_mesh = np.meshgrid(lats, lons, indexing='ij')

        # define inclination
        inc = np.pi / 2.

        # calculate substellar temperature (f=1)
        T_sub = T_dmax * (3. / 2.) ** (1. / 4.)

        # calculate cosine from substellar and temperature at each point in 2D map
        cos_angle_hotspot = np.cos(lats_mesh) * np.cos(lons_mesh)
        T_points = T_sub * cos_angle_hotspot ** 0.25
        T_points[cos_angle_hotspot <= 0] = 0  # set nightside to 0
        if plots:
            plt.imshow(T_points.value)
            plt.colorbar()
            plt.show()

        # calculate spectrum at each point
        planet_spectrum = (2 * h * c ** 2 / wavelength ** 5 / (
                    np.exp(h * c / (wavelength * k * T_points[:, :, np.newaxis])) - 1)) / u.steradian
        cos_angle_obs = np.sin(lats_mesh) * np.sin(np.pi / 2 - inc) + np.cos(lats_mesh) * np.cos(
            np.pi / 2 - inc) * np.cos(lons_mesh)
        planet_spectrum[cos_angle_obs <= 0] = 0
        unit_solid_angle = np.cos(lats_mesh) * dlat * dlon

        # calculate intensity and sum them
        intensity = planet_spectrum * (unit_solid_angle * cos_angle_obs)[:, :, np.newaxis]
        sum_intensity = intensity.sum(axis=(0, 1)) / np.pi
        if plots:
            plt.imshow(intensity[:, :, 0].value)
            plt.colorbar()
            plt.show()

        return sum_intensity.to("W/(m**2*nm*sr)")