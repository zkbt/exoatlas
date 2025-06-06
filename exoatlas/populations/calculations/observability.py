from ...imports import *
from ...telescopes import *


def angular_separation(self, distribution=False, **kw):
    """
    Maximum Angular Separation (a/D, arcseconds)

    Calculate a representative maximum angular separation
    from the semimajor axis and distance.

    TO-DO:
    Include eccentricity.

    Parameters
    ----------
    distribution : bool
        If False, return a simple array of values.
        If True, return an astropy.uncertainty.Distribution,
        which can be used for error propagation.
    """

    a = self.semimajoraxis(distribution=distribution)
    D = self.distance(distribution=distribution)

    theta = np.arctan(a / D).to(u.arcsec)

    return theta


def imaging_contrast(self, albedo=2 / 3, phase_function=0.25, distribution=False, **kw):
    """
    Reflected Light Imaging Contrast (unitless)

    Calculate an estimated imaging contrast for reflected light.

    (FIXME - check math + explanation on albedo + phase_function!)

    Parameters
    ----------
    albedo : float
        What's the geometric albedo? Default is 2/3, which (I think)
        is reasonably the maximum for Lambertian scattering from
        a face-on planet.
    phase_function : float
        What fraction of reflected light comes toward us, compared
        to the immediate dayside reflection from direct illumination.
    distribution : bool
        If False, return a simple array of values.
        If True, return an astropy.uncertainty.Distribution,
        which can be used for error propagation.
    """

    return (
        phase_function
        * albedo
        * 0.25
        * (
            self.kludge_radius(distribution=distribution)
            / self.semimajoraxis(distribution=distribution)
        ).decompose()
        ** 2
    )


def transmission_signal(self, mu=2.3, kludge=False, distribution=False, **kw):
    """
    Transmission Transit Signal (2*H*Rp/Rs**2, unitless)

    Calculate an estimate of the transit depth of
    one atmospheric scale height transiting in front
    of the star.

    TO-DO:
    Should we be fussy about allowing other equilibrium
    temperatures for calculating the scale height?

    Parameters
    ----------
    mu : float
        Mean molecular weight (default 2.3 for something
        like solar composition in chemical equilibrium)
    kludge : bool
        Should we include kludged estimates for mass (from msini and/or
        empirical mass-radius) and/or radius (from empircal mass-radius)
        when doing this calculation?
    distribution : bool
        If False, return a simple array of values.
        If True, return an astropy.uncertainty.Distribution,
        which can be used for error propagation.
    """
    H = self.scale_height(mu=mu, kludge=kludge, distribution=distribution)
    Rp = self.radius(distribution=distribution)
    Rs = self.stellar_radius(distribution=distribution)
    depth = (2 * H * Rp / Rs**2).decompose()
    return depth


def reflection_signal(self, albedo=0.1, distribution=False, **kw):
    """
    Reflection Eclipse Signal (unitless)

    Calculate an estimate of the reflected light eclipse depth
    for a particular albedo.

    TO-DO:
    Check reflected light math!!

    Parameters
    ----------
    albedo : float
        How much incoming light is reflected away to space?
        Default of 0.1 might be reasonable for rocky surfaces
        and/or very absorbing hot Jupiters?
    distribution : bool
        If False, return a simple array of values.
        If True, return an astropy.uncertainty.Distribution,
        which can be used for error propagation.
    """
    Rp = self.radius(distribution=distribution)
    a = self.semimajoraxis(distribution=distribution)
    return albedo * 0.25 * (Rp / a).decompose() ** 2


def emission_signal(
    self, wavelength=5 * u.micron, albedo=0, f=1 / 4, distribution=False, **kw
):
    """
    Thermal Emission Eclipse Signal (unitless)

    Calculate an estimate of the thermal emission eclipse depth
    at a particular wavelength, for a planet's equilibrium
    temperature calculated for a particular albedo and f.
    Calculations (falsely!) assume Planck spectra for both
    the star and planet. Default parameters assumed 0 albedo
    and that heat is redistributed uniformly over the entire
    sphere of the planet.

    Parameters
    ----------
    wavelength : astropy.unit.Quantity
        The wavelength at which it should be calculated.
    albedo : float
        Fraction of light reflected instead of emitted.
        See .teq() docstring for details + interpretation.
    f : float
        Heat distribution parameter for dayside temperature.
        See .teq() docstring for details and interpretation.
    distribution : bool
        If False, return a simple array of values.
        If True, return an astropy.uncertainty.Distribution,
        which can be used for error propagation.
    """

    # create thermal emission sources for both star and planet
    import rainbowconnection as rc

    star = rc.Thermal(
        teff=self.stellar_teff(distribution=distribution),
        radius=self.stellar_radius(distribution=distribution),
    )
    planet = rc.Thermal(
        teff=self.teq(albedo=albedo, f=f, distribution=distribution),
        radius=self.radius(distribution=distribution),
    )

    # calculate the depth as the luminosity ratio
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        depths = planet.spectrum(wavelength) / star.spectrum(wavelength)

    return depths


def stellar_brightness(self, wavelength=5 * u.micron, distribution=False, **kw):
    """
    Stellar Brightness (photons/s/m^2/micron)

    Calculate the brightness of the star as observed from the Earth,
    at a particular wavelength assuming Planck thermal emission,
    based on the star's distance, radius, and effective temperature.

    TO-DO:
    Check consistency between Rs/Teff and Ls.
    Consider implementing model spectra grids?!

    Parameters
    ----------
    wavelength : astropy.unit.Quantity
        The wavelength at which it should be calculated.
    distribution : bool
        If False, return a simple array of values.
        If True, return an astropy.uncertainty.Distribution,
        which can be used for error propagation.
    """

    # import some tools for easy cartoon spectra
    import rainbowconnection as rc

    # create source with right temperature, size, distance
    Teff = self.stellar_teff(distribution=distribution)
    Rs = self.stellar_radius(distribution=distribution)
    D = self.distance(distribution=distribution)
    star = rc.Thermal(teff=Teff, radius=Rs).at(distance=D)

    # calculate the energy flux
    flux_in_energy = star.spectrum(wavelength)

    # convert to photon flux
    photon_energy = con.h * con.c / wavelength / u.ph
    flux_in_photons = flux_in_energy / photon_energy

    # return photon flux in good units
    return flux_in_photons.to("ph s^-1 m^-2 micron^-1")


def stellar_brightness_in_telescope_units(
    self, telescope_name="JWST", distribution=False, **kw
):
    """
    Stellar Brightness (photons/telescope/hour)

    Calculate the brightness of the star as observed from the Earth,
    converted into units more relevant to a particular telescope.

    Parameters
    ----------
    wavelength : astropy.unit.Quantity
        The wavelength at which it should be calculated.
    R : float
        The spectral resolution at which the
        telescope will bin wavelengths.
    dt : astropy.units.quantity.Quantity
        The time over which the telescope exposes.
    distribution : bool
        If False, return a simple array of values.
        If True, return an astropy.uncertainty.Distribution,
        which can be used for error propagation.
    """

    # what counts as 1 "telescope unit" (e.g. JWST at R=20 at 5 microns for 1 hour)
    telescope_unit = define_telescope_unit_by_name(telescope_name, **kw)

    # what's the photon flux (photons/m**2/s)
    flux_in_photons = self.stellar_brightness(
        wavelength=telescope_unit.wavelength, distribution=distribution
    )

    # quote the brightness as (for example) gigaphotons/JWST at R=20 at 5 microns in 1 hour
    unit = lotsofphotons_unit / telescope_unit
    return flux_in_photons.to(unit)


def depth_uncertainty(
    self,
    telescope_name="JWST",
    per_transit=False,
    dt=1 * u.hour,
    distribution=False,
    **kw
):
    """
    Transit/Eclipse Depth Uncertainty (unitless)

    Calculate the estimated uncertainty expected for transit
    or eclipse depth measurements with a particular telescope.
    The depth uncertainty assumes photon-limited noise based on
    the telescope (with implied effective collecting area),
    the wavelength (which sets the stellar brightness),
    the spectral resolution (which sets the bandpass width),
    and the time being observed.

    By default, depth uncertainties will be calculated per constant time;
    optionally, depth uncertainty can be calculated per transit,
    where longer transits will have relatively smaller uncertainty.


    Parameters
    ----------
    telescope_name : str
        The name of the telescope.
    per_transit : bool
        If True, calculate the depth uncertainty for one transit,
        based on the calculated transit duration for the planet.
        If False, calculate the depth uncertainty for a certain amount
        of in-transit time. You likely want to specify `dt` as a
        keyword argument to set that amount of in-transit time.
        In either case, an out-of-transit baseline equal to the
        total in-transit time will be assumed. This means the actual
        time cost will be twice the transit duration or `dt` chosen,
        and the depth uncertainty will be a factor sqrt(2) larger
        than the pure photon noise binned to the relevant timescale.
    wavelength : astropy.unit.Quantity
        The wavelength at which it should be calculated.
    R : float
        The spectral resolution at which the
        telescope will bin wavelengths.
    dt : astropy.units.quantity.Quantity
        The time over which the telescope exposes. If `per_transit=True`,
        this will be ignored. Otherwise, it will set the total amount
        of in-transit time observed, assuming that an equal amount of
        time will *also* be observed out of transit.
    distribution : bool
        If False, return a simple array of values.
        If True, return an astropy.uncertainty.Distribution,
        which can be used for error propagation.
    """

    # what counts as 1 "telescope unit" (e.g. JWST at R=20 at 5 microns for 1 hour)
    telescope_unit = define_telescope_unit_by_name(telescope_name, dt=dt, **kw)

    # what's the photon flux (photons/m**2/s)
    flux_in_photons = self.stellar_brightness(
        wavelength=telescope_unit.wavelength, distribution=distribution
    )

    # what's the total collecting power?
    if per_transit:
        ratio_of_collecting_time = self.transit_duration(distribution=distribution) / dt
    else:
        ratio_of_collecting_time = 1
    collecting_power = 1 * telescope_unit * ratio_of_collecting_time

    # what's the total number of photons collected during transit
    N = (flux_in_photons * collecting_power).to(u.ph).value

    # what's the flux uncertainty on the time scale of one transit?
    sigma = 1 / np.sqrt(N)

    # inflate by a factor of sqrt(2) for equal out-of-transit
    oot = np.sqrt(2)
    sigma_depth = sigma * oot

    return u.Quantity(sigma_depth)


def _get_noise_and_unit(
    self, telescope_name="JWST", per_transit=False, distribution=False, **kw
):
    """
    Tiny helper to get the noise and the telescope_unit
    for a telescope observation of a planet.

    Parameters
    ----------
    telescope_name : str
        The name of the telescope.
    per_transit : bool
        If True, calculate the depth uncertainty for one transit,
        based on the calculated transit duration for the planet.
        If False, calculate the depth uncertainty for a certain amount
        of in-transit time. You likely want to specify `dt` as a
        keyword argument to set that amount of in-transit time.
        In either case, an out-of-transit baseline equal to the
        total in-transit time will be assumed. This means the actual
        time cost will be twice the transit duration or `dt` chosen,
        and the depth uncertainty will be a factor sqrt(2) larger
        than the pure photon noise binned to the relevant timescale.
    wavelength : astropy.unit.Quantity
        The wavelength at which it should be calculated.
    R : float
        The spectral resolution at which the
        telescope will bin wavelengths.
    dt : astropy.units.quantity.Quantity
        The time over which the telescope exposes. If `per_transit=True`,
        this will be ignored. Otherwise, it will set the total amount
        of in-transit time observed, assuming that an equal amount of
        time will *also* be observed out of transit.
    distribution : bool
        If False, return a simple array of values.
        If True, return an astropy.uncertainty.Distribution,
        which can be used for error propagation.
    """

    # figure out the noise
    noise = self.depth_uncertainty(
        telescope_name=telescope_name,
        per_transit=per_transit,
        distribution=distribution,
        **kw
    )

    # create a telescope unit (mostly to get a default wavelength)
    telescope_unit = define_telescope_unit_by_name(telescope_name, **kw)

    return noise, telescope_unit


def depth_snr(self, telescope_name="JWST", distribution=False, **kw):
    """
    Transit Depth S/N

    Calculate the S/N for detecting the transit
    with a particular telescope.

    Parameters
    ----------
    telescope_name : str
        The name of the telescope.
    per_transit : bool
        If True, calculate the depth uncertainty for one transit,
        based on the calculated transit duration for the planet.
        If False, calculate the depth uncertainty for a certain amount
        of in-transit time. You likely want to specify `dt` as a
        keyword argument to set that amount of in-transit time.
        In either case, an out-of-transit baseline equal to the
        total in-transit time will be assumed. This means the actual
        time cost will be twice the transit duration or `dt` chosen,
        and the depth uncertainty will be a factor sqrt(2) larger
        than the pure photon noise binned to the relevant timescale.
    wavelength : astropy.unit.Quantity
        The wavelength at which it should be calculated.
    R : float
        The spectral resolution at which the
        telescope will bin wavelengths.
    dt : astropy.units.quantity.Quantity
        The time over which the telescope exposes. If `per_transit=True`,
        this will be ignored. Otherwise, it will set the total amount
        of in-transit time observed, assuming that an equal amount of
        time will *also* be observed out of transit.
    distribution : bool
        If False, return a simple array of values.
        If True, return an astropy.uncertainty.Distribution,
        which can be used for error propagation.
    """

    noise, telescope_unit = self._get_noise_and_unit(
        telescope_name=telescope_name, distribution=distribution, **kw
    )
    signal = self.transit_depth(distribution=distribution)
    return signal / noise


def emission_snr(
    self, telescope_name="JWST", albedo=0, f=1 / 4, distribution=False, **kw
):
    """
    Thermal Emission Eclipse Depth S/N

    Calculate the S/N for detecting the
    planet's thermal emission eclipse
    with a particular telescope.

    Parameters
    ----------
    telescope_name : str
        The name of the telescope.
    per_transit : bool
        If True, calculate the depth uncertainty for one transit,
        based on the calculated transit duration for the planet.
        If False, calculate the depth uncertainty for a certain amount
        of in-transit time. You likely want to specify `dt` as a
        keyword argument to set that amount of in-transit time.
        In either case, an out-of-transit baseline equal to the
        total in-transit time will be assumed. This means the actual
        time cost will be twice the transit duration or `dt` chosen,
        and the depth uncertainty will be a factor sqrt(2) larger
        than the pure photon noise binned to the relevant timescale.
    wavelength : astropy.unit.Quantity
        The wavelength at which it should be calculated.
    R : float
        The spectral resolution at which the
        telescope will bin wavelengths.
    dt : astropy.units.quantity.Quantity
        The time over which the telescope exposes. If `per_transit=True`,
        this will be ignored. Otherwise, it will set the total amount
        of in-transit time observed, assuming that an equal amount of
        time will *also* be observed out of transit.
    albedo : float
        Fraction of light reflected instead of emitted.
        See .teq() docstring for details + interpretation.
    f : float
        Heat distribution parameter for dayside temperature.
        See .teq() docstring for details and interpretation.
    distribution : bool
        If False, return a simple array of values.
        If True, return an astropy.uncertainty.Distribution,
        which can be used for error propagation.
    """

    noise, telescope_unit = self._get_noise_and_unit(
        telescope_name=telescope_name, distribution=distribution, **kw
    )
    kw["wavelength"] = telescope_unit.wavelength
    signal = self.emission_signal(albedo=albedo, f=f, distribution=distribution, **kw)
    return signal / noise


def reflection_snr(self, telescope_name="JWST", albedo=0.1, distribution=False, **kw):
    """
    Reflected Light Eclipse Depth S/N

    Calculate the S/N for detecting the
    planet's reflected light eclipse
    with a particular telescope.

    Parameters
    ----------
    telescope_name : str
        The name of the telescope.
    per_transit : bool
        If True, calculate the depth uncertainty for one transit,
        based on the calculated transit duration for the planet.
        If False, calculate the depth uncertainty for a certain amount
        of in-transit time. You likely want to specify `dt` as a
        keyword argument to set that amount of in-transit time.
        In either case, an out-of-transit baseline equal to the
        total in-transit time will be assumed. This means the actual
        time cost will be twice the transit duration or `dt` chosen,
        and the depth uncertainty will be a factor sqrt(2) larger
        than the pure photon noise binned to the relevant timescale.
    wavelength : astropy.unit.Quantity
        The wavelength at which it should be calculated.
    R : float
        The spectral resolution at which the
        telescope will bin wavelengths.
    dt : astropy.units.quantity.Quantity
        The time over which the telescope exposes. If `per_transit=True`,
        this will be ignored. Otherwise, it will set the total amount
        of in-transit time observed, assuming that an equal amount of
        time will *also* be observed out of transit.
    albedo : float
        Fraction of light reflected instead of emitted.
        See .teq() docstring for details + interpretation.
    f : float
        Heat distribution parameter for dayside temperature.
        See .teq() docstring for details and interpretation.
    albedo : float
        How much incoming light is reflected away to space?
        Default of 0.1 might be reasonable for rocky surfaces
        and/or very absorbing hot Jupiters?
    distribution : bool
        If False, return a simple array of values.
        If True, return an astropy.uncertainty.Distribution,
        which can be used for error propagation.
    """

    noise, telescope_unit = self._get_noise_and_unit(
        telescope_name=telescope_name, distribution=distribution, **kw
    )
    signal = self.reflection_signal(albedo=albedo, distribution=distribution)
    return signal / noise


def transmission_snr(
    self, telescope_name="JWST", mu=2.3, kludge=False, distribution=False, **kw
):
    """
    Transmission Spectroscopy Transit Depth S/N

    Calculate the S/N for detecting the planet's
    transmission spectrum change in transit depths
    with a particular telescope.

    Parameters
    ----------
    telescope_name : str
        The name of the telescope.
    per_transit : bool
        If True, calculate the depth uncertainty for one transit,
        based on the calculated transit duration for the planet.
        If False, calculate the depth uncertainty for a certain amount
        of in-transit time. You likely want to specify `dt` as a
        keyword argument to set that amount of in-transit time.
        In either case, an out-of-transit baseline equal to the
        total in-transit time will be assumed. This means the actual
        time cost will be twice the transit duration or `dt` chosen,
        and the depth uncertainty will be a factor sqrt(2) larger
        than the pure photon noise binned to the relevant timescale.
    wavelength : astropy.unit.Quantity
        The wavelength at which it should be calculated.
    R : float
        The spectral resolution at which the
        telescope will bin wavelengths.
    dt : astropy.units.quantity.Quantity
        The time over which the telescope exposes. If `per_transit=True`,
        this will be ignored. Otherwise, it will set the total amount
        of in-transit time observed, assuming that an equal amount of
        time will *also* be observed out of transit.
    albedo : float
        Fraction of light reflected instead of emitted.
        See .teq() docstring for details + interpretation.
    f : float
        Heat distribution parameter for dayside temperature.
        See .teq() docstring for details and interpretation.
    mu : float
        Mean molecular weight (default 2.3 for something
        like solar composition in chemical equilibrium)
    kludge : bool
        Should we include kludged estimates for mass (from msini and/or
        empirical mass-radius) and/or radius (from empircal mass-radius)
        when doing this calculation?
    distribution : bool
        If False, return a simple array of values.
        If True, return an astropy.uncertainty.Distribution,
        which can be used for error propagation.
    """

    noise, telescope_unit = self._get_noise_and_unit(
        telescope_name=telescope_name, distribution=distribution, **kw
    )
    signal = self.transmission_signal(mu=mu, kludge=kludge, distribution=distribution)
    return signal / noise
