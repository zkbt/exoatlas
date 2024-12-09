from ...imports import *
from ...telescopes import *


@property
def angular_separation(self):
    """
    Calculate the angular separation,
    simply as theta = a/D
    """

    a = self.semimajoraxis
    D = self.distance

    theta = np.arctan(a / D).to(u.arcsec)

    return theta


@property
def imaging_contrast(self):
    """
    What is the reflected light eclipse depth,
    for an albedo of 100%?

    But use a kludged radius
    """
    return 0.25 * (self.kludge_radius / self.semimajoraxis).decompose() ** 2


def transmission_signal(self, mu=2.32, threshold=2):
    """
    What is the transit depth of 1 scale height of an
    atmosphere transiting in front of the star.

    Parameters
    ----------
    mu : float
        Mean molecular weight (default 2.2 for H/He)
    threshold : float
        By how many sigma must the planet mass be detected?

    """
    with np.errstate(invalid="ignore"):
        H = self.scale_height(mu)
        Rp = self.radius
        Rs = self.stellar_radius
        depth = (2 * H * Rp / Rs**2).decompose()

        dlnm = self.get_uncertainty("mass") / self.mass
        bad = dlnm > 1 / threshold
        depth[bad] = np.nan
        return depth


def reflection_signal(self, albedo=0.1):
    """
    What is the reflected light eclipse depth,
    for an albedo of 100%?
    """
    return albedo * 0.25 * (self.radius / self.semimajoraxis).decompose() ** 2


def emission_signal(self, wavelength=5 * u.micron):
    """
    What is the thermal emission eclipse depth,
    assuming Planck spectra for both star and planet?

    This calculation assumes a Bond albedo of 0
    and that heat is uniformly distributed over the planet.

    Parameters
    ----------
    wavelength : astropy.unit.Quantity
        The wavelength at which it should be calculated.
    """

    # create thermal emission sources for both star and planet
    import rainbowconnection as rc

    star = rc.Thermal(teff=self.stellar_teff, radius=self.stellar_radius)
    planet = rc.Thermal(teff=self.teq, radius=self.radius)

    # calculate the depth as the luminosity ratio
    depths = planet.spectrum(wavelength) / star.spectrum(wavelength)

    return depths


def stellar_brightness(self, wavelength=5 * u.micron):
    """
    How many photons/s/m^2/micron do we receive from the star?

    This is calculated from the distance, radius, and
    stellar effective temperature of the stars.

    (It could be potentially be improved with PHOENIX
    model grids and/or cleverness with photometry.)

    Parameters
    ----------
    wavelength : astropy.unit.Quantity
        The wavelength at which it should be calculated.
    """

    # import some tools for easy cartoon spectra
    import rainbowconnection as rc

    # create source with right temperature, size, distance
    teff, radius = self.stellar_teff, self.stellar_radius
    star = rc.Thermal(teff=teff, radius=radius).at(self.distance)

    # calculate the energy flux
    flux_in_energy = star.spectrum(wavelength)

    # convert to photon flux
    photon_energy = con.h * con.c / wavelength / u.ph
    flux_in_photons = flux_in_energy / photon_energy

    # return the
    return flux_in_photons.to("ph s^-1 m^-2 micron^-1")


def stellar_brightness_in_telescope_units(self, telescope_name="JWST", **kw):
    """
    The stellar brightness, converted to telescope units.

    Parameters
    ----------
    telescope_name : str
        The name of the telescope.

    wavelength : astropy.unit.Quantity
        The wavelength at which it should be calculated.

    R : float
        The spectral resolution at which the
        telescope will bin wavelengths.

    dt : astropy.units.quantity.Quantity
        The time over which the telescope exposes.
    """

    # what counts as 1 "telescope unit" (e.g. JWST at R=20 at 5 microns for 1 hour)
    telescope_unit = define_telescope_unit_by_name(telescope_name, **kw)

    # what's the photon flux (photons/m**2/s)
    flux_in_photons = self.stellar_brightness(telescope_unit.wavelength)

    # quote the brightness as (for example) gigaphotons/JWST at R=20 at 5 microns in 1 hour
    unit = lotsofphotons_unit / telescope_unit
    return flux_in_photons.to(unit)


def depth_uncertainty(
    self, telescope_name="JWST", per_transit=False, dt=1 * u.hour, **kw
):
    """
    What is the transit/eclipse depth uncertainty
    with a particular telescope
    at a particular wavelength
    at a particular resolution?

    By default, this will be calculated for one transit.
    Optionally, it can be calculated for a given amount of time instead.

    Parameters
    ----------
    telescope_name : str
        The name of the telescope.

    per_transit : bool
        If True, calculate the depth uncertainty for one transit.
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
    """

    # what counts as 1 "telescope unit" (e.g. JWST at R=20 at 5 microns for 1 hour)
    telescope_unit = define_telescope_unit_by_name(telescope_name, dt=dt, **kw)

    # what's the photon flux (photons/m**2/s)
    flux_in_photons = self.stellar_brightness(telescope_unit.wavelength)

    # what's the total collecting power?
    if per_transit:
        ratio_of_collecting_time = self.transit_duration / dt
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

    return sigma_depth


def _get_noise_and_unit(self, telescope_name="JWST", per_transit=False, **kw):
    """
    Tiny helper to get the noise and the telescope_unit
    for a telescope observation of a planet.
    """

    # figure out the noise
    noise = self.depth_uncertainty(
        telescope_name=telescope_name, per_transit=per_transit, **kw
    )

    # create a telescope unit (mostly to get a default wavelength)
    telescope_unit = define_telescope_unit_by_name(telescope_name, **kw)

    return noise, telescope_unit


def depth_snr(self, telescope_name="JWST", **kw):
    """
    What's the approximate S/N for the detection of the planet's transit?
    """

    noise, telescope_unit = self._get_noise_and_unit(
        telescope_name=telescope_name, **kw
    )
    signal = self.transit_depth
    return signal / noise


def emission_snr(self, telescope_name="JWST", **kw):
    """
    What's the approximate S/N for the detection of the
    thermal emission eclipse of a planet?
    """

    noise, telescope_unit = self._get_noise_and_unit(
        telescope_name=telescope_name, **kw
    )
    signal = self.emission_signal(wavelength=telescope_unit.wavelength)
    return signal / noise


def reflection_snr(self, telescope_name="JWST", albedo=1, **kw):
    """
    What's the approximate S/N for the detection of the
    reflected light eclipse of a planet?
    """

    noise, telescope_unit = self._get_noise_and_unit(
        telescope_name=telescope_name, **kw
    )
    signal = self.reflection_signal(albedo=albedo)
    return signal / noise


def transmission_snr(self, telescope_name="JWST", mu=2.32, threshold=2, **kw):
    """
    What's the approximate S/N for the detection of the
    reflected light eclipse of a planet?
    """

    noise, telescope_unit = self._get_noise_and_unit(
        telescope_name=telescope_name, **kw
    )
    signal = self.transmission_signal(mu=mu, threshold=threshold)
    return signal / noise
