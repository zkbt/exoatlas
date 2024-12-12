from ...imports import *


def semimajoraxis_from_period(self, distribution=False, **kw):
    """
    Planet Semi-major Axis (a, AU)

    Calculate "a" from the period and stellar mass.
    This might be used if no semi-major axis is
    available in the standardized table.

    Parameters
    ----------
    distribution : bool
        If False, return a simple array of values.
        If True, return an astropy.uncertainty.Distribution,
        which can be used for error propagation.
    """
    P = self.get("period", distribution=distribution)
    M = self.get("stellar_mass", distribution=distribution)
    G = con.G
    a = ((G * M * P**2 / 4 / np.pi**2) ** (1 / 3)).to("AU")
    return a


def semimajoraxis_from_transit_scaled_semimajoraxis(self, distribution=False, **kw):
    """
    Planet Semi-major Axis (a, AU)

    Calculate "a" from the transit-derived a/Rs.
    This might be used if no semi-major axis is
    available in the standardized table.

    Parameters
    ----------
    distribution : bool
        If False, return a simple array of values.
        If True, return an astropy.uncertainty.Distribution,
        which can be used for error propagation.
    """
    a_over_R = self.get("transit_scaled_semimajoraxis", distribution=distribution)
    R = self.get("stellar_radius", distribution=distribution)
    a = (a_over_R * R).to("AU")
    return a


def semimajoraxis(self, distribution=False, **kw):
    """
    Planet Semi-major Axis (a, AU)

    Retrieve "a" first from the standardized table,
    then from period/mass, then from transit a/R.

    Parameters
    ----------
    distribution : bool
        If False, return a simple array of values.
        If True, return an astropy.uncertainty.Distribution,
        which can be used for error propagation.
    """
    a = self._choose_calculation(
        methods=[
            "semimajoraxis_from_table",
            "semimajoraxis_from_period",
            "semimajoraxis_from_transit_scaled_semimajoraxis",
        ],
        distribution=distribution,
        **kw,
    )
    return a


def scaled_semimajoraxis_from_semimajoraxis(self, distribution=False, **kw):
    """
    Planet Scaled Semi-major Axis (a/R*, unitless)

    Calculate "a/Rs" from the semimajor axis a and
    radius Rs. This might be used if the table has
    no direct transit-derived value in its table.

    Parameters
    ----------
    distribution : bool
        If False, return a simple array of values.
        If True, return an astropy.uncertainty.Distribution,
        which can be used for error propagation.
    """
    a = self.get("semimajoraxis", distribution=distribution)
    R = self.get("stellar_radius", distribution=distribution)
    return (a / R).decompose()


def scaled_semimajoraxis(self, distribution=False, **kw):
    """
    Planet Scaled Semi-major Axis (a/R*, unitless)

    Retrieve "a/Rs" first from the transit-derived value
    in the standardized table, then from a and Rs.

    Parameters
    ----------
    distribution : bool
        If False, return a simple array of values.
        If True, return an astropy.uncertainty.Distribution,
        which can be used for error propagation.
    """
    a_over_Rs = self._choose_calculation(
        methods=[
            "transit_scaled_semimajoraxis_from_table",
            "scaled_semimajoraxis_from_semimajoraxis",
        ],
        distribution=distribution,
        **kw,
    )
    return a_over_Rs


def eccentricity(self, distribution=False, **kw):
    """
    Planet Orbital Eccentricity (e, unitless)

    The eccentricity of the planet's orbit. If eccentricity is
    `nan` in the standardized table, this method will quietly
    assume it to be zero.

    Parameters
    ----------
    distribution : bool
        If False, return a simple array of values.
        If True, return an astropy.uncertainty.Distribution,
        which can be used for error propagation.
    """

    # pull out the actual values from the table
    e = self.get_values_from_table("eccentricity", distribution=distribution)

    # replace nans with 0
    if distribution == False:
        bad = np.isfinite(e) == False
        e[bad] = 0
    # (by not doing anything to replace nan values with 0 if `distribution==True`,
    # planets missing original values will end up with `nan` uncertainties)

    return e


def argument_of_periastron(self, distribution=False, **kw):
    """
    Planet Orbital Argument of Periastron ($\omega$, degrees)

    The argument of periastron of the planet's orbit, $\omega$.
    If it is `nan` (often because eccentricity is zero), then
    this will be assumed to be 0.

    Parameters
    ----------
    distribution : bool
        If False, return a simple array of values.
        If True, return an astropy.uncertainty.Distribution,
        which can be used for error propagation.
    """
    # pull out the actual values from the table
    argument_of_periastron = self.standard["argument_of_periastron"].copy()

    # replace nans with 0
    if distribution == False:
        bad = np.isfinite(argument_of_periastron) == False
        argument_of_periastron[bad] = 0 * u.deg
    # (by not doing anything to replace nan values with 0 if `distribution==True`,
    # planets missing original values will end up with `nan` uncertainties)

    return argument_of_periastron


def transit_impact_parameter_from_inclination(self, distribution=False, **kw):
    """
    Planet Semi-major Axis (a, AU)

    Calculate "b" from the inclination, scaled semimajor axis,
    eccentricity, and argument of periastron. This might be
    used if no direct transit-derived impact parameter is
    available in the standardized table.

    Parameters
    ----------
    distribution : bool
        If False, return a simple array of values.
        If True, return an astropy.uncertainty.Distribution,
        which can be used for error propagation.
    """
    # extract necessary quantities
    a_over_Rs = self.scaled_semimajoraxis(distribution=distribution)
    i = self.inclination(distribution=distribution)
    e = self.eccentricity(distribution=distribution)
    omega = self.argument_of_periastron(distribution=distribution)

    # calculate impact parameter based on instantaneous distance at transit
    b = a_over_Rs * np.cos(i) * ((1 - e**2) / (1 + e * np.sin(omega)))
    return b


def transit_impact_parameter(self, distribution=False, **kw):
    """
    Planet Impact Parameter (b)

    Retrieve "b" first from the transit-derived value
    in the standardized table, then from inclination.

    Parameters
    ----------
    distribution : bool
        If False, return a simple array of values.
        If True, return an astropy.uncertainty.Distribution,
        which can be used for error propagation.
    """
    b = self._choose_calculation(
        methods=[
            "transit_impact_parameter_from_table",
            "transit_impact_parameter_from_inclination",
        ],
        distribution=distribution,
        **kw,
    )
    return b


# the 1360 W/m^2 that Earth receives from the Sun
earth_insolation = (1 * u.Lsun / 4 / np.pi / u.AU**2).to(u.W / u.m**2)


@property
def insolation(self, distribution=False, **kw):
    """
    Planet Insolation (S, W/m**2)

    Calculate the insolation the planet receives from its star,
    given the luminosity of the star and the semimajor axis,
    expressed in units of W/m**2. (For reference, Earth
    receives 1360 W/m**2).

    Parameters
    ----------
    distribution : bool
        If False, return a simple array of values.
        If True, return an astropy.uncertainty.Distribution,
        which can be used for error propagation.
    """

    # calculate the average insolation the planet receives
    L = self.stellar_luminosity(distribution=distribution)
    a = self.semimajoraxis(distribution=distribution)
    S = L / 4 / np.pi / a**2
    return S.to(u.W / u.m**2)


def relative_insolation(self, distribution=False, **kw):
    """
    Relative Planet Insolation  (S/S_Earth)

    Calculate the insolation the planet receives from its star,
    given the luminosity of the star and the semimajor axis,
    expressed relative to Earth's insolation.

    Parameters
    ----------
    distribution : bool
        If False, return a simple array of values.
        If True, return an astropy.uncertainty.Distribution,
        which can be used for error propagation.
    """
    return self.insolation(distribution=distribution) / earth_insolation


def log_relative_insolation(self, distribution=False, **kw):
    """
    log(Relative Planet Insolation)

    Calculate log10 of the insolation the planet receives from its star,
    given the luminosity of the star and the semimajor axis,
    expressed relative to Earth's insolation.

    Parameters
    ----------
    distribution : bool
        If False, return a simple array of values.
        If True, return an astropy.uncertainty.Distribution,
        which can be used for error propagation.
    """
    return np.log10(self.relative_insolation(distribution=distribution))


def relative_cumulative_xuv_insolation(self, distribution=False, **kw):
    """
    Relative Cumulative XUV Insolation

    Calculate a *very very very* approximate estimate for the cumulative XUV flux (J/m**2)
    felt by a planet over its lifetime. It comes from Zahnle + Catling (2017) Equation 27,
    where they say they did integrals over the Lammer et al. (2009) XUV flux relations.
    This effectively assumes that the early times dominate, so the time integral doesn't
    depend (?!?) on the age of the system. It's very rough!

    Parameters
    ----------
    distribution : bool
        If False, return a simple array of values.
        If True, return an astropy.uncertainty.Distribution,
        which can be used for error propagation.
    """
    xuv_proxy = (self.stellar_luminosity(distribution=distribution) / u.Lsun) ** -0.6
    return self.relative_insolation(distribution=distribution) * xuv_proxy


def teq(self, distribution=False, albedo=0, f=1 / 4, **kw):
    """
    Planet Equilibrium Temperature (K)

    Calculate the equilbrium temperature of the planet. With default keywords,
    this assumes zero-albedo and uniform heat redistribution.

    Parameters
    ----------
    albedo : float
        The planet's Bond albedo, meaning the fraction of incoming
        stellar light reflected away, integrated over wavelength.
        The remaining light is absorbed and needs to be emitted
        as thermal (often infrared) radiation, thus contributing
        to the equilibrium temperature.
    f : float
        The ratio of the effective area collecting light from the star
        (which is the cross sectional area of the planet's shadow pi*R**2)
        to the effective area emitting thermal radiation away to space
        (which depends on how efficiently the planet distributes heat).
        Higher values of f correspond to a smaller fraction of the planet's
        surface doing the emitting, and thus a higher temperature.
        Two limiting values are:
            f = 1/4 = planet uniformly distributes heat, thus radiating
                      from a spherical surface with 4*pi*R**2 area
            f = 2/3 = looking at the dayside of planet that instantly
                      reradiates absorbed heat, with stronger weighting
                      toward the substellar point meaning the effective
                      emission area is only 3/2*pi*R**2 (which is smaller
                      than the spherical area of the dayside hemisphere,
                      which would give f=1/2)
    distribution : bool
        If False, return a simple array of values.
        If True, return an astropy.uncertainty.Distribution,
        which can be used for error propagation.
    """
    S = self.insolation(distribution=distribution)
    sigma = con.sigma_sb
    teq = ((S * f * (1 - albedo) / sigma) ** (1 / 4)).to(u.K)
    return teq


def planet_luminosity(self, distribution=False, **kw):
    """
    Planet Luminosity (K)

    Calculate total thermal luminosity of the planet,
    due to reradiation of absorbed starlight.

    Parameters
    ----------
    albedo : float
        The planet's Bond albedo, meaning the fraction of incoming
        stellar light reflected away, integrated over wavelength.
        The remaining light is absorbed and needs to be emitted
        as thermal (often infrared) radiation, thus contributing
        to the equilibrium temperature.
    distribution : bool
        If False, return a simple array of values.
        If True, return an astropy.uncertainty.Distribution,
        which can be used for error propagation.
    """
    S = self.insolation(distribution=distribution)
    R = self.radius(distribution=distribution)
    L_in = (S * np.pi * R**2).to(u.W)
    L_out = L_in
    return L_out


def transit_depth_from_radii(self, distribution=False, **kw):
    """
    Transit Depth ((Rp/Rs)**2, unitless)

    Calculate the transit depth from the planet and star radius.
    This calculates simply (Rp/Rs)**2; it neglects the effects of
    limb-darkening and/or grazing transits.

    Parameters
    ----------
    distribution : bool
        If False, return a simple array of values.
        If True, return an astropy.uncertainty.Distribution,
        which can be used for error propagation.
    """

    Rp = self.radius(distribution=distribution)
    Rs = self.stellar_radius(distribution=distribution)
    depth = (Rp / Rs).decompose() ** 2
    return depth


def transit_depth(self, distribution=False, **kw):
    """
    Transit Depth ((Rp/Rs)**2, unitless)

    Retrieve transit depth first from the standardized table,
    then calculated from planet radius + stellar radius.

    Parameters
    ----------
    distribution : bool
        If False, return a simple array of values.
        If True, return an astropy.uncertainty.Distribution,
        which can be used for error propagation.
    """
    return self._choose_calculation(
        methods=[
            "transit_depth_from_table",
            "transit_depth_from_radii",
        ],
        distribution=distribution,
        **kw,
    )


def transit_duration_from_orbit(self, distribution=False, **kw):
    """
    Transit Duration (days)

    Calculate the total transit duration (= the time between 1st and
    4th contact, from when first touches the stellar disk to when it
    last leaves). This calculation might be used if there is no
    transit duration available in the standardized table.

    Parameters
    ----------
    distribution : bool
        If False, return a simple array of values.
        If True, return an astropy.uncertainty.Distribution,
        which can be used for error propagation.
    """

    P = self.period(distribution=distribution)
    scaled_semimajoraxis = self.scaled_semimajoraxis(distribution=distribution)
    b = self.transit_impact_parameter(distribution=distribution)
    k = self.scaled_radius(distribution=distribution)
    T0 = P / np.pi / scaled_semimajoraxis
    T = T0 * np.sqrt((1 + k**2) - b**2)

    e = self.eccentricity(distribution=distribution)
    omega = self.argument_of_periastron(distribution=distribution)
    factor = np.sqrt(1 - e**2) / (1 + e * np.sin(omega))

    duration = (T * factor).to(u.day)

    return duration


def transit_duration(self, distribution=False, **kw):
    """
    Transit Duration (days)

    Retrieve transit duration first from the standardized table,
    then calculated from the orbital parameters

    Parameters
    ----------
    distribution : bool
        If False, return a simple array of values.
        If True, return an astropy.uncertainty.Distribution,
        which can be used for error propagation.
    """
    return self._choose_calculation(
        methods=[
            "transit_duration_from_table",
            "transit_duration_from_orbit",
        ],
        distribution=distribution,
        **kw,
    )


@property
def kludge_mass(self):
    """
    Have a safe way to calculate the mass of planets,
    that fills in gaps as necessary. Basic strategy:

        First from table.
        Then from msini.
    """

    # pull out the actual values from the table
    M = self.standard["mass"].copy()

    # try to replace bad ones with NVK3L
    bad = np.isfinite(M) == False
    self._speak(f"{sum(bad)}/{len(self)} masses are missing")

    # estimate from the msini
    try:
        M[bad] = self.msini[bad]
    except (KeyError, AssertionError, AtlasError, AttributeError):
        pass

    # replace those that are still bad with the a/R*
    stillbad = np.isfinite(M) == False
    self._speak(f"{sum(stillbad)}/{len(self)} are still missing after msini")

    return M


@property
def kludge_radius(self):
    """
    Have a safe way to calculate the radii of planets,
    that fills in gaps as necessary. Basic strategy:

        First from table.
        Then from mass, via Chen & Kipping (2017).
    """

    # pull out the actual values from the table
    R = self.standard["radius"].copy()

    # try to replace bad ones with NVK3L
    bad = np.isfinite(R) == False
    self._speak(f"{sum(bad)}/{len(self)} radii are missing")

    # estimate from Chen and Kipping
    try:
        M = self.kludge_mass
        R[bad] = estimate_radius(M[bad])
    except (KeyError, AssertionError, AtlasError, AttributeError):
        pass

    # replace those that are still bad with the a/R*
    stillbad = np.isfinite(R) == False
    self._speak(
        f"{sum(stillbad)}/{len(self)} are still missing after Chen & Kipping (2017)"
    )

    return R


@property
def kludge_age(self):
    """
    Have a safe way to calculate the age of planets,
    that fills in gaps as necessary. Basic strategy:

        First from table.
        Then assume 5 Gyr.
    """

    # pull out the actual values from the table
    age = self.standard["stellar_age"].copy()

    # try to replace bad ones with NVK3L
    bad = np.isfinite(age) == False
    self._speak(f"{sum(bad)}/{len(self)} ages are missing")

    # estimate from the msini
    try:
        age[bad] = 5 * u.Gyr
    except (KeyError, AssertionError, AtlasError, AttributeError):
        pass

    # replace those that are still bad with the a/R*
    stillbad = np.isfinite(age) == False
    self._speak(
        f"{sum(stillbad)}/{len(self)} are still missing after blindly assuming 5Gyr for missing ages"
    )

    return age


@property
def surface_gravity(self):
    """
    (FIXME) -- make an assumption for planets without masses
    """

    G = con.G
    M = self.mass
    R = self.radius

    g = (G * M / R**2).to("m/s**2")
    return g


@property
def density(self):
    """
    The density of the planet.
    """
    mass = self.mass
    volume = 4 / 3 * np.pi * (self.radius) ** 3
    return (mass / volume).to("g/cm**3")


@property
def escape_velocity(self):
    """
    The escape velocity of the planet.
    """
    G = con.G
    M = self.mass
    R = self.radius
    return np.sqrt(2 * G * M / R).to("km/s")


@property
def orbital_velocity(self):
    return (2 * np.pi * self.semimajoraxis / self.period).to("km/s")


@property
def impact_velocity(self):
    return np.sqrt(self.orbital_velocity**2 + self.escape_velocity**2)


@property
def escape_parameter(self):
    """
    The Jeans atmospheric escape parameter for atomic hydrogen,
    at the equilibrium temperature of the planet.
    """
    k = con.k_B
    T = self.teq
    mu = 1
    m_p = con.m_p
    G = con.G
    M = self.mass
    R = self.radius

    e_thermal = k * T
    e_grav = G * M * m_p / R
    return (e_grav / e_thermal).decompose()


def scale_height(self, mu=2.32):
    """
    The scale height of the atmosphere, at equilibrium temperature.
    """
    k = con.k_B
    T = self.teq
    m_p = con.m_p
    g = self.surface_gravity
    return (k * T / mu / m_p / g).to("km")
