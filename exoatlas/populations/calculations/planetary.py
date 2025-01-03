from ...imports import *
from ...models.chen import *


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
    xuv_proxy = (
        self.stellar_luminosity(distribution=distribution).to_value(u.Lsun) ** -0.6
    )
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


def scaled_radius_from_radii(self, distribution=False, **kw):
    """
    Scaled Planet Radius ((Rp/Rs), unitless)

    Calculate the radius ratio from the planet and star radius.
    This simply calculates (Rp/Rs); it might be used if no
    radius ratio is provided in the table.

    Parameters
    ----------
    distribution : bool
        If False, return a simple array of values.
        If True, return an astropy.uncertainty.Distribution,
        which can be used for error propagation.
    """

    Rp = self.radius(distribution=distribution)
    Rs = self.stellar_radius(distribution=distribution)
    ratio = (Rp / Rs).decompose()
    return ratio


def scaled_radius(self, distribution=False, **kw):
    """
    Scaled Planet Radius ((Rp/Rs), unitless)

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
            "transit_scaled_radius_from_table",
            "scaled_radius_from_radii",
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


def mass_estimated_from_radius(self, distribution=False, **kw):
    """
    Estimated Planet Mass (Earth masses)

    Retrieve an estimate of the planet mass from the
    planet's radius, using a very approximate
    empirical relation.

    Parameters
    ----------
    distribution : bool
        If False, return a simple array of values.
        If True, return an astropy.uncertainty.Distribution,
        which can be used for error propagation.

    """
    if distribution:
        raise RuntimeError(
            "mass_estimated_from_radius does not propagate uncertainties yet; sorry!"
        )
    radius = self.get("radius", distribution=False)
    return use_chen_and_kipping_to_estimate_mass_from_radius(R=radius)


def radius_estimated_from_mass(self, distribution=False, **kw):
    """
    Estimated Planet Radius (Earth radii)

    Retrieve an estimate of the planet radius from the
    planet's mass, using a very approximate
    empirical relation.

    Parameters
    ----------
    distribution : bool
        If False, return a simple array of values.
        If True, return an astropy.uncertainty.Distribution,
        which can be used for error propagation.

    """
    if distribution:
        raise RuntimeError(
            "radius_estimated_from_mass does not propagate uncertainties yet; sorry!"
        )

    mass = self.get("mass", distribution=False)
    return use_chen_and_kipping_to_estimate_mass_from_radius(M=mass)


def kludge_mass(self, distribution=False, **kw):
    """
    Planet Mass or msini (Earth masses)

    Retrieve an estimate of the planet mass,
    starting first from an actual published table mass,
    and then from msini assuming sini=1,
    and then from an emprical radius-mass relation.

    Parameters
    ----------
    distribution : bool
        If False, return a simple array of values.
        If True, return an astropy.uncertainty.Distribution,
        which can be used for error propagation.

    """
    return self._choose_calculation(
        methods=["mass_from_table", "msini_from_orbit", "mass_estimated_from_radius"],
        distribution=distribution,
        **kw,
    )


def kludge_radius(self, distribution=False, **kw):
    """
    Planet Radius or Estimated Planet Radius (Earth radii)

    Retrieve an estimate of the planet radius,
    starting first from an actual published table radius,
    and then from a mass-radius relation.

    Parameters
    ----------
    distribution : bool
        If False, return a simple array of values.
        If True, return an astropy.uncertainty.Distribution,
        which can be used for error propagation.

    """
    return self._choose_calculation(
        methods=[
            "radius_from_table",
            "radius_estimated_from_mass",
        ],
        distribution=distribution,
        **kw,
    )


def kludge_stellar_age(self, distribution=False, **kw):
    """
    System Age (Gyr)

    Retrieve an estimate of the system age,
    starting first from an actual published table age,
    and then boldy/foolishly assuming it's 5 Gyr.

    Parameters
    ----------
    distribution : bool
        If False, return a simple array of values.
        If True, return an astropy.uncertainty.Distribution,
        which can be used for error propagation.

    """

    # pull out the actual values from the table
    age = self.get("stellar_age", distribution=distribution)

    # try to replace bad ones with NVK3L
    bad = np.isfinite(age) == False
    age[bad] = 5 * u.Gyr

    return age


def surface_gravity(self, kludge=False, distribution=False, **kw):
    """
    Planet Surface Gravity (m/s**2)

    Calculate the planet's surface gravity
    from its mass and radius.

    TO-DO:
    Explore if we should choose a different calculation more
    closely tied to observables to minimize uncertainties.

    Parameters
    ----------
    kludge : bool
        Should we include kludged estimates for mass (from msini and/or
        empirical mass-radius) and/or radius (from empircal mass-radius)
        when doing this calculation?
    distribution : bool
        If False, return a simple array of values.
        If True, return an astropy.uncertainty.Distribution,
        which can be used for error propagation.
    """

    # choose whether or not include estimated masses/radii
    if kludge:
        M = self.get("kludge_mass", distribution=distribution)
        R = self.get("kludge_radius", distribution=distribution)
    else:
        M = self.get("mass", distribution=distribution)
        R = self.get("radius", distribution=distribution)

    G = con.G
    g = (G * M / R**2).to("m/s**2")
    return g


def density(self, kludge=False, distribution=False, **kw):
    """
    Planet Density (m/s**2)

    Calculate the planet's bulk density
    from its mass and radius.

    TO-DO:
    Explore if we should choose a different calculation more
    closely tied to observables to minimize uncertainties.

    Parameters
    ----------
    kludge : bool
        Should we include kludged estimates for mass (from msini and/or
        empirical mass-radius) and/or radius (from empircal mass-radius)
        when doing this calculation?
    distribution : bool
        If False, return a simple array of values.
        If True, return an astropy.uncertainty.Distribution,
        which can be used for error propagation.
    """

    # choose whether or not include estimated masses/radii
    if kludge:
        M = self.get("kludge_mass", distribution=distribution)
        R = self.get("kludge_radius", distribution=distribution)
    else:
        M = self.get("mass", distribution=distribution)
        R = self.get("radius", distribution=distribution)

    V = 4 / 3 * np.pi * R**3
    density = (M / V).to("g/cm**3")
    return density


def escape_velocity(self, kludge=False, distribution=False, **kw):
    """
    Planet Escape Velocity (km/s)

    Calculate the planet's escape velocity
    from its mass and radius.

    TO-DO:
    Explore if we should choose a different calculation more
    closely tied to observables to minimize uncertainties.

    Parameters
    ----------
    kludge : bool
        Should we include kludged estimates for mass (from msini and/or
        empirical mass-radius) and/or radius (from empircal mass-radius)
        when doing this calculation?
    distribution : bool
        If False, return a simple array of values.
        If True, return an astropy.uncertainty.Distribution,
        which can be used for error propagation.
    """

    if kludge:
        M = self.get("kludge_mass", distribution=distribution)
        R = self.get("kludge_radius", distribution=distribution)
    else:
        M = self.get("mass", distribution=distribution)
        R = self.get("radius", distribution=distribution)

    G = con.G
    escape_velocity = np.sqrt(2 * G * M / R).to("km/s")
    return escape_velocity


def orbital_velocity(self, distribution=False, **kw):
    """
    Planet Orbital Velocity (km/s)

    Calculate the average tangential orbital speed of the
    the planet in its orbit around its star.

    Parameters
    ----------
    distribution : bool
        If False, return a simple array of values.
        If True, return an astropy.uncertainty.Distribution,
        which can be used for error propagation.
    """
    a = self.semimajoraxis(distribution=distribution)
    P = self.period(distribution=distribution)
    return (2 * np.pi * a / P).to("km/s")


def impact_velocity(self, distribution=False, **kw):
    """
    Very Approximately Estimated Impact Velocity (km/s)

    Calcuate a back-of-the-envelope estimate of the typical
    impact velocity which which something might hit the planet,
    as a combination of orbital and escape velocity.

    Parameters
    ----------
    distribution : bool
        If False, return a simple array of values.
        If True, return an astropy.uncertainty.Distribution,
        which can be used for error propagation.
    """
    v_orbital = self.orbital_velocity(distribution=distribution)
    v_escape = self.escape_velocity(distribution=distribution)
    return np.sqrt(v_orbital**2 + v_escape**2)


def escape_parameter(
    self,
    temperature="teq",
    mu=1,
    albedo=0,
    f=1 / 4,
    kludge=False,
    distribution=False,
    **kw
):
    """
    Jeans Atmospheric Escape Parameter (unitless)

    Calculate the Jeans escape parameter, as the ratio
    of a particle's gravitational binding energy to its
    thermal energy. This is a super-approximate tracer
    for susceptability to atmospheric loss, depending on
    lots of quantities that we often in truth don't know.

    TO-DO:
    Think about whether there is an interesting intermediate
    to assume for exospheric temperature, or whether we should
    just keep this as a simple relative scaling.

    TO-DO:
    Should we include a correction factor for different
    gravitational acceleration at high altitudes?

    Parameters
    ----------
    temperature : str, Quantity
        The temperature to use for the planet's exosphere,
        in units of Kelvin. The default string `temperature='teq'`
        will calculate an equilibrium temperature (assuming some
        albedo and heat redistribution parameter); otherwise
        something like `temperature=1000*u.K` will assume a
        constant exospheric temperature. Neither approximation
        is probably super great.
    mu : float
        Mean molecular weight of escaping particle,
        in atomic mass units. The default of mu=1
        corresponds to atomic hydrogen.
    albedo : float
        Albedo for calculating equilibrium temperature
        (see docstring for `.teq()` for details)
    f : float
        Heat redistribution for calculating equilibrium temperature
        (see docstring for `.teq()` for details)
    kludge : bool
        Should we include kludged estimates for mass (from msini and/or
        empirical mass-radius) and/or radius (from empircal mass-radius)
        when doing this calculation?
    distribution : bool
        If False, return a simple array of values.
        If True, return an astropy.uncertainty.Distribution,
        which can be used for error propagation.
    """
    if temperature == "teq":
        T = self.teq(albedo=albedo, f=f, distributio=distribution)
    else:
        T = temperature

    if kludge:
        M = self.get("kludge_mass", distribution=distribution)
        R = self.get("kludge_radius", distribution=distribution)
    else:
        M = self.get("mass", distribution=distribution)
        R = self.get("radius", distribution=distribution)

    k = con.k_B
    m_p = con.m_p
    G = con.G
    e_thermal = k * T
    e_grav = G * M * m_p / R
    return (e_grav / e_thermal).decompose()


def scale_height(self, mu=2.3, albedo=0, f=1 / 4, kludge=False, distribution=False):
    """
    Atmospheric Scale Height (km)

    Calculate the scale height of the planet's atmosphere,
    assuming a mean molecular weight and the

    TO-DO:
    Consider alternate ways of estimating mean-molecular weight
    as a function of planet properties, so that (for example),
    we could use this to calculate scale heights for primary
    and secondary atmospheres together. There are probably
    so many unknowns that would need to go into this that
    we should make people makes those choices themselves.

    Parameters
    ----------
    mu : float
        Mean molecular weight of the atmosphere,
        in atomic mass units. The default of mu=2.3
        corresponds approximately to chemical equilibrium
        of a solar composition atmosphere at fairly
    albedo : float
        Albedo for calculating equilibrium temperature
        (see docstring for `.teq()` for details)
    f : float
        Heat redistribution for calculating equilibrium temperature
        (see docstring for `.teq()` for details)
    kludge : bool
        Should we include kludged estimates for mass (from msini and/or
        empirical mass-radius) and/or radius (from empircal mass-radius)
        when doing this calculation?
    distribution : bool
        If False, return a simple array of values.
        If True, return an astropy.uncertainty.Distribution,
        which can be used for error propagation.
    """
    k = con.k_B
    T = self.teq(albedo=albedo, f=f, distribution=distribution)
    m_p = con.m_p
    g = self.surface_gravity(kludge=kludge, distribution=distribution)
    H = (k * T / mu / m_p / g).to("km")
    return H
