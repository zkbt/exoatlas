from ..imports import *


def define_telescope_unit(
    telescope_name="telescope",
    area=1 * u.m**2,
    wavelength=5 * u.micron,
    R=20,
    dt=1 * u.hr,
    **kw,
):
    """
    Create a custom astropy unit to represent
    the collecting area of a generic telscope
    with a particular collecting area
    at a particular wavelength
    at a particular resolution
    observing for a particular time.

    Parameters
    ----------
    telescope_name : str
        The name of the telescope.

    area : astropy.units.quantity.Quantity
        The collecting area of the telescope.

    wavelength : astropy.unit.Quantity
        The wavelength at which it should be calculated.

    R : float
        The spectral resolution at which the
        telescope will bin wavelengths..
        (Ignored if telescope is None.)

    dt : astropy.units.quantity.Quantity
        The time over which the telescope exposes.
        (Ignored if telescope is None.)
    """

    dw = wavelength / R
    unit = u.def_unit(
        f"[{telescope_name} | {wavelength} | {dt} | R={R}]",
        dt * area * dw,
        doc=f"""
                      This custom unit represents the
                      photon-collecting power of the
                      {telescope_name} telescope, when integrating
                      for a time of {dt}, at a spectral
                      resolution of R={R} (evaluated
                      at {wavelength}).
                      """,
    )
    unit.telescope_name = telescope_name
    unit.wavelength = wavelength
    unit.dt = dt
    unit.R = R
    return unit


def define_JWST_unit(wavelength=5 * u.micron, efficiency=0.5, **kw):
    """
    Create a JWST telescope unit.

    Parameters
    ----------
    wavelength : astropy.unit.Quantity
        The wavelength at which it should be calculated.

    R : float
        The spectral resolution at which the
        telescope will bin wavelengths..
        (Ignored if telescope is None.)

    dt : astropy.units.quantity.Quantity
        The time over which the telescope exposes.
        (Ignored if telescope is None.)
    """
    return define_telescope_unit(
        telescope_name="JWST",
        wavelength=wavelength,
        area=efficiency * 25 * u.m**2,
        **kw,
    )


def define_HST_unit(wavelength=1.4 * u.micron, efficiency=0.3, **kw):
    """
    Create a HST telescope unit.

    Parameters
    ----------
    wavelength : astropy.unit.Quantity
        The wavelength at which it should be calculated.

    R : float
        The spectral resolution at which the
        telescope will bin wavelengths..
        (Ignored if telescope is None.)

    dt : astropy.units.quantity.Quantity
        The time over which the telescope exposes.
        (Ignored if telescope is None.)
    """
    return define_telescope_unit(
        telescope_name="HST",
        wavelength=wavelength,
        area=efficiency * 4.5 * u.m**2,
        **kw,
    )


def define_APO_unit(wavelength=0.5 * u.micron, efficiency=0.1, **kw):
    """
    Create an APO telescope unit.

    Parameters
    ----------
    wavelength : astropy.unit.Quantity
        The wavelength at which it should be calculated.

    R : float
        The spectral resolution at which the
        telescope will bin wavelengths..
        (Ignored if telescope is None.)

    dt : astropy.units.quantity.Quantity
        The time over which the telescope exposes.
        (Ignored if telescope is None.)
    """
    radius = 3.5 * u.m / 2
    area = np.pi * radius**2
    return define_telescope_unit(
        telescope_name="APO", wavelength=wavelength, area=area * efficiency, **kw
    )


def define_SBO_unit(wavelength=0.5 * u.micron, efficiency=0.1, **kw):
    """
    Create an APO telescope unit.

    Parameters
    ----------
    wavelength : astropy.unit.Quantity
        The wavelength at which it should be calculated.

    R : float
        The spectral resolution at which the
        telescope will bin wavelengths..
        (Ignored if telescope is None.)

    dt : astropy.units.quantity.Quantity
        The time over which the telescope exposes.
        (Ignored if telescope is None.)
    """
    radius = 0.5 * u.m / 2
    area = np.pi * radius**2
    return define_telescope_unit(
        telescope_name="SBO", wavelength=wavelength, area=area * efficiency, **kw
    )


def define_TESS_unit(wavelength=0.8 * u.micron, R=2, **kw):
    """
    Create a TESS telescope unit.

    Parameters
    ----------
    wavelength : astropy.unit.Quantity
        The wavelength at which it should be calculated.

    R : float
        The spectral resolution at which the
        telescope will bin wavelengths..
        (Ignored if telescope is None.)

    dt : astropy.units.quantity.Quantity
        The time over which the telescope exposes.
        (Ignored if telescope is None.)
    """
    return define_telescope_unit(
        telescope_name="TESS", wavelength=wavelength, area=0.0086 * u.m**2, R=R, **kw
    )


def define_Kepler_unit(wavelength=0.65 * u.micron, R=2, **kw):
    """
    Create a Kepler telescope unit.

    Parameters
    ----------
    wavelength : astropy.unit.Quantity
        The wavelength at which it should be calculated.

    R : float
        The spectral resolution at which the
        telescope will bin wavelengths..
        (Ignored if telescope is None.)

    dt : astropy.units.quantity.Quantity
        The time over which the telescope exposes.
        (Ignored if telescope is None.)
    """
    return define_telescope_unit(
        telescope_name="Kepler", wavelength=wavelength, area=0.708 * u.m**2, R=R, **kw
    )


# create a dictionary of unit factories
telescope_units = dict(
    Kepler=define_Kepler_unit,
    TESS=define_TESS_unit,
    HST=define_HST_unit,
    JWST=define_JWST_unit,
    APO=define_APO_unit,
    SBO=define_SBO_unit,
)


def define_telescope_unit_by_name(telescope_name, wavelength=None, **kw):
    """
    Wrapper to create a telescope unit
    by passing the telescope name.

    Parameters
    ----------
    telescope_name : str
        The name of the telescope.

    wavelength : astropy.unit.Quantity
        The wavelength at which it should be calculated.

    R : float
        The spectral resolution at which the
        telescope will bin wavelengths..
        (Ignored if telescope is None.)

    dt : astropy.units.quantity.Quantity
        The time over which the telescope exposes.
        (Ignored if telescope is None.)
    """
    if wavelength is None:
        return telescope_units[telescope_name](**kw)
    else:
        return telescope_units[telescope_name](wavelength=wavelength, **kw)


# define a useful big number of photons (=30ppm)
lotsofphotons_unit = u.def_unit("Gigaphotons", 1e9 * u.ph)
