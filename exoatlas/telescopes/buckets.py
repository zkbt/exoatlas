from ..imports import *

def define_telescope_unit(name='telescope',
                          area=1*u.m**2,
                          wavelength=5*u.micron,
                          R=20,
                          dt=1*u.hr):
    '''
    Create a custom astropy unit to represent
    the collecting area of JWST
    at a R=20 spectral resolution
    integrating over one hour.

    Parameters
    ----------
    name : str
        The name of the telescope.

    area : astropy.units.quantity.Quantity
        The collecting area of the telescope.

    wavelength : astropy.unit.Quantity
        The wavelength at which it should be calculated.

    R : float
        The spectral resolution at which the
        telescope will be binned.
        (Ignored if telescope is None.)

    dt : astropy.units.quantity.Quantity
        The time over which the telescope exposes.
        (Ignored if telescope is None.)
    '''


    dw = wavelength/R
    unit = u.def_unit(f'[{name} | {dt} | R={R}]', #@{wavelength}
                      dt*area*dw,
                      doc=f'''
                      This custom unit represents the
                      photon-collecting power of the
                      {name} telescope, when integrating
                      for a time of {dt}, at a spectral
                      resolution of R={R} (evaluated
                      at {wavelength}).
                      ''')
    unit.wavelength = wavelength
    return unit

def define_JWST_unit(wavelength=5*u.micron, **kw):
    return define_telescope_unit(name='JWST',
                                 wavelength=wavelength,
                                 area=25*u.m**2,
                                 **kw)

def define_HST_unit(wavelength=1.4*u.micron, **kw):
    return define_telescope_unit(name='HST',
                                 wavelength=wavelength,
                                 area=4.5*u.m**2,
                                 **kw)

def define_TESS_unit(wavelength=0.8*u.micron, R=2, **kw):
    return define_telescope_unit(name='TESS',
                                 wavelength=wavelength,
                                 area=0.0086*u.m**2,
                                 R=R,
                                 **kw)

def define_Kepler_unit(wavelength=0.65*u.micron, R=2, **kw):
    return define_telescope_unit(name='Kepler',
                                 wavelength=wavelength,
                                 area=0.708*u.m**2,
                                 R=R,
                                 **kw)

telescope_units = dict(Kepler=define_Kepler_unit,
                       TESS=define_TESS_unit,
                       HST=define_HST_unit,
                       JWST=define_JWST_unit)

def define_telescope_unit_by_name(name, wavelength=None, **kw):
    if wavelength is None:
        return telescope_units[name](**kw)
    else:
        return telescope_units[name](wavelength=wavelength, **kw)

photon_unit = u.def_unit('Gigaphotons', 1e9*u.ph)
