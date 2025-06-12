import astropy.constants as cc
import astropy.units as u
import numpy as np

def calc_error(telescope_params, 
               photons_raw, 
               tE=30*u.s, 
               dλ=550*u.nm, 
               npix=0, 
               atm_transmission=1, 
               airmass=1, 
               h0=8000*u.m, 
               nE=1, 
               N_sky=0 / u.s / u.m**2, 
               debug=False):

    """
    Calculates the expected error for a given set
    of telescope / stellar properties.

    Parameters:
    -----------
    telescope_params :: dict
        Contains all the import operations properties
        of the telescope you are observing with. Should
        include...
            - diameter
            - collecting_area
            - QE (quantum efficiency)
            - gain
            - RN (read noise)
            - R_dark (dark rate)
            - altitude
    photons_raw :: int/float/array
        Number of photons recieved at Earth
    tE :: astropy.Quantity
        Exposure time   
    dλ :: astropy.Quantity
        Wavelength bin length
    npix :: int
        Number of pixels in extraction aperture
    atm_transmission :: float
        Amount of light that passes through the atmosphere
        (must be between 0 and 1)
    airmass :: float
        Airmass of observation
    h0 :: astropy.Quantity
        Scale height of the atmosphere
    nE :: int
        Number of uncorrelated reference stars
    debug :: bool
        Allows for diagnostic values to be plotted
    """

    ### Loads telescope-specific data
    D = telescope_params["diameter"]
    A = telescope_params["collecting_area"]
    QE = telescope_params["QE"]
    gain = telescope_params["gain"]
    RN = telescope_params["RN"]
    R_dark = telescope_params["R_dark"]
    h = telescope_params["altitude"]

    ### Provides assumed units if none are provided
    if not isinstance(D, u.Quantity): D *= u.m;
    if not isinstance(A, u.Quantity): A *= u.m**2;
    if not isinstance(R_dark, u.Quantity): R_dark *= 1/u.s;
    if not isinstance(h, u.Quantity): h *= u.m;
    
    ### Ensures that single-value calls can run / handles annoying photon units
    if photons_raw.isscalar: photons_raw = np.array([photons_raw.value]) * photons_raw.unit;
    if photons_raw.unit.is_equivalent(u.ph / (u.s * u.m**2 * u.nm)): photons_raw /= u.ph;    

    ### Calculates the number of photons recieved at CCD
    N_star = photons_raw*A*tE*dλ*QE*atm_transmission

    ### Calculates components of error
    σN = np.sqrt(N_star + N_sky * tE * A)/N_star
    σDark = np.sqrt(R_dark * tE * npix)/N_star
    σScint = ((0.135 * u.cm**(2/3) * u.s**(1/2)) * D**(-2/3) * airmass**1.75 / np.sqrt(2*tE) * np.exp(-h/h0) * np.sqrt(1 + 1/nE)).decompose()
    σRN = np.sqrt(npix * RN**2)/N_star

    ### Calculates the total error
    σtotal = np.sqrt(σN**2 + σDark**2 + σScint**2 + σRN**2)

    if debug:

        print(
f"""
 ------------------------------
|    Error Source Breakdown    |
 ------------------------------

     Photon Noise (σN):
         min: {np.round(min(σN), 6)}
         max: {np.round(max(σN), 6)}

     Dark Noise (σDark):
         min: {np.round(min(σDark), 6)}
         max: {np.round(max(σDark), 6)}

     Scint Noise (σScint):
            ~ {np.round(σScint, 6)}

     Read Noise (σRN):
         min: {np.round(min(σRN), 6)}
         max: {np.round(max(σRN), 6)}

- - - - - - - - - - - - - - - -

     Total Noise (σtotal):
         min: {np.round(min(σtotal), 6)}
         max: {np.round(max(σtotal), 6)}

-------------------------------
""")
        
    return σtotal

def calc_AB_photon_flux(magnitude, ref_wavelength=639.07*u.nm):

    """
    Calculates the photons/(s m2 nm) corresponding
    to a given AB magnitude

    Parameters:
    -----------
    magnitude :: float
        AB magnitude of a star
    ref_wavelength :: astropy.Quantity
        Reference wavelength of bandpass

    Returns:
    --------
    AB_phot_counts :: astropy.Quantity
        Number of photons/(s m2 nm)
    """

    ### Ensures that all variables of some astropy unit
    if not isinstance(ref_wavelength, u.Quantity): ref_wavelength *= u.nm;

    ### Calculates the energy of a single reference photon
    E_phot = (cc.h*cc.c/u.photon)/ref_wavelength

    ### Calculates the photon flux in the AB system
    AB_fluxes = (1 * u.W/u.m**2/u.Hz) * 10**(-(magnitude+56.1)/2.5)
    AB_phot_counts = (AB_fluxes * cc.c/ref_wavelength**2/E_phot).to(u.photon/(u.s * u.m**2 * u.nm))
    
    return AB_phot_counts

def calc_AB_mag(photon_flux, ref_wavelength=639.07*u.nm):

    """
    Calculates AB magnitude corresponding to
    a given photon flux.

    Parameters:
    -----------
    photon_flux :: astropy.Quantity
        Photon flux in photons/(s m2 nm)
    ref_wavelength :: astropy.Quantity
        Reference wavelength of bandpass

    Returns:
    --------
    AB_mag :: astropy.Quantity
        AB magnitude corresponding to the given flux
    """

    ### Ensures that all variables of some astropy unit
    if not isinstance(photon_flux, u.Quantity): photon_flux *= u.photon/(u.s * u.m**2 * u.nm);
    if not isinstance(ref_wavelength, u.Quantity): ref_wavelength *= u.nm;

    ### Calculates the energy of a single reference photon
    E_phot = (cc.h*cc.c/u.photon)/ref_wavelength

    ### Calculates the photon flux in the AB system
    AB_fluxes = photon_flux * (ref_wavelength**2 * E_phot / cc.c)
    AB_mag = -2.5*np.log10((AB_fluxes / (1 * u.W/u.m**2/u.Hz)).decompose()) - 56
    
    return AB_mag