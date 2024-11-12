from ...imports import *
from .plottable import *
from ...telescopes import *


class Flux(PlottableAxis):
    source = "relative_insolation"
    label = "Bolometric Flux Received (relative to Earth)"
    scale = "log"
    lim = [6e4, 2e-4]


class CumulativeXUVFlux(Flux):
    source = "relative_cumulative_xuv_insolation"
    label = "Time-Integrated XUV Flux Received (relative to Earth)"


class ImpactVelocity(PlottableAxis):
    source = "impact_velocity"
    label = "Estimated Impact Velocity (km/s)"
    scale = "log"


class LogFlux(PlottableAxis):
    source = "relative_insolation"
    label = "Bolometric Flux Received (relative to Earth)"
    scale = "log"
    lim = [6e4, 2e-4]


class Radius(PlottableAxis):
    source = "radius"
    label = "Planet Radius (Earth radii)"
    scale = "log"
    lim = [0.3, 30]

    def value_lowerupper(self):
        return self.panel.pop.get_uncertainty_lowerupper("radius")


class Mass(PlottableAxis):
    source = "mass"
    label = "Planet Mass\n(Earth masses)"
    scale = "log"
    lim = [0.03, 3000]

    def value_lowerupper(self):
        return self.panel.pop.get_uncertainty_lowerupper("mass")


class SemimajorAxis(PlottableAxis):
    source = "semimajoraxis"
    label = "Semimajor Axis (AU)\n"
    scale = "log"
    lim = [0.001, 1000]


class AngularSeparation(PlottableAxis):
    source = "angular_separation"
    label = "Angular Separation (arcsec)\n"
    scale = "log"
    lim = [0.001, 10]


class Contrast(PlottableAxis):
    source = "imaging_contrast"
    label = "Planet-to-Star Contrast"
    scale = "log"
    lim = [1e-10, 1e-3]

    def __init__(self, phase_function=0.25, albedo=0.25, **kw):
        """
        Initialize for a particular wavelength, because the
        eclipse depth will depend on the thermal emission
        spectrum of the planet and star.
        """
        PlottableAxis.__init__(self, **kw)
        self.phase_function = phase_function
        self.albedo = albedo
        self.label = f"Reflected Light Planet-to-Star Contrast\n(albedo = {albedo:.0%}, phase function = {phase_function:.0%})"

    def value(self):
        return self.panel.pop.imaging_contrast * self.albedo * self.phase_function


class KludgedMass(PlottableAxis):
    source = "kludge_mass"
    label = "Planet Mass or msini\n(Earth masses)"
    scale = "log"
    lim = [0.03, 4131]


class StellarTeff(PlottableAxis):
    source = "stellar_teff"
    label = "Stellar Temperature (K)"
    scale = "linear"
    lim = [2000, 7000]


class StellarLuminosity(PlottableAxis):
    source = "stellar_luminosity"
    label = "Stellar Luminosity (L_$\odot$)"
    scale = "log"
    lim = [None, None]


class Distance(PlottableAxis):
    source = "distance"
    label = "Distance\n(parsecs)"
    scale = "log"
    lim = [5, 1000]


class EscapeVelocity(PlottableAxis):
    source = "escape_velocity"
    label = "Escape\nVelocity\n(km/s)"
    scale = "log"
    lim = [2, 500]


class Escape(PlottableAxis):
    source = "escape_parameter"
    label = "$\lambda = E_{grav}/E_{thermal}$"
    scale = "log"
    lim = [None, None]

    def value_lowerupper(self):
        """
        Calculate a rough estimate of the uncertainty on the escape
        velocity, simply by propagating the mass and radius uncertainties.
        This could probably be done more cleverly by deriving it from
        the planetary surface gravity uncertainty, which as fewer hidden
        covariances because it can be derived directly from transit
        and RV observables.
        """

        m = self.panel.pop.mass
        sigma_m = self.panel.pop.get_uncertainty("mass")

        r = self.panel.pop.radius
        sigma_r = self.panel.pop.get_uncertainty("radius")

        dlnm = sigma_m / m
        dlnr = sigma_r / r

        dlne = np.sqrt((dlnm / 2) ** 2 + (dlnr / 2) ** 2)
        sigma = dlne * self.panel.pop.escape_velocity

        return sigma, sigma


class Density(PlottableAxis):
    source = "density"
    label = "Planet Density\n(g/cm$^3$)"
    scale = "log"
    lim = [0.01, 100]


class StellarRadius(PlottableAxis):
    source = "stellar_radius"
    label = "Stellar Radius\n(solar radii)"
    scale = "linear"
    lim = [0.0, 2.0]


class Period(PlottableAxis):
    source = "period"
    label = "Period (days)\n"
    scale = "log"
    lim = [0.15, 365]


class Jmag(PlottableAxis):
    source = "Jmag"
    label = "J (magnitude)\n"
    scale = "linear"
    lim = [3.5, 14.5]


class Depth(PlottableAxis):
    source = "transit_depth"
    label = "Transit Depth\nof Planet"
    scale = "log"
    lim = [2e-6, 2e-1]


class StellarBrightness(PlottableAxis):
    # source='stellar_brightness'
    scale = "log"
    lim = [None, None]  # [1e2, 1e8]

    def __init__(self, wavelength=1 * u.micron, **kw):
        """
        Initialize for a particular wavelength, because the
        eclipse depth will depend on the thermal emission
        spectrum of the planet and star.
        """
        PlottableAxis.__init__(self, **kw)
        self.wavelength = wavelength

        # set up the units
        self.unit = u.Unit(("ph s^-1 m^-2 micron^-1"))
        self.unit_string = "photons/s/m$^2$/$\mu$m"

        # set the label
        self.label = f'Stellar Brightness at Earth at $\lambda={self.wavelength.to("micron").value:.1f}\mu$m\n({self.unit_string})'

    def value(self):
        """
        What data value to plot?
        """
        return self.panel.pop.stellar_brightness(self.wavelength).to(self.unit)


class StellarBrightnessTelescope(PlottableAxis):
    scale = "log"
    lim = [None, None]  # [1e-3, 1e3]

    def __init__(
        self, panel=None, orientation=None, telescope_name="JWST", wavelength=None, **kw
    ):
        """
        Initialize the StellarBrightness plottable.
        It depends on wavelength, due to the thermal
        emission Planck spectrum.

        Parameters
        ----------
        telescope_name : None, str
            The telescope unit in which to express the
            stellar brightness. Options include:
            'JWST', 'Hubble', 'Kepler', 'TESS'

        wavelength : astropy.units.quantity.Quantity
            The wavelength at which we want the
            stellar brightness to be computed.

        R : float
            The spectral resolution at which the
            telescope will bin wavelengths..
            (Ignored if telescope is None.)

        dt : astropy.units.quantity.Quantity
            The time over which the telescope exposes.
            (Ignored if telescope is None.)
        """

        # initialize the basic plottable axis
        PlottableAxis.__init__(self, panel=panel, orientation=orientation, **kw)

        # do all the basic telescope + wavelength setup
        self.setup_telescope(telescope_name=telescope_name, wavelength=wavelength, **kw)

        # specify the unit to be used for plotting
        self.setup_unit()

        # reset the label, because it depends on the telescope inputs
        self.setup_label()

    def setup_telescope(self, telescope_name=None, wavelength=None, **kw):
        """
        Setup the basics of the telescope, wavelength, and units.
        """

        # break if there's not telescope
        assert telescope_name is not None

        # keep track of the telescope (if any)
        self.telescope_name = telescope_name

        # define a unit of collecting power for the telescope
        self.telescope_unit = define_telescope_unit_by_name(
            self.telescope_name, wavelength=wavelength, **kw
        )

        # note the wavelength the telescope is set at
        self.wavelength = self.telescope_unit.wavelength
        self.R = self.telescope_unit.R

    def setup_unit(self):
        """
        Define the unit for plotting.
        """
        # what's the astropy unit associated with the axis?
        self.unit = lotsofphotons_unit / self.telescope_unit

        # what should go in the () in the axis label?
        self.unit_string = self.unit.to_string()

    def setup_label(self):
        """
        How should this plottable be labled on an axis?
        """
        # define the label, based on the wavelength and telescope
        w = self.wavelength.to(u.micron).value
        self.label = (
            rf"Stellar Brightness at Earth at $\lambda={w}\mu$m\n({self.unit_string})"
        )

    def value(self):
        """
        What data value to plot?
        """
        return self.panel.pop.stellar_brightness(self.wavelength).to(self.unit)


class DepthSNR(StellarBrightnessTelescope):
    scale = "log"
    size_normalization = 10

    def __init__(self, **kw):
        """
        Initialize S/N for the basic transit depth.

        Parameters
        ----------
        wavelength : astropy.units.quantity.Quantity
            The wavelength at which we want the
            eclipse depth to be computed.
        """
        StellarBrightnessTelescope.__init__(self, **kw)
        self.kw = kw

    def setup_label(self):
        """
        How should this plottable be labled on an axis?
        """
        # define the label, based on the wavelength and telescope
        w = self.wavelength.to(u.micron).value
        self.label = rf"S/N for Transit Depth\nat $\lambda={self.wavelength.to(u.micron).value}\mu m$\n(R={self.R})"

    def value(self):
        """
        What data value to plot?
        """
        return self.panel.pop.depth_snr(**self.kw)


class Transmission(Depth):
    def __init__(self, mu=2.32, threshold=2, **kw):
        """
        Initialize for a particular mean molecular weight.

        Parameters
        ----------
        mu : float
            Mean molecular weight (default 2.2 for H/He)
        threshold : float
            By how many sigma must the planet mass be detected?
        """
        PlottableAxis.__init__(self, **kw)
        self.mu = mu
        self.threshold = threshold
        self.label = rf"Transit Depth\nof 1 Scale Height\n for $\mu$={mu} Atmosphere"

    def value(self):
        return self.panel.pop.transmission_signal(mu=self.mu, threshold=self.threshold)


class TransmissionSNR(StellarBrightnessTelescope):
    scale = "log"
    size_normalization = 10

    def __init__(self, mu=2.32, threshold=2, **kw):
        """
        Initialize for a particular mean molecular weight.

        Parameters
        ----------
        mu : float
            Mean molecular weight (default 2.2 for H/He)
        threshold : float
            By how many sigma must the planet mass be detected?
        """
        self.mu = mu
        self.threshold = threshold
        StellarBrightnessTelescope.__init__(self, **kw)
        self.kw = kw

    def setup_label(self):
        """
        How should this plottable be labled on an axis?
        """
        # define the label, based on the wavelength and telescope
        w = self.wavelength.to(u.micron).value
        self.label = rf"S/N for Transit Depth\nof 1 Scale Height\n for $\mu$={self.mu} Atmosphere\nat $\lambda={w}\mu$m (R={self.R})"

    def value(self):
        """
        What data value to plot?
        """
        return self.panel.pop.transmission_snr(
            mu=self.mu, threshold=self.threshold, **self.kw
        )


class Reflection(Depth):
    def __init__(self, albedo=0.1, **kw):
        """
        Initialize for a particular albedo.

        Parameters
        ----------
        albedo : float
            What fraction of starlight does the planet reflect?
        """
        PlottableAxis.__init__(self, **kw)
        self.albedo = albedo
        self.label = f"Reflected Light\nEclipse Depth\n({albedo:.0%} albedo)"

    def value(self):
        return self.panel.pop.reflection_signal(self.albedo)


class ReflectionSNR(StellarBrightnessTelescope):
    scale = "log"
    size_normalization = 10

    def __init__(self, albedo=0.1, **kw):
        """
        Initialize for a particular albedo.

        Parameters
        ----------
        albedo : float
            What fraction of starlight does the planet reflect?
        """
        self.albedo = albedo
        StellarBrightnessTelescope.__init__(self, **kw)
        self.kw = kw

    def setup_label(self):
        """
        How should this plottable be labled on an axis?
        """
        # define the label, based on the wavelength and telescope
        w = self.wavelength.to(u.micron).value
        self.label = f"S/N for Reflected Light\nEclipse Depth\n({self.albedo:.0%} albedo)\nat $\lambda={w}\mu$m (R={self.R})"

    def value(self):
        """
        What data value to plot?
        """
        return self.panel.pop.reflection_snr(albedo=self.albedo, **self.kw)


class Emission(Depth):
    source = "emission_signal"

    def __init__(self, wavelength=5 * u.micron, **kw):
        """
        Initialize for a particular wavelength, because the
        eclipse depth will depend on the thermal emission
        spectrum of the planet and star.
        """
        PlottableAxis.__init__(self, **kw)
        self.wavelength = wavelength
        self.label = f"Thermal Emission\nEclipse Depth\nat $\lambda={self.wavelength.to(u.micron).value}\mu m$"

    def value(self):
        return self.panel.pop.emission_signal(self.wavelength)


class EmissionSNR(StellarBrightnessTelescope):
    scale = "log"
    size_normalization = 10

    def __init__(self, **kw):
        """
        Initialize for a particular albedo.

        Parameters
        ----------
        wavelength : astropy.units.quantity.Quantity
            The wavelength at which we want the
            eclipse depth to be computed.
        """
        StellarBrightnessTelescope.__init__(self, **kw)
        self.kw = kw

    def setup_label(self):
        """
        How should this plottable be labled on an axis?
        """
        # define the label, based on the wavelength and telescope
        w = self.wavelength.to(u.micron).value
        self.label = f"S/N for Thermal Emission\nEclipse Depth\nat $\lambda={self.wavelength.to(u.micron).value}\mu m$\n(R={self.R})"

    def value(self):
        """
        What data value to plot?
        """
        return self.panel.pop.emission_snr(**self.kw)


# FIXME
# -have a clearer way of tracking the
