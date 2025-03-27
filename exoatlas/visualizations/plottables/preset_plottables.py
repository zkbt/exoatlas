from ...imports import *
from .plottable import *
from ...telescopes import *


class Flux(Plottable):
    source = "relative_insolation"
    label = "Planet Bolometric Flux Received\n(relative to Earth)"
    scale = "log"
    lim = [6e4, 2e-4]


class Teq(Plottable):
    source = "teq"
    scale = "log"
    lim = [200, 2000]

    def _update_label(self):

        self.label = (
            "Planet Equilibrum Temperature (K, f={f}, {albedo:.0%} albedo)".format(
                **self.kw
            )
        )


class CumulativeXUVFlux(Flux):
    source = "relative_cumulative_xuv_insolation"
    label = "Time-Integrated XUV Flux Received\n(relative to Earth)"


class ImpactVelocity(Plottable):
    source = "impact_velocity"
    label = "Estimated Impact Velocity (km/s)"
    scale = "log"


class Radius(Plottable):
    source = "radius"
    label = "Planet Radius (Earth radii)"
    scale = "log"
    lim = [0.3, 30]


class Mass(Plottable):
    source = "mass"
    label = "Planet Mass\n(Earth masses)"
    scale = "log"
    lim = [0.03, 3000]


class SemimajorAxis(Plottable):
    source = "semimajoraxis"
    label = "Semimajor Axis (AU)\n"
    scale = "log"
    lim = [0.001, 1000]


class AngularSeparation(Plottable):
    source = "angular_separation"
    label = "Angular Separation (arcsec)\n"
    scale = "log"
    lim = [0.001, 10]


class Contrast(Plottable):
    source = "imaging_contrast"
    label = "Planet-to-Star Contrast"
    scale = "log"
    lim = [1e-10, 1e-3]


class KludgedMass(Plottable):
    source = "kludge_mass"
    label = "Planet Mass or msini\n(Earth masses)"
    scale = "log"
    lim = [0.03, 4131]


class StellarTeff(Plottable):
    source = "stellar_teff"
    label = "Stellar Temperature (K)"
    scale = "linear"
    lim = [2000, 7000]


class StellarLuminosity(Plottable):
    source = "stellar_luminosity"
    label = "Stellar Luminosity (L_$\odot$)"
    scale = "log"
    lim = [None, None]


class Distance(Plottable):
    source = "distance"
    label = "Distance\n(parsecs)"
    scale = "log"
    lim = [5, 1000]


class Distance(Plottable):
    source = "distance"
    label = "Distance\n(parsecs)"
    scale = "log"
    lim = [5, 1000]


class EscapeVelocity(Plottable):
    source = "escape_velocity"
    label = "Escape\nVelocity\n(km/s)"
    scale = "log"
    lim = [2, 500]


class EscapeParameter(Plottable):
    source = "escape_parameter"
    label = "$\lambda = E_{grav}/E_{thermal}$"
    scale = "log"
    lim = [None, None]


class Density(Plottable):
    source = "density"
    label = "Planet Density\n(g/cm$^3$)"
    scale = "log"
    lim = [0.01, 100]


class StellarRadius(Plottable):
    source = "stellar_radius"
    label = "Stellar Radius\n(solar radii)"
    scale = "linear"
    lim = [0.0, 2.0]


class Period(Plottable):
    source = "period"
    label = "Period (days)\n"
    scale = "log"
    lim = [0.15, 365]


class Gmag(Plottable):
    source = "magnitude_gaia"
    label = "G (magnitude)\n"
    scale = "linear"
    lim = [3.5, 14.5]


class Depth(Plottable):
    source = "transit_depth"
    label = "Transit Depth\nof Planet"
    scale = "log"
    lim = [2e-6, 2e-1]


class StellarBrightness(Plottable):
    source = "stellar_brightness"
    scale = "log"
    lim = [1e2, 1e8]
    unit = u.Unit("ph s^-1 m^-2 micron^-1")

    def _update_label(self):
        self.label = f'Stellar Brightness at Earth at $\lambda={self.kw['wavelength'].to("micron").value:.1f}\mu$m\n({self.unit.to_string("latex_inline")})'


class StellarBrightnessTelescope(Plottable):
    source = "stellar_brightness_in_telescope_units"
    scale = "log"
    lim = [1e-3, 1e3]

    def _update_label(self):
        self.label = (
            f"Stellar Brightness at Earth at $\lambda={w}\mu$m\n({self.unit_string})"
        )

    def __init__(self, telescope_name="JWST", **kw):
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

        # do all the basic telescope + wavelength setup
        self.setup_telescope(telescope_name=telescope_name, **kw)

        # specify the unit to be used for plotting
        self.setup_unit()

        # initialize the basic plottable axis
        Plottable.__init__(self, **kw)

    def setup_telescope(self, telescope_name=None, **kw):
        """
        Setup the basics of the telescope, wavelength, and units.
        """

        # break if there's not telescope
        assert telescope_name is not None

        # keep track of the telescope (if any)
        self.telescope_name = telescope_name

        # define a unit of collecting power for the telescope
        self.telescope_unit = define_telescope_unit_by_name(self.telescope_name, **kw)

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
        self.unit_string = self.unit.to_string("latex_inline")

    def _update_label(self):
        """
        How should this plottable be labled on an axis?
        """
        # define the label, based on the wavelength and telescope
        w = self.wavelength.to(u.micron).value
        self.label = (
            f"Stellar Brightness at Earth at $\lambda={w}\mu$m\n({self.unit_string})"
        )


class DepthSNR(StellarBrightnessTelescope):
    source = "depth_snr"
    scale = "log"
    lim = [1e-3, 1e3]

    def setup_unit(self):
        # kludge to undo telescope brightness unit
        self.unit = u.Unit()
        self.unit_string = ""

    def _update_label(self):
        """
        How should this plottable be labled on an axis?
        """
        # define the label, based on the wavelength and telescope
        w = self.wavelength.to(u.micron).value
        self.label = f"S/N for Transit Depth\n{self.telescope_unit}"


class Transmission(Depth):
    source = "transmission_signal"

    def _update_label(self):
        mu = self.kw["mu"]
        self.label = f"Transit Depth\nof 1 Scale Height\n for $\mu$={mu} Atmosphere"


class TransmissionSNR(DepthSNR):
    scale = "log"
    size_normalization = 10
    source = "transmission_snr"

    def _update_label(self):
        mu = self.kw["mu"]
        w = self.wavelength.to(u.micron).value
        R = self.R
        self.label = f"S/N for Transit Depth\nof 1 Scale Height\n for $\mu$={mu} Atmosphere\n{self.telescope_unit}"


class Reflection(Depth):
    source = "reflection_signal"

    def _update_label(self):
        albedo = self.kw["albedo"]
        self.label = f"Reflected Light\nEclipse Depth\n({albedo:.0%} albedo)"


class ReflectionSNR(DepthSNR):
    scale = "log"
    size_normalization = 10
    source = "reflection_snr"

    def _update_label(self):
        albedo = self.kw["albedo"]
        w = self.wavelength.to(u.micron).value
        R = self.R
        self.label = f"S/N for Reflected Light\nEclipse Depth\n({albedo:.0%} albedo)\n{self.telescope_unit}"


class Emission(Depth):
    source = "emission_signal"

    def _update_label(self):
        albedo = self.kw["albedo"]
        f = self.kw["f"]
        w = self.kw["wavelength"].to_value("micron")
        self.label = f"Thermal Eclipse Depth\nat $\lambda={w}\mu m$\n(f={f:.2f}, {albedo:.0%} albedo)"


class EmissionSNR(DepthSNR):
    scale = "log"
    size_normalization = 10
    source = "emission_snr"

    def _update_label(self):
        albedo = self.kw["albedo"]
        f = self.kw["f"]
        w = self.wavelength.to(u.micron).value
        R = self.R
        self.label = f"S/N for Thermal\nEclipse Depth\n( f={f:.2f}, {albedo:.0%} albedo)\n{self.telescope_unit}"


class RightAscension(Plottable):
    source = "ra"
    label = "Right Ascension (hours)"
    scale = "linear"
    lim = [24, 0]
    unit = u.hourangle


class Declination(Plottable):
    source = "dec"
    label = "Declination (˚)"
    scale = "linear"
    lim = [-90, 90]
    unit = u.deg


preset_plottables = {}
local_variables = dict(**locals())
for k, v in local_variables.items():
    try:
        assert issubclass(v, Plottable)
        preset_plottables[k] = v
    except (AssertionError, TypeError):
        continue
preset_plottables.pop("Plottable")
