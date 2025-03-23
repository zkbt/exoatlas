from .setup_tests import *
from exoatlas.imports import *
from exoatlas import *
from exoatlas.telescopes import * 
from exoatlas.visualizations import *


def test_telescope_units():
    for t in ["Kepler", "TESS", "JWST", "HST"]:
        define_telescope_unit_by_name(t)
        define_telescope_unit_by_name(t, wavelength=0.7 * u.micron)


def test_buckets():
    """
    Make sure that our photon-counting tools plottables work,
    with their various possible telescope units.
    """

    with mock.patch("builtins.input", return_value=""):
        t = TransitingExoplanets()
    DepthBrightness().build([t])
    BubblePanel(StellarBrightness, Depth).build([t])
    BubblePanel(StellarBrightness(5 * u.micron), Depth).build([t])
    for k in telescope_units:
        BubblePanel(StellarBrightnessTelescope(telescope_name=k), Depth).build([t])
    BubblePanel(
        StellarBrightnessTelescope(
            telescope_name="JWST", wavelength=10 * u.micron, R=10
        ),
        Depth,
    ).build([t])


if __name__ == "__main__":  # pragma: no cover
    outputs = {k.split("_")[-1]: v() for k, v in locals().items() if "test_" in k}
