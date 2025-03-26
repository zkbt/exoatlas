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

    t = TransitingExoplanets()[:20]
    Depth_x_Brightness().build([t])
    BubblePanel(xaxis=StellarBrightness, yaxis=Depth).build([t])
    BubblePanel(xaxis=StellarBrightness(wavelength=5 * u.micron), yaxis=Depth).build(
        [t]
    )
    for k in telescope_units:
        BubblePanel(
            xaxis=StellarBrightnessTelescope(telescope_name=k), yaxis=Depth
        ).build([t])
    BubblePanel(
        xaxis=StellarBrightnessTelescope(
            telescope_name="JWST", wavelength=10 * u.micron, R=10
        ),
        yaxis=Depth,
    ).build([t])


if __name__ == "__main__":  # pragma: no cover
    outputs = {k.split("_")[-1]: v() for k, v in locals().items() if "test_" in k}
