from .setup_tests import *

from exoplanet_atlas.imports import *
from exoplanet_atlas.populations import TransitingExoplanets
from exoplanet_atlas.whatsup import Plan


def test_whatsup():
    with mock.patch("builtins.input", return_value=""):
        planets = TransitingExoplanets()

    p = Plan(
        planets,  # what should we try?
        start=Time("2020-03-14"),  # when should we start?
        finish=Time("2020-03-14"),  # when should we stop?
        observatory="SBO",  # what observatory are starting from?
        directory="upcoming-transits",
    )  # where should we


if __name__ == "__main__":  # pragma: no cover
    outputs = {k.split("_")[-1]: v() for k, v in locals().items() if "test_" in k}
