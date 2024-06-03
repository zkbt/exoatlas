from .setup_tests import *

from exoatlas.imports import *
from exoatlas.populations.downloaders import *


def test_exoplanets():
    with mock.patch("builtins.input", return_value=""):
        planetary_systems_downloader.get()


# def test_composite():
#    with mock.patch("builtins.input", return_value=""):
#        composite_planetary_systems_downloader.get()


# def test_tess():
#    with mock.patch("builtins.input", return_value=""):
#        toi_exofop.get()


if __name__ == "__main__":  # pragma: no covers
    a = test_exoplanets()
    c = test_composite()
    t = test_tess()
