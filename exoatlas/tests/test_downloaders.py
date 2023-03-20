from exoatlas.imports import *
from exoatlas.populations.downloaders import *


def test_exoplanets():
    with mock.patch("builtins.input", return_value=""):
        return exoplanets.get()


def test_composite():
    with mock.patch("builtins.input", return_value=""):
        return merged_exoplanets.get()


def test_tess():
    with mock.patch("builtins.input", return_value=""):
        return toi_exofop.get()


if __name__ == "__main__":  # pragma: no cover
    a = test_exoplanets()
    c = test_composite()
    t = test_tess()
