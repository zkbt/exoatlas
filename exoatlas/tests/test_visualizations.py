from .setup_tests import *

from exoatlas import * 
from exoatlas.imports import *
from exoatlas.visualizations import * 

import exoatlas as ex
import matplotlib.pyplot as plt
from exoatlas.visualizations.panels.preset_panels import predefined_panels


def test_panels():
    pops = {}
    pops["solarsystem"] = SolarSystem()
    for p in predefined_panels:
        print(p)
        plt.figure()
        p().build(pops=pops)
        plt.close()


def test_panel_types():
    pops = {}
    pops["solarsystem"] = SolarSystem()

    fr = FluxRadius()
    fr.build(pops=pops)

    fr = BubblePanel("insolation", "radius")
    fr.build(pops=pops)

    fr = ErrorPanel("insolation", "radius")
    fr.build(pops=pops)


def test_multipanel_presets():
    with mock.patch("builtins.input", return_value=""):
        t = TransitingExoplanets()
        s = SolarSystem()

    observable_summary([t, s])
    physical_summary([t, s])


def test_colors():
    with mock.patch("builtins.input", return_value=""):
        t = TransitingExoplanets()
        s = SolarSystem()

    physical_summary([t, s])
    t.color = None
    s.color = None
    physical_summary([t, s])


def test_fourpanels():
    pops = {}
    pops["solarsystem"] = SolarSystem()
    f = FourPanels()
    f.build(pops)


if __name__ == "__main__":  # pragma: no cover
    outputs = {k.split("_")[-1]: v() for k, v in locals().items() if "test_" in k}
