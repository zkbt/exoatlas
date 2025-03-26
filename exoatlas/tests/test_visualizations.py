from .setup_tests import *

from exoatlas import *
from exoatlas.imports import *
from exoatlas.visualizations import *

import exoatlas as ex
import matplotlib.pyplot as plt
from exoatlas.visualizations.panels.preset_panels import preset_panels


def test_plottables():
    e = TransitingExoplanets()[:3]
    for k, v in preset_plottables.items():
        x = v()
        print(k)
        print(x.kw)
        print(x.label)
        print(x.value(e))
        print()


def test_panels():
    pops = {}
    pops["exo"] = TransitingExoplanets()
    pops["solarsystem"] = SolarSystem()
    for k, p in preset_panels.items():
        print(k, p)
        p().build(pops=pops)
        plt.title(f"Panel={k}")


def test_panel_types():
    pops = {}
    pops["solarsystem"] = SolarSystem()
    pops["exo"] = TransitingExoplanets()

    fr = Flux_x_Radius()
    fr.build(pops=pops)

    fr = BubblePanel(
        xaxis=Plottable(source="insolation"), yaxis=Plottable(source="radius")
    )
    fr.build(pops=pops)

    fr = ErrorPanel(
        xaxis=Plottable(source="insolation"), yaxis=Plottable(source="radius")
    )
    fr.build(pops=pops)


def test_galleries():
    pops = {}
    pops["solarsystem"] = SolarSystem()
    TransitGallery().build(pops)
    physical_summary().build(pops)


if __name__ == "__main__":  # pragma: no cover
    outputs = {k.split("_")[-1]: v() for k, v in locals().items() if "test_" in k}
