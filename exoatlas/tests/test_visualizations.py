from .setup_tests import *

from exoatlas import *
from exoatlas.imports import *
from exoatlas.visualizations import *

import exoatlas as ex
import matplotlib.pyplot as plt
from exoatlas.visualizations.maps.preset_maps import preset_maps


def test_plottables():
    e = TransitingExoplanets()[:3]
    for k, v in preset_plottables.items():
        x = v()
        print(k)
        print(x.kw)
        print(x.label)
        print(x.value(e))
        print()


def test_maps():
    pops = {}
    pops["exo"] = TransitingExoplanets()
    pops["solarsystem"] = SolarSystem()
    for k, p in preset_maps.items():
        print(k, p)
        p().build(pops=pops)
        plt.title(f"Map={k}")


def test_map_types():
    pops = {}
    pops["solarsystem"] = SolarSystem()
    pops["exo"] = TransitingExoplanets()

    fr = Flux_x_Radius()
    fr.build(pops=pops)

    fr = BubbleMap(
        xaxis=Plottable(source="insolation"), yaxis=Plottable(source="radius")
    )
    fr.build(pops=pops)

    fr = ErrorMap(
        xaxis=Plottable(source="insolation"), yaxis=Plottable(source="radius")
    )
    fr.build(pops=pops)


def test_galleries():
    pops = {}
    pops["solarsystem"] = SolarSystem()
    FourPanelTransitGallery().build(pops)
    PlanetGallery().build(pops)


if __name__ == "__main__":  # pragma: no cover
    outputs = {k.split("_")[-1]: v() for k, v in locals().items() if "test_" in k}
