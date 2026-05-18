from .setup_tests import *
from exoatlas import *
from exoatlas.visualizations import * 
from exoatlas.calculations.shoreline import *


def test_shoreline_probability_of_atmosphere():
    # need to add automatic download for posteriors!
    s = SolarSystem()
    shoreline = Shoreline()

    p = s.probability_of_atmosphere(shoreline=shoreline)
    p_uncertainty = s.probability_of_atmosphere_uncertainty(shoreline=shoreline)

    plt.figure(figsize=(8, 3))
    plt.errorbar(
        s.semimajoraxis(), p, p_uncertainty, linewidth=0, elinewidth=2, marker="o"
    )
    plt.xscale("log")
    plt.xlabel("Planet Semimajor Axis (AU)")
    plt.ylabel("Probability of Atmosphere")


def test_shoreline_visualizations():

    # create populations
    s = SolarSystem()
    e = TransitingExoplanets()

    # pick the name of the planet to highlight
    planet_name = 'LTT1445Ab'

    # create a subset population to highlight that one planet
    highlight = e[planet_name]

    # define a single panel of the visualization
    m = ShorelineStandardMap()

    # construct a grid of multiple slices
    g = SliceGridGallery(m, N=4)

    # add the planets to the panels
    g.build([s, e, highlight])

    # add embellishments to the plot, including probability 
    g.refine()

    # colorbar for probability
    g.add_colorbar()