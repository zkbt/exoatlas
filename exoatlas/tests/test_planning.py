from .setup_tests import *

from exoatlas import *


def test_airmass(visualize=False):
    e = TransitingExoplanets()

    positions = e.altaz()
    airmass = e.airmass()


def test_show_upcoming_transits(visualize=False):
    e = TransitingExoplanets()
    nearby = e[e.distance() < 15 * u.pc]
    transits = nearby.show_upcoming_transits(window=10 * u.day, visualize=visualize)
