from .setup_tests import *

from exoatlas import *


def test_airmass():
    e = TransitingExoplanets()

    positions = e.altaz()

    fi, ax = plt.subplots(2, 1, constrained_layout=True, figsize=(8, 8))
    plt.sca(ax[0])

    kw = dict(c=e.airmass(), cmap="magma", vmax=2)
    plt.title(f"{positions.obstime.iso}\n{positions.location}")
    plt.scatter(e.ra(), e.dec(), **kw)
    plt.xlabel("Right Ascension (degrees)")
    plt.ylabel("Declination (degrees)")
    plt.colorbar()
    plt.sca(ax[1])

    plt.scatter(positions.az, positions.alt, **kw)
    plt.xlabel("Azimuth (degrees)")
    plt.ylabel("Altitude (degrees)")
    plt.colorbar()
    plt.xlim(0, 360)
    plt.ylim(-90, 90)


def test_show_upcoming_transits():
    e = TransitingExoplanets()
    nearby = e[e.distance() < 15 * u.pc]
    transits = nearby.show_upcoming_transits(window=10 * u.day)
