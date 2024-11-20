from .setup_tests import *

from exoatlas.imports import *
from exoatlas.populations import *


def test_population():
    """
    Can we make a population from scratch from a table?
    """

    fake = QTable({x: [0] * 3 for x in core_descriptions}, masked=True)
    p = Population(standard=fake, label="fake")
    p._validate_columns()


def test_solarsystem():
    """
    Can we make a population of Solar System planets?
    """
    for p in [
        SolarSystem,
        SolarSystemDwarfPlanets,
        SolarSystemMoons,
        SolarSystemMinorPlanets,
    ]:
        p()._validate_columns()


def test_exoplanets():
    """
    Can we make a population of confirmed exoplanets?
    """
    p = Exoplanets(remake=False)
    p.load_individual_references()
    p._validate_columns()

    # check that the background reference population for Exoplanets filters too
    f = p[:10]
    for x in f.hostname:
        assert x in f.individual_references.hostname
    for x in f.individual_references.hostname:
        assert x in f.hostname


def test_transitingexoplanets():
    """
    Can we make a population of confirmed transiting exoplanets?
    """
    p = TransitingExoplanets(remake=False)
    p._validate_columns()


def test_subsets():
    """
    Can we make a population of confirmed Kepler planets?
    """
    for x in [Kepler, NonKepler, TESS, NonTESS, Space, Ground, GoodMass, BadMass]:
        p = x(remake=False)
        p._validate_columns()

    with pytest.raises(NotImplementedError):
        q = TransitingExoplanetsSubset()


def test_indexing():
    """
    Can we make a population of confirmed transiting exoplanets?
    """
    p = Exoplanets(remake=False)

    # try different subsets
    b = p[5]
    c = p[0:10]
    with np.errstate(invalid="ignore"):
        d = p[p.stellar_radius < 1.0 * u.Rsun]
    e = p["GJ 1214b"]
    f = p[["GJ 1214b", "LHS 1140b", "GJ 1132b"]]
    g = p[p.discovery_facility == "Kepler"]
    h = p["TRAPPIST-1b"]
    i = p["TRAPPIST-1"]
    j = e + h
    k = f - e
    l = p.create_subset_by_name("GJ1214b")
    m = p.create_subset_by_hostname("GJ1214")
    coordinates = SkyCoord(e.ra, e.dec)
    n = p.create_subset_by_position(coordinates)


def test_table():
    """
    Can we make a custom astropy table out of a Population?
    """

    p = SolarSystem()
    p.create_table()


def test_attributes():
    """
    Test the interaction with attributes.
    """
    p = SolarSystem()

    print(p.color)
    p.alpha = 0.5
    assert p._plotkw["alpha"] == 0.5

    for k in basic_columns:
        getattr(p, k)


def test_transiting(planet="GJ1214b"):
    with mock.patch("builtins.input", return_value=""):
        t = TransitingExoplanets()
        p = t[planet]

    summarize_planet(p)

    p.transmission_signal(mu=5, threshold=3)
    p.transmission_snr(
        mu=5, threshold=3, telescope_name="JWST", wavelength=3 * u.micron
    )

    p.emission_signal(wavelength=5 * u.micron)
    p.emission_snr(wavelength=5 * u.micron, telescope_name="JWST")

    p.reflection_signal(albedo=0.5)
    p.reflection_snr(albedo=0.5, wavelength=5 * u.micron, telescope_name="JWST")
