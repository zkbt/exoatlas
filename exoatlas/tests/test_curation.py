from .setup_tests import *

from exoatlas import *


def test_load_references():
    e = Exoplanets()
    e.load_individual_references()
    p = e["GJ1132b"]
    assert len(p.individual_references.standard) > len(p.standard)


def test_check_and_update_references():
    e = Exoplanets()
    e.load_individual_references()
    e.display_individual_references(planets="HD189733b")
    e.update_reference(references="Ivshina + Winn 2022")
    assert (
        e["HD189733b"].get_values_from_table("period_reference")[0]
        == "Ivshina + Winn 2022"
    )


def test_update_values():
    e = Exoplanets()
    e.update_values(
        planets="GJ1214b", mass=1 * u.M_earth, mass_uncertainty=0.1 * u.M_earth
    )
    one = e["GJ1214b"]
    assert one.get("mass") == 1 * u.M_earth
    assert one.get_uncertainty("mass") == 0.1 * u.M_earth
    assert one.get_uncertainty_lowerupper("mass")[0] == 0.1 * u.M_earth
