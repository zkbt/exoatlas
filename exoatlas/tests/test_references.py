from .setup_tests import *

from exoatlas import *


def test_load_references():
    e = Exoplanets()
    e.load_individual_references()
    planet = "GJ1132b"
    p = e["GJ1132b"]
    assert len(p.individual_references.standard) > len(p.standard)


def test_check_and_update_references():
    e = Exoplanets()
    e.load_individual_references()
    e.check_individual_references(planets="HD189733b")
    e.update_reference(references="Ivshina + Winn 2022")
    assert e["HD189733b"].period_reference[0] == "Ivshina + Winn 2022"
