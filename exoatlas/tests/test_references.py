from .setup_tests import *

from exoatlas import *


def test_references():
    e = Exoplanets()
    e.load_individual_references()
    planet = "GJ1132b"
    p = e["GJ1132b"]
    assert len(p.individual_references.standard) > len(p.standard)
