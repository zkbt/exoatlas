from .setup_tests import *

from exoatlas import *


def test_exoplanets_downloader():
    planet_name = "GJ 1132 b"

    for t in ["ps", "pscomppars"]:
        fresh = ExoplanetArchiveDownloader(which_table=t, planet=planet_name).get(
            remake=True
        )
        loaded = ExoplanetArchiveDownloader(which_table=t, planet=planet_name).get(
            remake=False
        )
        assert np.all(loaded["ra"] == fresh["ra"])
