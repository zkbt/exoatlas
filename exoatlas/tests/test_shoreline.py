from .setup_tests import *
from exoatlas import *
from exoatlas.populations.calculations.shoreline import *


def test_probability_of_atmosphere():
    # need to change for other people to be able to run tests
    posterior_filename = "/Users/zabe0091/Dropbox/zach/code/exoatlas/sandbox/shoreline/all-any-uncertainties=True-numpyro.nc"
    # need to add automatic download for posteriors!
    s = SolarSystem()
    import arviz as az

    posterior = az.from_netcdf(posterior_filename)
    shoreline = Shoreline(posterior)

    p = s.probability_of_atmosphere(shoreline=shoreline)
    p_uncertainty = s.probability_of_atmosphere_uncertainty(shoreline=shoreline)

    plt.figure(figsize=(8, 3))
    plt.errorbar(
        s.semimajoraxis(), p, p_uncertainty, linewidth=0, elinewidth=2, marker="o"
    )
    plt.xscale("log")
    plt.xlabel("Planet Semimajor Axis (AU)")
    plt.ylabel("Probability of Atmosphere")
