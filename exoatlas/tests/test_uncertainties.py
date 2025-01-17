from .setup_tests import *

from exoatlas.imports import *
from exoatlas.populations import *
from exoatlas.populations.pineda_skew import *


def test_skew(N_data=25):
    mu = np.zeros(N_data)
    sigma_lower = 10 ** np.linspace(-1, 1, N_data)
    sigma_upper = 10 ** np.linspace(1, -1, N_data)
    x = make_skew_samples_from_lowerupper(
        mu=mu, sigma_lower=sigma_lower, sigma_upper=sigma_upper
    )
    measured_lower, measured_upper = np.percentile(
        x,
        np.array([0.5 - gaussian_central_1sigma / 2, 0.5 + gaussian_central_1sigma / 2])
        * 100,
        axis=-1,
    )

    x_axis = np.log10(sigma_upper / sigma_lower)  # np.arange(N_data) #
    expected_kw = dict(color="gray", zorder=-1, alpha=0.3)
    fi, ax = plt.subplots(2, 1, sharex=True)
    plt.sca(ax[0])
    plt.plot(
        x_axis,
        mu,
        linestyle="-",
        label=r"requested $\mu^{+\sigma_{upper}}_{-\sigma_{lower}}$",
        **expected_kw
    )
    plt.plot(x_axis, mu - sigma_lower, linestyle="--", **expected_kw)
    plt.plot(x_axis, mu + sigma_upper, linestyle="--", **expected_kw)
    if N_data == 1:
        w = 1
    else:
        w = np.median(np.gradient(x_axis)) * 0.9
    plt.violinplot(
        x.T,
        positions=x_axis,
        widths=w,
        showextrema=False,
        quantiles=[
            [0.5 - gaussian_central_1sigma / 2, 0.5, 0.5 + gaussian_central_1sigma / 2]
        ]
        * N_data,
    )
    plt.legend(frameon=False)
    plt.ylabel("simulated distributions")

    plt.sca(ax[1])
    plt.plot(x_axis, (mu - measured_lower) / sigma_lower, marker="o", label="upper")
    plt.plot(x_axis, (measured_upper - mu) / sigma_upper, marker="o", label="lower")
    plt.legend(frameon=False)
    plt.xlabel(r"$\log_{10}(\sigma_{upper}/\sigma_{lower})$")
    plt.ylabel(r"$\sigma_{68\%, measured}/\sigma_{injected}$")


def test_make_astropy_distribution(key="radius"):
    e = TransitingExoplanets()
    mu = e.get_values_from_table(key)
    lower, upper = e.get_uncertainty_lowerupper_from_table(key)
    inject_kw = dict(color="gray", linewidth=5, alpha=0.2)
    fi, ax = plt.subplots(2, 1, sharex=True)
    plt.sca(ax[0])
    plt.plot(-lower, **inject_kw)
    plt.plot(upper, **inject_kw)
    plt.title(key)
    plt.ylabel("upper + lower\n" + r"$\sigma_{samples}$ and $\sigma_{table}$")
    d = e.get_values_from_table(key, distribution=True)
    sample_lower, sample_upper = d.pdf_percentiles([15.8, 84.2])
    sample_kw = dict(color="black")
    plt.plot(sample_lower - mu, **sample_kw)
    plt.plot(sample_upper - mu, **sample_kw)
    plt.sca(ax[1])
    plt.plot((sample_lower - mu) / lower, **sample_kw)
    plt.plot((sample_upper - mu) / upper, **sample_kw)
    plt.ylabel("upper + lower\n" + r"$\sigma_{samples}/\sigma_{table}$")


def test_uncertainties():
    """
    Can we estimate uncertainties for lots of quantities?
    """

    p = SolarSystem()

    uncertainty = p.get_uncertainty("radius")
    assert np.all(uncertainty == 0 * u.Rearth)

    upper, lower = p.get_uncertainty_lowerupper("radius")
    assert np.all(lower == 0 * u.Rearth)
    assert np.all(upper == 0 * u.Rearth)

    bad = p.get_uncertainty("distance")
    assert np.all(np.isnan(bad))


if __name__ == "__main__":  # pragma: no cover
    outputs = {k.split("_")[-1]: v() for k, v in locals().items() if "test_" in k}
