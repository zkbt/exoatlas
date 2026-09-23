from .setup_tests import *
from exoatlas.models import *
from exoatlas.imports import *


def test_tang():
    t = Tang()
    t.plot()


def test_baraffe():
    b = Baraffe()
    b.plot()
    b.imshow()


def test_spectrum():
    N = 100
    w = np.linspace(1, 1000, N) * u.nm
    f = np.random.normal(100, 1, N) * u.W / u.m**2 / u.nm
    s = Spectrum(w, f)
    s.integrate()

    s = Thermal()
    s.integrate()

    s = SumOfThermal()
    s.integrate()


def test_mamajek():
    m = Mamajek()
    m.plot(["logT", "logL", "logR", "logM"])
    logL = np.linspace(-6, 6)
    logT = m.tofrom("logT")("logL")(logL)


def test_seager():
    """
    Test the Seager mass-radius relations.
    """
    plt.cla()
    plot_both_seager()


def test_zeng():
    plt.cla()
    plot_three_zeng()


def test_chen():
    plt.cla()
    plot_chen("mass")
    plot_chen("radius", alpha=0.5, linewidth=4)


def test_kopparapu():
    # test a default HZ
    T = np.linspace(2600, 7000)
    f = make_hz()
    f(T)

    # test that it breaks when it should
    with pytest.raises(Exception):
        make_hz("something-that-does-not-exist")


if __name__ == "__main__":  # pragma: no cover
    outputs = {k.split("_")[-1]: v() for k, v in locals().items() if "test_" in k}
