from exoatlas.models import *
from exoatlas.imports import *


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
    plot_chen()


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
