from exoatlas.models import *
from exoatlas.imports import *

def test_seager():
    plt.cla()
    plot_both_seager()

def test_zeng():
    plt.cla()
    plot_three_zeng()

def test_chen():
    plt.cla()
    plot_chen()

if __name__ == '__main__':
    outputs = {k.split('_')[-1]:v()
               for k, v in locals().items()
               if 'test_' in k}
