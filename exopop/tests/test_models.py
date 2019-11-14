from exopop.models import *

def test_seager():
    plot_both_seager()

if __name__ == '__main__':
    outputs = {k.split('_')[-1]:v()
               for k, v in locals().items()
               if 'test_' in k}
