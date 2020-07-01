from exoatlas.imports import *

def test_directory():
    print(directories['data'])

def test_name2color():
    one = name2color('black')
    another = name2color('#000000')
    assert(one == another)

if __name__ == '__main__': # pragma: no cover
    outputs = {k.split('_')[-1]:v()
               for k, v in locals().items()
               if 'test_' in k}
