from exoatlas.imports import *

def test_directory():
    print(directories['data'])

if __name__ == '__main__':
    outputs = {k.split('_')[-1]:v()
               for k, v in locals().items()
               if 'test_' in k}
