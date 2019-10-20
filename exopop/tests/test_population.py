from exopop.populations.Population import *
from exopop.imports import *

def test_populationfromstandard():

    fake = Table({x:[0]*3 for x in necessary_columns}, masked=True)
    p = PopulationFromStandard(standard=fake, label='fake')
    p.validate_columns()
    return p

if __name__ == '__main__':
    outputs = {k.split('_')[-1]:v()
               for k, v in locals().items()
               if 'test_' in k}
