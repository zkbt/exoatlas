import exopop as ex
import matplotlib.pyplot as plt

def test_():
    pops = {}
    pops['solarsystem'] = ex.SolarSystem()
    for p in [  ex.MassRadius,
                ex.FluxRadius,
                ex.StellarRadius,
                ex.DistanceRadius]:
        plt.figure()
        p().build(pops=pops)
        plt.draw()

def test_fourpanels():
    pops = {}
    pops['solarsystem'] = ex.SolarSystem()
    f = ex.FourPanels()
    f.build(pops)

if __name__ == '__main__':
    outputs = {k.split('_')[-1]:v()
               for k, v in locals().items()
               if 'test_' in k}
