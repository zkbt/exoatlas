import exopop as ex
import matplotlib.pyplot as plt

def test_somepanels():
    pops = {}
    pops['solarsystem'] = ex.SolarSystem()
    for p in [  ex.MassRadius,
                ex.FluxRadius,
                ex.StellarRadius,
                ex.DistanceRadius]:
        plt.figure()
        p().build(pops=pops)
        plt.draw()

def test_presets():
    pops = {}
    pops['solarsystem'] = ex.SolarSystem()
    ex.physical_summary(pops)
    ex.observable_summary(pops)


def test_fourpanels():
    pops = {}
    pops['solarsystem'] = ex.SolarSystem()
    f = ex.FourPanels()
    f.build(pops)

if __name__ == '__main__':
    outputs = {k.split('_')[-1]:v()
               for k, v in locals().items()
               if 'test_' in k}
