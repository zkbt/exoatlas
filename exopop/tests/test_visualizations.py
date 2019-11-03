import exopop as ex

def test_():
    pops = {}
    pops['solarsystem'] = ex.SolarSystem()
    f = ex.MultiPanelPlot(pops,
                          panels=[ex.MassRadius, ex.FluxRadius])
    f.build()

def test_fourpanels():
    pops = {}
    pops['solarsystem'] = ex.SolarSystem()
    f = ex.FourPanels(pops)
    f.build()

if __name__ == '__main__':
    outputs = {k.split('_')[-1]:v()
               for k, v in locals().items()
               if 'test_' in k}
