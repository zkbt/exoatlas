from exopop.populations.downloaders import *

def test_exoplanets():
    return exoplanets.get(skip_update=True)

def test_composite():
    return composite_exoplanets.get(skip_update=True)

def test_composite():
    return merged_exoplanets.get(skip_update=True)

def test_tess():
    return toi_exofop.get(skip_update=True)

if __name__ == '__main__':
    a = test_exoplanets()
    c = test_composite()
    t = test_tess()
