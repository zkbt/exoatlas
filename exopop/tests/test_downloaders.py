from exopop.populations.downloaders import *

def test_exoplanets():
    return exoplanets.get()

def test_composite():
    return composite_exoplanets.get()

def test_composite():
    return merged_exoplanets.get()


def test_tess():
    return toi_exofop.get()

if __name__ == '__main__':
    a = test_exoplanets()
    c = test_composite()
    t = test_tess()
