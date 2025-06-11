from exoatlas import *
from exoatlas.visualizations import *

e = TransitingExoplanets()
s = SolarSystem()
h = e["HD209458b"]

PlanetGallery().build([e, s, h])
plt.savefig("paper/joss-exoatlas-example.png")
