# exopop
Tools for compiling and plotting populations of transiting exoplanets. It's still a bit of a work in progress, but can be generally useful for downloading and working with exoplanet populations.

### Installation
To install, this simplest way is probably simply to install it directly via `pip` from any UNIX prompt:
```
pip install git+https://github.com/zkbt/exopop
```

If you want to be able to modify the code youself, please also feel free to fork/clone this repository and install it into you local Python library via the command (from the main `exopop/` repository directory):
```
python setup.py install
```

You probably also need [craftroom](https://github.com/zkbt/craftroom) and [thistothat](https://github.com/zkbt/thistothat).

### Usage
Here's a quick preview:

```
import matplotlib.pyplot as plt
import numpy as np

# pull the confirmed exoplanets from the NASA Exoplanet Archive
from exopop.Confirmed import Confirmed
confirmed = Confirmed()

# there's a standardized table of properties
print(confirmed.standard)

# and there are quantities that can be derived from these
print(confirmed.distance)

# you can make your own plot
plt.scatter(confirmed.planet_mass, confirmed.planet_radius)
plt.ylim(0.5, 30)
plt.xlim(0.5, 1000)
plt.yscale('log')
plt.xscale('log')
plt.show()

# or rely on some prepackaged ones
from exopop.PlotPanels import MassRadius
import matplotlib.pyplot as plt
pops = dict(goodmass=confirmed)
mr = MassRadius(pops=pops)
mr.build()
mr.ax.set_ylim(0.5, 30)
mr.ax.set_xlim(0.5, 1000)
mr.ax.set_yscale('log')
mr.ax.set_xscale('log')
plt.show()
```

### Authors
This toolkit was made by [Zach Berta-Thompson](http://casa.colorado.edu/~bertathompson/). It relies heavily on the incredible work done by the folks over at the [NASA Exoplanet Archive](https://exoplanetarchive.ipac.caltech.edu), and their generously designed API.
