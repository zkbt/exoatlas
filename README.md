# exoplanet-atlas
Tools for compiling and plotting populations of transiting exoplanets. This package still a bit of a work in progress, but can be generally useful for downloading and working with exoplanet populations. For draft documentation explaing how to use `exoatlas`, please [read the documentation](https://zkbt.github.io/exoplanet-atlas/build/html/index.html).

### Installation
If you want the most recent stable version, the simplest way is to install it from PyPI directly via `pip` from any UNIX prompt:
```bash
pip install exoplanet-atlas
```

Or, if you want the very-most-up-to-date version, you can install directly from this repository via:
```bash
pip install git+https://github.com/zkbt/exoplanet-atlas
```

Or, if you want to be able to modify the code yourself, please also feel free to fork/clone this repository onto your own computer and install directly from that editable package. For example, this might look like:
```bash
git clone https://github.com/zkbt/exoplanet-atlas.git
cd exoplanet-atlas
pip install -e .
```
The `pip install -e .` command will link the installed version of the package to the directory of your local repository. Changes you make to the code in that directory should be reflected in the version Python sees when it tries to `import exoatlas`.

### Usage
Here's a quick preview:


```python
# import some population definitions and plotting tool
from exoatlas import *

# create a dictionary of populations
exo = TransitingExoplanets()
solar = SolarSystem()
pops = {'solar':solar,
        'exo':exo}

# use a default visualization to summarize these populations
physical_summary(pops)
```

### Authors
This toolkit was made by [Zach Berta-Thompson](http://casa.colorado.edu/~bertathompson/). It relies heavily on the incredible work done by the folks over at the [NASA Exoplanet Archive](https://exoplanetarchive.ipac.caltech.edu), and their generously designed API.
