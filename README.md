# exoatlas
Tools for compiling and plotting populations of transiting exoplanets. This package still a bit of a work in progress, but can be generally useful for downloading and working with exoplanet populations. For draft documentation explaing how to use `exoatlas`, please [read the documentation](https://zkbt.github.io/exoatlas/).

### Installation
If you want the most recent stable version, the simplest way is to install it from PyPI directly via `pip` from any UNIX prompt:
```bash
pip install exoatlas
```

Or, if you want the very-most-up-to-date version, you can install directly from this repository via:
```bash
pip install git+https://github.com/zkbt/exoatlas
```

Or, if you want to be able to modify the code yourself, please also feel free to fork/clone this repository onto your own computer and install directly from that editable package. For example, this might look like:
```bash
git clone https://github.com/zkbt/exoatlas.git
cd exoatlas
pip install -e .
```
The `pip install -e .` command will link the installed version of the package to the directory of your local repository. Changes you make to the code in that directory should be reflected in the version Python sees when it tries to `import exoatlas`.

### Usage
Here's a very quick preview:


```python
# import some population definitions and plotting tool
from exoatlas import TransitingExoplanets, SolarSystem
from exoatlas.visualizations import PlanetGallery

# create a dictionary of populations
exo = TransitingExoplanets()
solar = SolarSystem()

# use a default visualization to summarize these populations
PlanetGallery().build([solar, exo])
```
For a slightly less quick preview, which will hopefully entice you to keep reading through the rest of the documentation, please 

### Authors
This toolkit was made by [Zach Berta-Thompson](http://casa.colorado.edu/~bertathompson/). It relies heavily on the incredible work done by the folks over at the [NASA Exoplanet Archive](https://exoplanetarchive.ipac.caltech.edu), and their generously designed API.
