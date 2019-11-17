# exoatlas
Tools for compiling and plotting populations of transiting exoplanets. It's still a bit of a work in progress, but can be generally useful for downloading and working with exoplanet populations.

### Installation
To install, this simplest way is probably simply to install it directly via `pip` from any UNIX prompt:
```
pip install git+https://github.com/zkbt/exoatlas
```

If you want to be able to modify the code yourself, please also feel free to fork/clone this repository onto your own computer and install directly from that editable package. For example, this might look like:
```
git clone https://github.com/zkbt/exoatlas.git
cd exoatlas
pip install -e .
```
This will link the installed version of the `exoatlas` package to your local repository. Changes you make to the code in the repository should be reflected in the version Python sees when it tries to `import exoatlas`.

### Usage
Here's a quick preview:

```
from exoatlas import *
exo = TransitingExoplanets()
solar = SolarSystem()
pops = {'solar':solar,
        'exo':exo}
physical_summary(pops)
```

### Authors
This toolkit was made by [Zach Berta-Thompson](http://casa.colorado.edu/~bertathompson/). It relies heavily on the incredible work done by the folks over at the [NASA Exoplanet Archive](https://exoplanetarchive.ipac.caltech.edu), and their generously designed API.
