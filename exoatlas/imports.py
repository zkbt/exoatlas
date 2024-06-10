from .version import *

# imports that are need by many exoatlas subsections
from ast import Import
import os, sys, time, shutil, warnings, copy, glob
from tqdm import tqdm

# (possibly different on Mac, Linux, Windows, even for Python versions 3.8-3.12)
try:
    from importlib.resources import files
except (ModuleNotFoundError, AttributeError, ImportError):
    from importlib_resources import files

code_directory = files(import_name)


import numpy as np, matplotlib.pyplot as plt, matplotlib.animation as animation
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter, LogLocator


from astropy.utils.exceptions import AstropyDeprecationWarning

warnings.simplefilter("ignore", category=AstropyDeprecationWarning)

# this function downloads a file and returns its filepath
from astropy.utils.data import download_file
from astropy.io import ascii
from astropy.table import Table, vstack, join, setdiff
from astropy.visualization import quantity_support
from astropy.coordinates import SkyCoord

quantity_support()

# some general custom utilities from Zach
from .talker import Talker

# units and constants from astropy
import astropy.units as u, astropy.constants as con
from astropy.time import Time


def mkdir(path):
    """A mkdir that doesn't complain if it fails."""
    try:
        os.mkdir(path)
    except:
        pass


import matplotlib.colors as co


def name2color(name):
    """
    Return the 3-element RGB array of a given color name.

    Parameters
    ----------
    name : str
            The name ('black') or hex ('#000000')
            for a particular color.

    Returns
    -------
    rgb : tuple
            Three-elements array of floats,
            expressing brightness in RGB.
    """
    if "#" in name:
        h = name
    else:
        h = co.cnames[name].lower()
    return co.hex2color(h)


# create a directory structure ()
try:
    # search for an environment variable
    base = os.getenv("exoatlas_DATA")
    assert base is not None
except AssertionError:
    # otherwise put it in the local directory
    cwd = os.getcwd()
    base = os.path.join(cwd, "downloads-for-exoatlas")
mkdir(base)


def locate_local_data():
    print("💾 `exoatlas` archive data will be stored in:")
    print(base)


directories = dict(data=os.path.join(base, "data/"))
for k in directories.keys():
    mkdir(directories[k])


def reset_standardized_data():
    files = glob.glob(os.path.join(directories["data"], "standardized-*"))
    if "y" in input(f"Are you sure you want to wipe {files}? [y/N]"):
        for f in files:
            os.remove(f)


def reset_local_data():
    if "y" in input(
        "Are you sure you want to wipe all " "local exoatlas data files? [y/N]"
    ):
        shutil.rmtree(directories["data"])
        mkdir(directories["data"])
        print(f"Removed all local data from {directories['data']}")


# some kludge for dealing with Python 3 vs 2?
try:
    FileNotFoundError
except NameError:
    FileNotFoundError = IOError


def clean(s):
    """
    A wrapper function to clean up complicated strings.
    """
    bad = """ !@#$%^&*()+-_'",./<>?"""
    cleaned = str(s) + ""
    for c in bad:
        cleaned = cleaned.replace(c, "")
    return cleaned


def time_from_modified(filename):
    """
    How long ago was this file last modified?
    """
    try:
        dt = Time.now().unix - os.path.getmtime(filename)
        return dt / 60 / 60 / 24
    except FileNotFoundError:
        return np.inf


def check_if_needs_updating(filename, maximum_age=1.0):
    """
    Do a (possibly interactive) check to see if this file
    is so old that it needs to be updated.

    Returns
    -------
    old : bool
            True if the file is so old it needs updating
            False if the file doesn't need updating
    """

    # how long ago was the data updated?
    dt = time_from_modified(filename)

    old = False
    if dt == np.inf:
        old = True
    elif dt > maximum_age:
        print(f"{filename} is {dt:.3f} days old.")
        old = "y" in input("Should it be updated? [y/N]").lower()

    return old


def one2another(bottom="white", top="red", alphabottom=1.0, alphatop=1.0, N=256):
    """
    Create a cmap that goes smoothly (linearly in RGBA) from "bottom" to "top".
    """
    try:
        rgb_bottom = bottom[0:3] + 0
    except:
        rgb_bottom = name2color(bottom or "white")
    try:
        rgb_top = top[0:3] + 0
    except:
        rgb_top = name2color(top or "black")

    r = np.linspace(rgb_bottom[0], rgb_top[0], N)
    g = np.linspace(rgb_bottom[1], rgb_top[1], N)
    b = np.linspace(rgb_bottom[2], rgb_top[2], N)
    a = np.linspace(alphabottom, alphatop, N)
    colors = np.transpose(np.vstack([r, g, b, a]))
    cmap = co.ListedColormap(colors, name="{bottom}2{top}".format(**locals()))
    return cmap


class AtlasError(ValueError):
    pass


# KLUDGE

import warnings

warnings.catch_warnings()
warnings.simplefilter("ignore")
