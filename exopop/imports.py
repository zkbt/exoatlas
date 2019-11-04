# imports that are need by many exopop subsections
import os, sys, time, shutil, warnings
from tqdm import tqdm

import numpy as np, matplotlib.pyplot as plt, matplotlib.animation as animation
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter, LogLocator



from astropy.utils.exceptions import AstropyDeprecationWarning
warnings.simplefilter('ignore', category=AstropyDeprecationWarning)

# this function downloads a file and returns its filepath
from astropy.utils.data import download_file
from astropy.io import ascii
from astropy.table import Table, vstack, join
from astropy.visualization import quantity_support
quantity_support()

# some general custom utilities from Zach
import craftroom.strings,  craftroom.utils, craftroom.color
from .talker import Talker

# units and constants from astropy
import astropy.units as u, astropy.constants as con
from astropy.time import Time

# create a directory structure in the user's home directory
base = os.path.join(os.getenv('HOME'), '.exopop')
craftroom.utils.mkdir(base)
directories = dict(
                    data=os.path.join(base, 'data/'),
                    plots=os.path.join(base, 'plots/'),
                    models=os.path.join(base, 'models/')
                    )
for k in directories.keys():
    craftroom.utils.mkdir(directories[k])

# some kludge for dealing with Python 3 vs 2?
try:
    FileNotFoundError
except NameError:
    FileNotFoundError = IOError

def clean(s):
    '''
    A wrapper function to clean up complicated strings.
    '''
    bad = ' !@#$%^&*()-,./<>?'
    cleaned = s + ''
    for c in bad:
        cleaned = cleaned.replace(c, '')
    return cleaned

def time_from_modified(filename):
    '''
    How long ago was this file last modified?
    '''
    try:
        dt = Time.now().unix - os.path.getmtime(filename)
        return dt/60/60/24
    except FileNotFoundError:
        return np.inf


def check_if_needs_updating(filename, maximum_age=1.0):
    '''
    Do a (possibly interactive) check to see if this file
    is so old that it needs to be updated.

    Returns
    -------
    old : bool
        True if the file is so old it needs updating
        False if the file doesn't need updating
    '''

    # how long ago was the data updated?
    dt = time_from_modified(filename)

    old = False
    if dt == np.inf:
        old = True
    elif dt > maximum_age:
        print(f'{filename} is {dt:.3f} days old.')
        old = 'y' in input('Should it be updated? [y/N]').lower()

    return old
