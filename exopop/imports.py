# imports that are need by many exopop subsections
import os, sys, time, shutil

import numpy as np, matplotlib.pyplot as plt, matplotlib.animation as animation

# this function downloads a file and returns its filepath
from astropy.utils.data import download_file
from astropy.io import ascii
from astropy.table import Table, vstack

# some general custom utilities from Zach
import thistothat, craftroom.strings, craftroom.units, craftroom.utils, craftroom.color
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
