'''
These are some common modules to import,
local to the tiny whatsup part of exoatlas.
'''

"""# zachopy contains lots of odds and ends
from craftroom.Talker import Talker
import thistothat
import exoatlas as ea

import os
def mkdir(path):
	'''A mkdir that doesn't complain if it fails.'''
	try:
		os.mkdir(path)
	except:
		pass


import numpy as np, matplotlib.pyplot as plt
import matplotlib.dates as mdates
import astropy.units, astropy.time, astropy.coordinates, astropy.io.ascii, astropy.table
from astropy.coordinates import SkyCoord

from astropy.time import Time
import astropy.units as u

from tqdm import tqdm"""

from ..imports import *

def clean(s):
    '''strip messy characters from a string'''
    new = s + ''
    for k in ' -!@#$%^&*':
        new = new.replace(k, '')
    return new
