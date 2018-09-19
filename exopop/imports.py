# imports that are need by many exopop subsections

from craftroom.Talker import Talker
import os, sys
import numpy as np, matplotlib.pyplot as plt, matplotlib.animation as animation
import astropy.io.ascii, astropy.table

if sys.version_info[0] == 2:
    import urllib as request
else:
    from urllib import request

import thistothat, craftroom.strings, craftroom.units, craftroom.utils, craftroom.color
import astropy.units as u, astropy.constants as con
base = os.path.join(os.getenv('HOME'), '.exopop')
craftroom.utils.mkdir(base)
directories = dict(
                    data=os.path.join(base, 'data/'),
                    plots=os.path.join(base, 'plots/'),
                    models=os.path.join(base, 'models/')
                    )

for k in directories.keys():
    craftroom.utils.mkdir(directories[k])

try:
    FileNotFoundError
except NameError:
    FileNotFoundError = IOError
