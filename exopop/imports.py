# imports that are need by many exopop subsections

from zachopy.Talker import Talker
import numpy as np, matplotlib.pyplot as plt, matplotlib.animation as animation
import astropy.io.ascii, astropy.table
import urllib
import zachopy.relations, zachopy.strings, zachopy.units, zachopy.utils, zachopy.color

directories = dict(data='data/', plots='plots/', models='models/')
for k in directories.keys():
    zachopy.utils.mkdir(directories[k])
