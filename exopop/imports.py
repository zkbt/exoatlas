# imports that are need by many exopop subsections

from zachopy.Talker import Talker
import numpy as np, matplotlib.pyplot as plt
import astropy.io.ascii, astropy.table
import zachopy.relations, zachopy.strings, zachopy.units, zachopy.utils

directories = dict(data='data/', plots='plots/', models='models/')
for k in directories.keys():
    zachopy.utils.mkdir(directories[k])
