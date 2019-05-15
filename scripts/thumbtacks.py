'''
This script makes a zoom-out animation of
RA vs distance, comparing Kepler, TESS, and
other known transiting exoplanets.

It still has some kludges in it.
''''

from exopop.ThumbtackPlot import ThumbtackPlot
from exopop.Confirmed import Kepler, TESS, Others
import numpy as np, matplotlib.pyplot as plt



# make versions with different planets labeled out to different radii (in pc)
for label in [0, 20, 30, 50, np.inf]:

	# construct a dictionary of populations of confirmed planets
	pops = dict(kepler=Kepler(),
				tess=TESS(),
				others=Others())

	# create a thumbtack plot (RA vs distance, with scaling for uniform density)
	t = ThumbtackPlot(pops=pops, lightyears=False)

	# define some plotting styles
	pops['kepler'].label='Kepler'
	pops['tess'].label='TESS'

	pops['kepler'].color = 'royalblue'
	pops['tess'].color = 'crimson'
	pops['others'].color = 'black'

	pops['kepler'].zorder=1
	pops['tess'].zorder=999
	pops['others'].zorder=0

	pops['kepler'].alpha=0.3
	pops['tess'].alpha = 0.5
	pops['others'].alpha=0.3

	# kludge clean up the system names (group planets in same star together)
	for k in [ 'kepler', 'others', 'tess']:
		p = pops[k]

		# group together the system names
		systemnames = np.array([str(n)[:-1] + '' for n in p.name])
		planetnames = [n[:-1] + '' for n in p.name]
		unique = np.unique(systemnames)
		for i, u in enumerate(unique):
			match = (systemnames == u).nonzero()[0]
			planets = ''
			for m in match:
				#print('---', u, planets)
				planets += str(p.name[m][-1])
				planetnames[m] = ''
			first = int(match[0])
			planetnames[first] = str('' + str(u) + planets)
			assert(planetnames[first] == str(u) + planets)
			#print(planetnames[first], str(u) + planets, first)

		# redefine the planet names to have only one (grouped one) per system
		pops[k].standard['name'] = planetnames
		pops[k].propagate()


	# specifiy the radius out to which planets should be labeled
	t.maxlabeldistance = label

	# create a movie
	t.movie(keys=['kepler', 'others', 'tess'], fileprefix='thumbtacks_{}pc'.format(label), maxdistance=1250.0, step=1.02)
