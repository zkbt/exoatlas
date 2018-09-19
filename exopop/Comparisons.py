# various kinds of bubble plot comparisons
# (probably needs a bit of work)

import matplotlib.pyplot as plt, numpy as np
import matplotlib.animation as animation
from craftroom.Talker import Talker

from .setup import *

# set up the grid of circles to include on the plots
circlegrid = [3,10,30,100,300,1000]
circlegrid = [10,30,50]

def test(t):
    t.plot('tess')
    t.zoom(100)
    t.highlight((t.pop.planet_radius > 2)*(t.pop.planet_radius < 4))
    plt.draw()

# set the size scale of the points
scale=1000
size=7
angle = 67.5*np.pi/180


# set the aspect ratios
aspect = 768/1024.0

# decide how to stretch radial coordinates
def stretch(d):
    return np.array(d)**1.5

plt.ion()
class BubblePlot(Talker):
    def __init__(self):
        Talker.__init__(self)
        self.named = []
        self.pops = pops

    def set(self, key):
        self.pop = self.pops[key]
        self.key = key

    def setup(self, ax=None):
        self.aspect = aspect

        if ax is None:
            plt.figure(self.title, figsize=(size,size*self.aspect), dpi=100)
            gs = plt.matplotlib.gridspec.GridSpec(1,1,wspace=0,hspace=0)
            self.ax = plt.subplot(gs[0])
        else:
            self.ax = ax

        #self.ax.set_title(self.title)
        #self.ax_scale = plt.subplot(gs[0,1])
        #self.ax_colors = plt.subplot(gs[1,1])
        #self.ax_scale.set_visible(True)
        #self.ax_colors.set_visible(False)

        #n = 5
        #self.ax_scale.scatter(np.zeros(n), np.logspace(0,np.log10(scale),n)[::-1], np.logspace(0,np.log10(scale),n))
        #self.ax_scale.text('rotate=90)

    def label(self, key):
        return self.title.replace(' ','') + '_' + key.title()

    def build(self, interactive=False, **kw):
        plt.cla()
        for key in self.pops.keys():
            self.plot(key, **kw)
            plt.savefig(self.label(key)+ '.pdf')
            if interactive:
                self.input(key)

    def kw(self, key):
        return dict(s=self.size()*scale, marker='o', linewidth=2, facecolors='none', edgecolors=colors[key], alpha=0.5, zorder=zorders[key], label=names[key])

    def plot(self, key, labels=False, ax=None, **kw):
        self.set(key)
        try:
            self.ax
        except:
            self.setup(ax=ax)
        self.ax.scatter(self.x, self.y, **self.kw(key))
        self.ax.set_ylim(*self.ylim)
        self.ax.set_xlim(*self.xlim)
        self.ax.set_xscale(self.xscale)
        self.ax.set_yscale(self.yscale)
        self.ax.set_xlabel(self.xlabel)
        self.ax.set_ylabel(self.ylabel)
        t = [1,2,3,4,5,6,7,8,9,10,20,30]
        s = ['1','2','3','4','5',' ',' ',' ',' ','10','20','30']

        plt.yticks(t,s)
        #self.ax.yaxis.set_major_formatter(plt.matplotlib.ticker.ScalarFormatter('{}'))
        #self.ax.yaxis.set_minor_formatter(plt.matplotlib.ticker.ScalarFormatter())

        if labels:
            best = np.argsort(self.size())
            for i in best[-10:]:
                try:
                    plt.text(self.x[i], self.y[i], self.pop.name[i])
                except:
                    print('!#!!@$!@$!@$')
                    print(self.pop.name[i])


        plt.draw()



class thumbtacks(BubblePlot):
    def __init__(self):
        BubblePlot.__init__(self)
        self.title = 'Thumbtacks'
        self.xlabel = ''
        self.ylabel = ''
        self.ylim=[0.8, 30.0]
        self.xlim=[200,3000]
        self.xscale='linear'
        self.yscale='log'

    def toclock(self, hour):
        return (6 - hour)*2*np.pi/24.0

    def setup(self):

        # set the figure's aspect ratio
        self.aspect = aspect

        # maek the figure
        plt.figure(self.title, figsize=(size,size*self.aspect), dpi=100)

        # populate with an axis
        gs = plt.matplotlib.gridspec.GridSpec(1,1,left=0,right=1,bottom=0,top=1,hspace=0,wspace=0)
        self.ax = plt.subplot(gs[:,:])
        self.ax.set_aspect('equal')

        # make the lines go away
        #self.ax.set_frame_on(False)
        self.ax.get_xaxis().set_visible(False)
        self.ax.get_yaxis().set_visible(False)

        # add RA labels
        hours = np.arange(0, 24, 3)
        theta = self.toclock(hours)
        months = ['September', 'December', 'March', 'June']
        clockr = 0.925
        for i in range(len(hours)):
            x =clockr*(np.cos(theta[i])*0.5)*self.aspect+0.5
            y = clockr*np.sin(theta[i])*0.5 + 0.5
            self.ax.text(x, y, '{0}h'.format(hours[i]), rotation=-90+theta[i]*180/np.pi, transform=self.ax.transAxes, va='center', ha='center', size=10, alpha=0.5)
            if (i % 2) == 0:
                self.ax.text(x, y, '\n\n' + months[i/2], size=6,  rotation=-90+theta[i]*180/np.pi, transform=self.ax.transAxes, va='center', ha='center', alpha=0.35)

    def circles(self, radii):


        ax = self.ax


        theta = np.linspace(0,2*np.pi,1000)
        gridkw = dict( color='gray', alpha=0.8, zorder=-100)
        self.circlelabels = {}
        nudge=1.15
        for originalradius in radii:
            radii = stretch(originalradius)
            nudge = 0.0
            ax.plot(radii*np.cos(theta), radii*np.sin(theta), linewidth=4, **gridkw)
            self.circlelabels[originalradius] = ax.text((nudge+radii)*np.cos(angle), (nudge+radii)*np.sin(angle), '{0:.0f} pc'.format(originalradius), rotation=-90+ angle*180/np.pi, va='center', ha='center', size=13, weight='extra bold', **gridkw)

    @property
    def size(self):
        return 50.0/1000.0

    def build(self, interactive=False):
        plt.cla()
        for key in ['new', 'confirmed', 'tess']:
            self.plot(key)
            for z in [10,30,100,300,1000]:
                self.zoom(z)
                self.clearnames()
                self.namestars('known')
                if key == 'kepler':
                    self.namestars('kepler')
                plt.draw()
                plt.savefig(self.label(key)+ '.pdf')

            if interactive:
                self.input(key)


    def movie(self, step=1.02, interactive=False, bitrate=10000, highlight=''):
        metadata = dict(title='Exoplanets Zoom', artist='Zach Berta-Thompson (zkbt@mit.edu)')
        self.writer = animation.FFMpegWriter(fps=15, metadata=metadata, bitrate=bitrate)

        plt.cla()
        for key in ['new', 'kepler', 'tess']:
            self.plot(key)
            if highlight == 'completeness':
                self.highlight((self.pop.planet_radius > 2)*(self.pop.planet_radius < 4)*(self.pop.period < 10), '2-4 Earth radii\nP<10 days')
            if highlight == 'small':
                self.highlight((self.pop.planet_radius < 2), '<2 Earth radii')
            if highlight == 'habitable':
                self.highlight((self.pop.planet_radius > 0.7)*(self.pop.planet_radius < 1.6)*(self.pop.teq < 310)*(self.pop.teq > 200), 'Potentially Habitable Planets')

            f = plt.gcf()
            filename = 'exoplanets_zoom_{0}{1}.mp4'.format(key, highlight)
            self.speak('writing movie to {0}'.format(filename))
            z = 2.0

            with self.writer.saving(f, filename, 1024.0/size):
                # loop over exposures
                while(z < 1100):
                    self.zoom(z)

                    self.clearnames()
                    if highlight != 'habitable':
                        self.namestars('known')
                        if key == 'kepler':
                            self.namestars('kepler')

                    self.writer.grab_frame()
                    self.speak('zoomed to {0}'.format(z))
                    z *= step
                    #plt.draw()
                #plt.savefig(self.label(key))


    def zoom(self, distance):
        self.outer = distance
        height = stretch(self.outer)*1.2
        width = height/self.aspect
        self.ax.set_ylim(-height, height)
        self.ax.set_xlim(-width, width)

        nudge = 0.05*stretch(self.outer)
        for d in self.circlelabels.keys():
            self.circlelabels[d].set_visible(d > 0.2*distance)
            self.circlelabels[d].set_x((nudge+stretch(d))*np.cos(angle))
            self.circlelabels[d].set_y((nudge+stretch(d))*np.sin(angle))


    def label(self, key):
        return self.title.replace(' ','') + '_' + key.title()  + '_{0:04.0f}pc'.format(self.outer)+'.pdf'

    def namestars(self, key):

        if True:
                old = self.key + ''
                self.set(key)

                xlim = self.ax.get_xlim()
                ylim = self.ax.get_ylim()
                onplot = (self.x > np.min(xlim))*(self.x < np.max(xlim))* (self.y > np.min(ylim))*(self.y < np.max(ylim))
                nottooclose = self.pop.distance > self.outer*0.2
                nottoofar = self.pop.distance <= 30
                #print np.sum(onplot*nottooclose*nottoofar)
                tolabel = (nottooclose*onplot*nottoofar).nonzero()[0]
                #tolabel = self.pop.find('WASP94Ab')
                if tolabel.size() > 1:
                    tolabel = tolabel[np.unique(self.x[tolabel], return_index=True)[1]]

                for c in tolabel:
                    #print self.pop.name[c]
                    self.named.append(plt.text(self.x[c], self.y[c] + stretch(self.outer)*0.02,  self.pop.name[c].replace('Kepler-444 b', 'Kepler-444 bcdef') , color=colors[key], alpha=0.75, va='bottom', ha='center', weight='bold', size=10))
                self.set(old)

    def clearnames(self):

        while(len(self.named) > 0):
            self.named.pop().remove()

    def plot(self, key, labels=False):
        self.set(key)
        try:
            self.ax
        except:
            self.setup()
            self.circles(circlegrid)
        kw = self.kw(key)
        kw['facecolors'] = colors[key]
        kw['edgecolors'] = 'none'
        try:
            self.signature
        except:
            self.signature = self.ax.text(0.02, 0.02, 'zoom by Zach Berta-Thompson, early 2015', transform=self.ax.transAxes, alpha=0.5, size=8)
        self.ax.scatter(self.x, self.y, **kw)
        self.ax.scatter(0,0, marker='x', s=100, alpha=1, color='lightgray', linewidth=4, zorder=-90)

        self.leg = plt.legend(bbox_to_anchor=(1, 1), fontsize=10, framealpha=0.0, scatterpoints=1, markerscale=1.5, title='Transiting Exoplanets')

    def highlight(self, indices, label='special!'):
        print(self.pop.standard[indices])
        kw = self.kw(self.key)
        kw['marker'] = '*'
        kw['edgecolors'] = colors[self.key].replace('black', 'gray')
        kw['facecolors'] = 'white'#colors[self.key]
        kw['s'] *= 8
        kw['alpha'] = 1
        kw['label'] = None
        kw['zorder'] = 20000
        handle = self.ax.scatter(self.x[indices], self.y[indices], **kw)
        #self.highlightleg = plt.legend(label, fontsize=10, framealpha=0, markerscale=1.5)
        self.speak('highlighting {0} points'.format(len(self.x[indices])))
        t = self.pop.teq[indices]
        r = self.pop.planet_radius[indices]
        #kic = self.pop.standard['kepid'][indices]
        for i in range(len(t)):
            self.ax.text(self.x[indices][i], self.y[indices][i], r'{0:.0f}K, {1:.1f}R$_\oplus$'.format(t[i], r[i]) +'\n', size=8, ha='center', va='bottom', color=colors[self.key].replace('black', 'gray'))
            print(r[i], t[i])#, self.pop.distance[indices][i]
        #(self.name == 'WASP-94A b').nonzero()[0]

    @property
    def r(self):
        return stretch(self.pop.distance)

    @property
    def theta(self):
        return self.toclock(self.pop.ra*24.0/360.0)

    @property
    def x(self):
        return stretch(self.pop.distance)*np.cos(self.theta)

    @property
    def y(self):
        return stretch(self.pop.distance)*np.sin(self.theta)




class trans(BubblePlot):
    def __init__(self):
        BubblePlot.__init__(self)
        self.title = 'Transmission Spectroscopy'
        self.xlabel = 'Planet Equilibrium Temperature (K)'
        self.ylabel = 'Planet Radius (Earth radii)'
        self.ylim=[0.8, 30.0]
        self.xlim=[200,3000]
        self.xscale='linear'
        self.yscale='log'

        self.set('known')
        self.normalization = self.unnormalizedsize()[self.pop.find('GJ1214b')]

    def unnormalizedsize(self):
        return (self.pop.transmissionsignal/self.pop.noisepertime)**2

    @property
    def size(self):
        return self.unnormalizedsize()/self.normalization

    @property
    def x(self):
        return self.pop.teq

    @property
    def y(self):
        return self.pop.planet_radius

class empty(trans):
    def __init__(self):
        trans.__init__(self)
        self.title = 'Empty'

    def build(self):
        self.setup()
        self.plot('known')
        plt.savefig(self.label('known') + '.pdf')

    def plot(self, key, labels=False):
        self.set(key)
        try:
            self.ax
        except:
            self.setup()
        #self.ax.scatter(self.x, self.y, **self.kw(key))
        self.ax.set_ylim(*self.ylim)
        self.ax.set_xlim(*self.xlim)
        self.ax.set_xscale(self.xscale)
        self.ax.set_yscale(self.yscale)
        self.ax.set_xlabel(self.xlabel)
        self.ax.set_ylabel(self.ylabel)
        t = [1,2,3,4,5,6,7,8,9,10,20,30]
        s = ['1','2','3','4','5',' ',' ',' ',' ','10','20','30']

        plt.yticks(t,s)
        plt.draw()



class summ(trans):
    def __init__(self):
        trans.__init__(self)
        self.title = 'Summary'

    def kw(self, key):
        return dict(s=self.size()*scale, marker='o', linewidth=2, edgecolors='none', facecolors=colors[key], alpha=0.5, zorder=zorders[key], label=names[key])


    @property
    def size(self):
        return 40.0/1000.0

class emis(trans):
    def __init__(self):
        trans.__init__(self)
        self.title = 'Emission Spectroscopy'

        self.set('known')
        self.normalization = self.unnormalizedsize()[self.pop.find('WASP43b')]

    def unnormalizedsize(self):
        return (self.pop.emissionsignal/self.pop.noisepertime)**2


class dept(emis):
    def __init__(self):
        trans.__init__(self)
        self.title = 'Transit Depth'

        self.set('known')
        self.normalization = self.unnormalizedsize()[self.pop.find('HD209458b')]

    def unnormalizedsize(self):
        return (self.pop.depth/self.pop.noisepertime)**2


class refl(emis):
    def __init__(self):
        trans.__init__(self)
        self.title = 'Reflected Light'

        self.set('known')
        self.normalization = self.unnormalizedsize()[self.pop.find('HD209458b')]

    def unnormalizedsize(self):
        return (self.pop.depth/self.pop.a_over_r**2/self.pop.noisepertime)**2
