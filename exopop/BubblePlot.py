# general class for plotting exoplanet populations,
# painting bubbles with x + y + area onto plots

import matplotlib.pyplot as plt, numpy as np
import matplotlib.animation as animation
from zachopy.Talker import Talker

# set the aspect ratios
aspect = 768/1024.0
scale=1000
size=7

class BubblePlot(Talker):
    def __init__(self, pops=None):
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
            plt.sca(self.ax)

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
        #plt.cla()
        for key in self.pops.keys():
            self.plot(key, **kw)
            #plt.savefig(self.label(key)+ '.pdf')
            if interactive:
                self.input(key)

    def kw(self, key, **kwargs):
        pop = self.pops[key]

        default = dict(s=self.size*scale, marker='o', linewidth=1, facecolors='none', edgecolors=pop.color, alpha=0.5, zorder=pop.zorder, label=pop.label)
        for k, v in kwargs.iteritems():
            default[k] = v
        return default

    def plot(self, key, labels=False, ax=None, withmass=None, **kwargs):
        self.set(key)
        try:
            self.ax
        except AttributeError:
            self.setup(ax=ax)

        if (withmass is not None)&(key.lower()=='confirmed'):
            if len(self.x) > 1:
                ok = withmass# == True
                kw = self.kw(key,**kwargs)
                kw['alpha'] = 1
                kw['s'] = (self.size*scale)[ok]
                self.ax.scatter(self.x[ok], self.y[ok],**kw)

                ok = np.arange(len(self.size))#withmass == False
                kw['alpha'] = 0.4
                kw['color'] = 'royalblue'
                kw['s'] = (self.size*scale)[ok]
                kw['zorder'] = -100000

                self.ax.scatter(self.x[ok], self.y[ok],**kw)
            else:
                self.ax.scatter(self.x, self.y, **self.kw(key,**kwargs))

        else:
            self.ax.scatter(self.x, self.y, **self.kw(key,**kwargs))

        self.ax.set_xscale(self.xscale)
        self.ax.set_yscale(self.yscale)
        self.ax.set_ylim(*self.ylim)
        self.ax.set_xlim(*self.xlim)
        self.ax.set_xlabel(self.xlabel)
        self.ax.set_ylabel(self.ylabel)
        t = [1,2,3,4,5,6,7,8,9,10,20,30]
        s = ['1','2','3','4','5',' ',' ',' ',' ','10','20','30']

        plt.yticks(t,s)
        #self.ax.yaxis.set_major_formatter(plt.matplotlib.ticker.ScalarFormatter('{}'))
        #self.ax.yaxis.set_minor_formatter(plt.matplotlib.ticker.ScalarFormatter())

        if labels:
            best = np.argsort(self.size)
            for i in best[-10:]:
                try:
                    plt.text(self.x[i], self.y[i], self.pop.name[i])
                except:
                    print '!#!!@$!@$!@$'
                    print self.pop.name[i]


        plt.draw()
