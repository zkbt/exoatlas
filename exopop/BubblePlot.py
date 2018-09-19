# general class for plotting exoplanet populations,
# painting bubbles with x + y + area onto plots

import matplotlib.pyplot as plt, numpy as np
import matplotlib.animation as animation
from craftroom.Talker import Talker

# set the aspect ratios
aspect = 768/1024.0
scale=1000
size=7

class BubblePlot(Talker):
    title = ''
    xlabel = '?'
    ylabel = '?'
    xscale = 'linear'
    yscale = 'linear'
    xlim = [None, None]
    ylim = [None, None]
    def __init__(self, pops=None):
        self.scattered = {}
        self.labels = {}

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

    def color(self):
        return None

    def label(self, key):
        return self.title.replace(' ','') + '_' + key.title()

    def build(self, interactive=False, output=False,  **kw):
        #plt.cla()
        for key in self.pops.keys():
            self.plot(key, **kw)
            if output:
                plt.savefig(self.label(key)+ '.pdf')
            if interactive:
                self.input(key)

    def kw(self, key, **kwargs):
        pop = self.pops[key]
        self.set(key)
        default = dict(s=self.size()*scale, c=self.color(), cmap='plasma', marker='o', linewidth=1, edgecolors=pop.color, alpha=pop.alpha, zorder=pop.zorder, label=pop.label)
        if default['c'] is None:
            default['facecolors']='none'
        for k, v in kwargs.items():
            default[k] = v
        return default

    def plot(self, key, labels=False, ax=None, withmass=None, custom=False, labelkw={}, **kwargs):
        self.set(key)
        try:
            self.ax
        except AttributeError:
            self.setup(ax=ax)

        '''
        if (withmass is not None)&(key.lower()=='confirmed'):
            if len(self.x) > 1:
                ok = withmass# == True
                kw = self.kw(key,**kwargs)
                try:
                    kw['alpha'] = self.pop.alpha
                except AttributeError:
                    kw['alpha'] = 1
                kw['s'] = (self.size()*scale)[ok]
                self.ax.scatter(self.x[ok], self.y[ok],**kw)

                ok = np.arange(len(self.size()))#withmass == False
                kw['alpha'] = 0.4
                kw['color'] = 'royalblue'
                kw['s'] = (self.size()*scale)[ok]
                kw['zorder'] = -100000

                self.ax.scatter(self.x[ok], self.y[ok],**kw)
            else:
                self.scattered[key] = self.ax.scatter(self.x, self.y, **self.kw(key,**kwargs))

        else:
        '''
        self.set(key)
        self.scattered[key] = self.ax.scatter(self.x, self.y, **self.kw(key,**kwargs))

        self.ax.set_xscale(self.xscale)
        self.ax.set_yscale(self.yscale)
        self.ax.set_ylim(*self.ylim)
        self.ax.set_xlim(*self.xlim)
        self.ax.set_xlabel(self.xlabel)
        self.ax.set_ylabel(self.ylabel)
        t = [1,2,3,4,5,6,7,8,9,10,20,30]
        s = ['1','2','3','4','5',' ',' ',' ',' ','10','20','30']

        if custom:
            plt.yticks(t,s)
        #self.ax.yaxis.set_major_formatter(plt.matplotlib.ticker.ScalarFormatter('{}'))
        #self.ax.yaxis.set_minor_formatter(plt.matplotlib.ticker.ScalarFormatter())

        if self.pop.labelplanets:
            self.labelplanets(**labelkw)

    def labelplanets(self, before='', after='', restrictlimits=False, **kw):
        '''
        Label the planets in a particular population.
        '''

        plt.sca(self.ax)
        for i in range(len(self.x)):
                x, y, name = self.x[i], self.y[i], self.pop.name[i]
                #name = name.replace('TRAPPIST-', 'T').replace('T1e', 'Trappist-1e')
                if restrictlimits:
                    if x < np.min(self.xlim) or x > np.max(self.xlim):
                        if y < np.min(self.ylim) or y > np.max(self.ylim):
                            continue

                textkw = dict(ha='center', va='top', fontsize=6, color=self.pop.color, alpha=self.pop.alpha)
                # think this is just as Python 3 thing
                textkw.update(**kw)

                # store the text plot, so it can be modified
                self.labels[name] = plt.text(x, y, before + name + after, **textkw)
                #self.labels[name] = plt.text(x, y, name, **textkw)
