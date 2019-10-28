from ..imports import *
from .BubblePlot import BubblePlot
from craftroom.painting import ink_errorbar
import craftroom.units as u

class DistanceRadius(BubblePlot):
    title = ''
    xlabel = 'Distance\n(parsecs)'
    ylabel = 'Planet Radius (Earth radii)'
    xscale = 'log'
    yscale = 'log'
    xlim = [4,1000]
    ylim = [0.4, 40]

    def __init__(self, lightyears=False, **kw):
        BubblePlot.__init__(self, **kw)

        self.lightyears = lightyears
        if lightyears:
            self.xlabel = "Distance from Earth (lightyears)"
            self.xlim = np.array([4,1000])*(u.pc/u.ly)

        self.normalization = 0.02#self.unnormalizedsize()[self.pop.standard['name'] == 'GJ 1214b']

    def unnormalizedsize(self):
        # set the symbol size to be transit depth
        #return (self.pop.planet_radius/self.pop.stellar_radius)**2
        return self.pop.depth

    def size(self):
        return 0.5*self.unnormalizedsize()/self.normalization

    @property
    def x(self):
        if self.lightyears:
            return self.pop.distance*(u.pc/u.ly)
        else:
            return self.pop.distance

    @property
    def y(self):
        return self.pop.planet_radius

    def build(self, interactive=False, output=False, **kw):
        #plt.cla()
        for key in self.pops.keys():
            self.plot(key, **kw)

class EscapeRadius(DistanceRadius):
    xlabel = '$\lambda = E_{grav}/E_{thermal}$'
    xscale = 'log'
    xlim = [None, None]

    @property
    def x(self):
        return self.pop.escape_parameter

class FluxRadius(DistanceRadius):
    xlabel = 'Flux Received\n(relative to Earth)'
    xscale = 'log'
    xlim = [5e4, 0.5]

    @property
    def x(self):
        a_over_rearth = u.au/u.Rsun
        teqearth = u.Tsun/np.sqrt(2*u.au/u.Rsun)
        return (self.pop.teq/teqearth)**4
Teq = FluxRadius

class PlanetDensityRadius(DistanceRadius):

    title = ''
    xlabel = 'Planet Density\n(g/cm$^3$)'
    xscale = 'log'
    xlim = [0.3, 12]

    @property
    def x(self):
        mass = self.pop.planet_mass*u.Mearth
        volume =4*np.pi*(self.pop.planet_radius*u.Rearth)**3/3.0
        return mass/volume

class StellarRadius(DistanceRadius):
    xlabel = 'Stellar Radius\n(solar radii)'
    xscale = 'linear'
    xlim = [0.1, 1.1]

    @property
    def x(self):
        return self.pop.stellar_radius

class JRadius(DistanceRadius):
    xlabel = 'J (magnitude)\n'
    xscale = 'linear'
    xlim = [3.5, 14.5]

    @property
    def x(self):
        return self.pop.J

class PeriodRadius(DistanceRadius):
    xlabel = 'Period (days)\n'
    xscale = 'log'
    xlim = [0.15, 365]

    @property
    def x(self):
        return self.pop.period

class MassRadius(BubblePlot):
    title = ''
    xlabel = 'Planet Mass\n(Earth masses)'
    ylabel = 'Planet Radius (Earth radii)'
    xscale = 'log'
    yscale = 'log'

    xlim = [None, None]
    ylim = [None, None]

    @property
    def x(self):
        return self.pop.planet_mass

    @property
    def y(self):
        return self.pop.planet_radius

    def plot(self, key, labels=False, ax=None, **kw):
        self.set(key)
        try:
            self.ax
        except AttributeError:
            self.setup(ax=ax)

        plt.sca(self.ax)
        self.ax.set_xscale(self.xscale)
        self.ax.set_yscale(self.yscale)

        rlower = np.abs(self.pop.planet_radius_lower)
        rupper = np.abs(self.pop.planet_radius_upper)
        rissmallenough = 0.5*(rlower+rupper)/self.pop.planet_radius < 1.0/1.0#2.5#/2.5
        risntzero = 0.5*(rlower+rupper) > 0.0

        mlower = np.abs( self.pop.planet_mass_lower)
        mupper = np.abs(self.pop.planet_mass_upper)
        missmallenough = 0.5*(mlower+mupper)/self.pop.planet_mass < 1.0/1.0#2.5#/2.5
        misntzero = 0.5*(mlower+mupper) > 0.0

        # solar system kludge!
        exact = (risntzero == False).all() and (misntzero == False).all()
        if exact:
            self.ok = np.arange(len(self.pop.planet_radius))
        else:
            self.ok = (rissmallenough*risntzero*missmallenough*misntzero)
        #assert(len(self.ok) ==1)
        #self.ok = ok.nonzero()

        #self.ok = np.arange(len(self.pop.planet_mass))
        x = self.pop.planet_mass[self.ok]
        try:
            xerr = np.vstack([-self.pop.planet_mass_lower[self.ok].filled(), self.pop.planet_mass_upper[self.ok].filled()])
        except:
            xerr = np.vstack([-self.pop.planet_mass_lower[self.ok], self.pop.planet_mass_upper[self.ok]])

        y = self.pop.planet_radius[self.ok]
        try:
            yerr = np.vstack([-self.pop.planet_radius_lower[self.ok].filled(), self.pop.planet_radius_upper[self.ok].filled()])
        except:
            yerr = np.vstack([-self.pop.planet_radius_lower[self.ok], self.pop.planet_radius_upper[self.ok]])


        #print self.ok
        width=1
        if exact:
            plotkw = dict(label=self.pop.name, color=self.pop.color, edgecolor=self.pop.color, **kw)
            plotkw['alpha'] = 1
            plotkw['zorder'] = 1e9
            self.scattered[key] = plt.scatter(x, y, **plotkw)
        else:
            kw = dict(marker='o', linewidth=0, elinewidth=width, alpha=1.0,  label=self.pop.name, color=self.pop.color, markeredgecolor=self.pop.color, capthick=width, capsize=2, markersize=3)

            #    self.ax.errorbar(x, y, xerr=xerr, yerr=yerr, **kw)

            color = self.pop.color
            r, g, b = craftroom.color.name2color(color.lower())
            n = len(x)
            rgba = np.zeros((n, 4))
            rgba[:,0] = r
            rgba[:,1] = g
            rgba[:,2] = b

            mass = x
            radius = y
            density = mass/radius**3
            merr = 0.5*(mlower + mupper)[self.ok]
            rerr = 0.5*(rlower + rupper)[self.ok]

            #density_uncertainty = np.sqrt((merr/mass)**2 + 9*(rerr/radius)**2)
            #density = mass/radius**3
            even_uncertainty = np.sqrt((merr/mass)**2 + (rerr/radius)**2)#*density
            density_uncertainty = np.sqrt((merr/mass)**2 + 3*(rerr/radius)**2)#*density
            gravity_uncertainty = np.sqrt((merr/mass)**2 + 2*(rerr/radius)**2)

            #assert(False)
            weights = np.minimum(1.0/density_uncertainty**2, 50.0)/50
            #over = weights > 1
            #weights[over] = 1
            kw['zorder'] = weights


            #rgba[:,3] = 1.0*weights
            rgba[:,:3] = 1.0 - (1.0 - rgba[:,:3])*weights.reshape(len(weights), 1)
            rgba[:,3] = 1.0#*weights
            #print rgba

            if (len(x) > 1)&(self.pop.plotkw.get('ink', True)):
                self.speak(key)
                ink_errorbar(x, y, xerr=xerr, yerr=yerr, colors=rgba, **kw)
                #assert(False)
            else:
                kw['zorder'] = 1e6
                self.ax.errorbar(x, y, xerr=xerr, yerr=yerr, **kw)

        #l = []
        #l.extend(self.pop.name)

        #if len(self.ok) > 1:
        """
        if len(x) > 0:
            for i, name in enumerate(self.pop.name[self.ok]):
                self.ax.text(x[i], y[i], name)"""



        self.ax.set_ylim(*self.ylim)
        self.ax.set_xlim(*self.xlim)
        self.ax.set_xlabel(self.xlabel)
        self.ax.set_ylabel(self.ylabel)


        #self.ax.yaxis.set_major_formatter(plt.matplotlib.ticker.ScalarFormatter('{}'))
        #self.ax.yaxis.set_minor_formatter(plt.matplotlib.ticker.ScalarFormatter())

        if labels:
            for i in range(len(self.x)):
                #try:
                plt.text(self.x[i], self.y[i], self.pop.name[i])
                #except:
                #    print('!#!!@$!@$!@$')
                #    print(self.pop.name[i])


        plt.draw()

    def fusswithticks(self):
        plt.sca(self.ax)
        t = [1,2,3,4,5,6,7,8,9,10,20,30]
        s = ['1','2','3','4',' ',' ',' ',' ',' ','10','20']
        plt.yticks(t,s)
        t = [.7,.8,.9,1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 0.5]
        s = ['{0}'.format(v) for v in t]

        s[0] = ' '
        s[1] = ' '
        s[2] = ' '
        s[8] = ' '
        s[9] = ' '
        s[10] = ' '
        s[11] = ' '


        plt.xticks(t,s)







class MultiPanelPlot(Talker):
    '''
    Make and modify row or a column of connected exopop plot panels.
    '''


    def __init__(self, pops,
                       panels=[MassRadius, FluxRadius, StellarRadius, DistanceRadius],
                       horizontal=True,
                       figsize=(12,8),
                       gridspec_kw=dict(hspace=0.1,
                                            left=0.15,
                                            right=0.95,
                                            bottom=0.15,
                                            wspace=0.05)):
        '''
        Set up the plotting panels,
        and pause for modifcations.
        '''

        # store list of panel names
        self.names = [p.__name__ for p in panels]
        self.pops = pops

        # create a panel plotting object
        self.panels = {p.__name__:p(pops=pops) for p in panels}

        # set up the geometry
        self.horizontal = horizontal
        if horizontal:
            nrows, ncols = 1, len(self.panels)
            sharex, sharey = False, True
        else:
            nrows, ncols = len(self.panels), 1
            sharex, sharey = True, False

        # set up the plotting axes, and store in dictionary
        self.figure, axgrid = plt.subplots(nrows, ncols,
                                           sharex=sharex, sharey=sharey,
                                           figsize=figsize,
                                           gridspec_kw=gridspec_kw)
        self.ax = {}

        for i, k in enumerate(self.names):
            if nrows*ncols == 1:
                self.ax[k] = axgrid
            else:
                self.ax[k] = axgrid[i]

    def build(self, **kw):
        '''
        Actually make the plot.
        '''
        # plot each population in each panel
        for i, k in enumerate(self.names):
            self.panels[k].build(ax=self.ax[k], **kw)

        # clean up unnecessary labels
        if self.horizontal:
            for k in self.names[1:]:
                plt.setp(self.ax[k].get_yticklabels(), visible=False)
                self.ax[k].set_ylabel('')
        else:
            for k in self.names[:-1]:
                plt.setp(self.ax[k].get_xticklabels(), visible=False)
                self.ax[k].set_xlabel('')

FourPanels = MultiPanelPlot
