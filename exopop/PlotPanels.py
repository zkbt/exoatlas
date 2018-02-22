from exopop.imports import *
from exopop import BubblePlot
from zachopy.painting import ink_errorbar
import zachopy.units as u

ylim = [0.3, 5.0]




class DistanceRadius(BubblePlot):
    def __init__(self, lightyears=False, **kw):
        BubblePlot.__init__(self, **kw)
        self.title = 'Context'
        self.lightyears = lightyears
        if lightyears:
            self.xlabel = "Distance from Earth (lightyears)"
            xlim = np.array([4,1000])*u.pc/u.ly
        else:
            self.xlabel = 'Distance\n(parsecs)'
            xlim = np.array([4,1000])

        self.ylabel = 'Planet Radius (Earth radii)'
        self.ylim=ylim
        self.xlim=xlim
        self.xscale='log'
        self.yscale='linear'
        self.set('goodmass')
        self.normalization = self.unnormalizedsize[self.pop.standard['name'] == 'GJ 1214b']

    @property
    def unnormalizedsize(self):
        # set the symbol size to be transit depth
        return (self.pop.planet_radius/self.pop.stellar_radius)**2

    @property
    def size(self):
        return 0.5*self.unnormalizedsize/self.normalization

    @property
    def x(self):
        if self.lightyears:
            return self.pop.distance*u.pc/u.ly
        else:
            return self.pop.distance

    @property
    def y(self):
        return self.pop.planet_radius

    def build(self, interactive=False, output=False, **kw):
        #plt.cla()
        for key in self.pops.keys():
            if key == 'goodmass':
                kw['alpha'] = 1
            elif key == 'badmass':
                kw['alpha'] = 0.5
            elif key == 'unconfirmed':
                kw['alpha'] = 0.25
            elif key == 'tess':
                kw['alpha'] = 0.5
            elif key == 'kepler':
                kw['alpha'] = 0.5
            else:
                kw['alpha'] = 1
            self.plot(key, **kw)

class Teq(DistanceRadius):
    def __init__(self, **kw):
        DistanceRadius.__init__(self, **kw)
        self.title = 'More Context'
        self.xlabel = 'Insolation\n(relative to Earth)'
        self.xlim=[5e4, 0.5]
        self.xscale='log'
        self.set('goodmass')
        self.normalization = self.unnormalizedsize[self.pop.standard['name'] == 'GJ 1214b']

    @property
    def x(self):
        a_over_rearth = u.au/u.Rsun
        teqearth = u.Tsun/np.sqrt(2*u.au/u.Rsun)
        return (self.pop.teq/teqearth)**4

class PlanetDensity(DistanceRadius):
    def __init__(self, **kw):
        DistanceRadius.__init__(self, **kw)
        self.title = 'More Context'
        self.xlabel = 'Planet Density\n(g/cc)'
        self.xlim=[0.3, 12]
        self.xscale='log'
        self.set('goodmass')
        self.normalization = self.unnormalizedsize[self.pop.standard['name'] == 'GJ 1214b']

    @property
    def x(self):
        mass = self.pop.planet_mass*u.Mearth
        volume =4*np.pi*(self.pop.planet_radius*u.Rearth)**3/3.0
        return mass/volume


class StellarRadius(DistanceRadius):
    def __init__(self, **kw):
        DistanceRadius.__init__(self, **kw)
        self.title = 'More Context'
        self.xlabel = 'Stellar Radius\n(solar radii)'
        self.xlim=[0.1, 1.1]
        self.xscale='linear'
        self.set('goodmass')
        self.normalization = self.unnormalizedsize[self.pop.standard['name'] == 'GJ 1214b']

    @property
    def x(self):
        return self.pop.stellar_radius



class MassRadius(BubblePlot):
    def __init__(self, **kw):
        BubblePlot.__init__(self, **kw)
        self.title = 'Context'
        self.xlabel = 'Planet Mass\n(Earth masses)'
        self.ylabel = 'Planet Radius (Earth radii)'
        self.ylim=ylim
        self.xlim=[0.7, 30]
        self.xscale='log'
        self.yscale='linear'
        self.set('goodmass')
        #self.pop.find('GJ1214b')
        self.normalization = self.unnormalizedsize[self.pop.standard['name'] == 'GJ 1214b']

    @property
    def unnormalizedsize(self):
        return (1.0/self.pop.stellar_radius)**2

    @property
    def size(self):
        return 0.25*self.unnormalizedsize/self.normalization

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
        kw = dict(marker='o', linewidth=0, elinewidth=width, alpha=1.0,  label=self.pop.name, color=self.pop.color, markeredgecolor=self.pop.color, capthick=width, capsize=2, markersize=3)

        #    self.ax.errorbar(x, y, xerr=xerr, yerr=yerr, **kw)

        color = self.pop.color
        r, g, b = zachopy.color.name2color(color.lower())
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

        print("!$!@%*!#%!@$!@#")
        print(self.pop)
        print(weights)
        #rgba[:,3] = 1.0*weights
        rgba[:,:3] = 1.0 - (1.0 - rgba[:,:3])*weights.reshape(len(weights), 1)
        rgba[:,3] = 1.0#*weights
        #print rgba
        if (len(x) > 1)&(self.pop.ink):
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
            best = np.argsort(self.size)
            for i in best[-10:]:
                try:
                    plt.text(self.x[i], self.y[i], self.pop.name[i])
                except:
                    print('!#!!@$!@$!@$')
                    print(self.pop.name[i])


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
