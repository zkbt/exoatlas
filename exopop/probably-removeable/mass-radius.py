class MassRadius(BubblePanel):
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
        self.point_at(key)
        try:
            self.ax
        except AttributeError:
            self.setup(ax=ax)

        plt.sca(self.ax)
        self.ax.set_xscale(self.xscale)
        self.ax.set_yscale(self.yscale)

        rlower = np.abs(self.pop.planet_radius_uncertainty_lower)
        rupper = np.abs(self.pop.planet_radius_uncertainty_upper)
        rissmallenough = 0.5*(rlower+rupper)/self.pop.planet_radius < 1.0/1.0#2.5#/2.5
        risntzero = 0.5*(rlower+rupper) > 0.0

        mlower = np.abs( self.pop.planet_mass_uncertainty_lower)
        mupper = np.abs(self.pop.planet_mass_uncertainty_upper)
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
            xerr = np.vstack([-self.pop.planet_mass_uncertainty_lower[self.ok].filled(), self.pop.planet_mass_uncertainty_upper[self.ok].filled()])
        except:
            xerr = np.vstack([-self.pop.planet_mass_uncertainty_lower[self.ok], self.pop.planet_mass_uncertainty_upper[self.ok]])

        y = self.pop.planet_radius[self.ok]
        try:
            yerr = np.vstack([-self.pop.planet_radius_uncertainty_lower[self.ok].filled(), self.pop.planet_radius_uncertainty_upper[self.ok].filled()])
        except:
            yerr = np.vstack([-self.pop.planet_radius_uncertainty_lower[self.ok], self.pop.planet_radius_uncertainty_upper[self.ok]])


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
            r, g, b = name2color(color.lower())
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
