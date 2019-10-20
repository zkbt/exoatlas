# general class for exoplanet populations
from ..imports import *
import string

# the mamajek relation is needed for estimating distances to stars, for those without
mamajek = thistothat.Mamajek()
mamajek._pithy = True

necessary_columns = [
'name',
'period',
'transit_epoch',
'transit_duration',
'stellar_teff',
'stellar_mass',
'stellar_radius',
'stellar_distance',
'planet_radius',
'planet_mass',
'a_over_r',
'b',
'ra', 'dec',
'J',
'discover']

desired_columns = [
'radius_ratio',
'planet_mass_upper',
'planet_mass_lower',
'planet_radius_upper',
'planet_radius_lower',
'stellar_distance_upper',
'stellar_distance_lower']




class PopulationFromStandard(Talker):
    '''
    Create a population from a standardized table.
    '''

    def __init__(self, standard,  **plotkw):
        '''
        Initialize a Population of exoplanets from a standardized table.

        Parameters
        ----------
        standard : astropy.table.Table
            A table that contains all the necessary columns.
        **plotkw : dict
            All other keyword arguments wil
        '''

        # a standardized table with a minimum set of columns we can expect
        self.standard = standard

        # keywords to use for plotting
        self.plotkw = plotkw

        # make sure it's searchable via planet name
        self.standard.add_index('name')


    def __getitem__(self, key):
        '''
        Create a subpopulation by indexing or masking this one.
        '''
        subset = self.standard[key]
        cls = type(self)
        return cls(standard=subset, label=self.label, **self.plotkw)

    def single(self, name):
        '''
        Create a subpopulation of a single planet.
        '''

        subset = self.standard.loc[name]
        cls = type(self)
        return cls(standard=subset, label=name, **self.plotkw)

    def validate_columns(self):
        '''
        Make sure this standardized table has all the necessary columns.
        Summarize the amount of good data in each.
        '''

        N = len(self.standard)
        for k in necessary_columns:
            n = sum(self.standard[k].mask == False)
            self.speak(f'{k:>25} | {n}/{N} rows = {n/N:.0%} are OK!')


class Population(PopulationFromStandard):
    '''
    Population object keeps track of an exoplanet population.
    '''


    def __init__(self, label='exoplanets', remake=False, **plotkw):

        '''Initialize a population, by trying the following steps:
                #) load a standardized .npy
                1) load a standardized ASCII table
                2) ingest a raw table, and standardize it

        Parameters
        ----------
        
        '''

        # set the name for this population
        self.label = label


        try:
            # try to load the standardized table
            assert(remake == False)
            standard = self.load_standard()
            self.propagate()
        except (IOError,FileNotFoundError,AssertionError):
            # or create a new standardized table and save it
            self.ingestNew(**kwargs)

        PopulationFromStandard()


    @property
    def fileprefix(self):
        return clean(self.label)

    def loadRaw(self, remake=False):
        raw_numpy = directories['data'] + self.fileprefix + '_raw.npy'
        try:
            assert(remake==False)
            self.table = Table(np.load(raw_numpy))
            self.speak('loaded raw (but pre-saved) population from {0}'.format(raw_numpy))
        except IOError:
            self.loadFromScratch()
            np.save(raw_numpy, self.table)
            self.speak('saved raw population to {0} for faster future loading'.format(raw_numpy))

    def ingestNew(self, **kwargs):
        '''Ingest a new population table of arbitrary format,
            and then standardize it, using the tools defined in
            inherited population classes.'''


        # load the raw table
        self.loadRaw()

        # trim elements from raw table as necessary
        self.trimRaw()

        # create a standardized table from the array
        self.createStandard(**kwargs)

        # propagate into the attribute names
        self.propagate()

        # save the standardized table
        self.saveStandard()


    def load_standard(self, remake=False):
        '''Load a standardized population table, attempting...
            ...first from an .npy file (fast)
            ...then from a text file.'''

        standard_numpy = directories['data'] + self.fileprefix + '.npy'
        '''try:
            # first try to load the saved numpy table
            self.standard = Table(np.load(standard_numpy))
            self.speak('loaded standardized table from {0}'.format(standard_numpy))
        except IOError:'''
        # then try to load the table from a text file
        standard_ascii = directories['data'] + self.fileprefix + '.ascii'
        edited_ascii = standard_ascii.replace('.ascii', '.edited')

        kw = dict( delimiter='|', fill_values=[('',np.nan), ('--', np.nan)])#, format='basic', fill_values=[('--', 'nan'), (' ', 'nan')], data_start=1)
        try:
            self.speak('attempting to load {0}'.format(edited_ascii))
            self.standard = ascii.read(edited_ascii, **kw)
            self.speak('success!')
        except (IOError, FileNotFoundError):
            self.speak('failed!')
            self.speak('attempting to load {0}'.format(standard_ascii))
            self.standard = ascii.read(standard_ascii,**kw)
            self.speak('loaded (hopefully) standardized table from {0}'.format(standard_ascii))

        #    # and resave it as a numpy table (for faster loading next time)
        #    np.save(standard_numpy, self.standard)


        self.standard = np.ma.filled(self.standard, fill_value=np.nan)
        # add all of the table columns as attributes of the object
        self.propagate()

    def saveStandard(self):

        # KLUDGE!
        from .curation.Confirmed import correct
        correct(self)
        #standard_numpy = self.fileprefix + '.npy'
        standard_ascii = directories['data'] + self.fileprefix + '.ascii'

        # save it as an ascii table for humans to read
        self.standard.write(standard_ascii, format='ascii.fixed_width', bookend=False, delimiter='|', )#, include_names=columns)
        self.speak('wrote a standardized ascii table to {0}'.format(standard_ascii))

        # and resave it as a numpy table (for faster loading next time)
        #np.save(standard_numpy, self.standard)
        #self.speak('and saved standardized table to {0}'.format(standard_numpy))

    def removeRows(self, indices):
        self.standard.remove_rows(indices.nonzero()[0])
        self.propagate()


    def propagate(self):
        '''Link the columns of the standardized table as object attributes.'''
        for k in self.standard.colnames:
            self.__dict__[k] = self.standard[k]


    def find(self, name):
        '''Return index of a particular planet in the population.'''
        return np.array([clean(name) in clean(x) for x in self.name]).nonzero()[0]

    def correct(self, name, **kwargs):

        # find the entry to replace
        match = self.find(name)
        if len(match) != 1:
            self.speak('FAILED WHEN TRYING TO REPLACE {0}'.format(name))
            return
        assert(len(match) == 1)

        # loop over the keys
        for k, v in kwargs.items():
            self.speak('{0} used to be {1}'.format(k, v))
            self.standard[k][match] = v
            self.speak('  now it is {0}'.format(v))

    def __str__(self):
        return '<{0} | {1} planets>'.format(self.label, self.n)

    @property
    def n(self):
        return len(self.standard)

    @property
    def semimajoraxis(self):
        P = self.period*u.day
        G = con.G


        M = self.stellar_mass*u.Msun
        # kludge
        bad = np.isfinite(M) == False
        M[bad] = self.stellar_radius[bad]*u.Msun

        return ((G*P**2*M/4/np.pi**2)**(1./3.)).to('AU').value


    @property
    def absoluteJ(self):
        # use the effective temperature to select a J-band bolometric correction, then

        # pull out a proxy for the Sun from the Mamajek table
        #sun = mamajek.table['SpT'] == 'G2V'
        solar_stellar_teff = 5780

        # figure out the bolometric luminosities
        stellar_teffratio = self.stellar_teff/solar_stellar_teff
        radiusratio = self.stellar_radius
        luminosities = stellar_teffratio**4*radiusratio**2

        # SUUUUUUPER KLUDGE
        for i in range(2):
            try:
                test = self.stellar_teff.data
                assert(len(test) == len(self.stellar_teff))
                stellar_teff = np.array(test)
            except:
                pass

        # figure out the absolute J magnitude of a dwarf with this stellar_teff
        dwarf_absolutej = mamajek.tofrom('M_J')('stellar_teff')(stellar_teff)

        # figure out how bright a dwarf star of the same effective temperature would be
        dwarf_luminosities = 10**mamajek.tofrom('logL')('stellar_teff')(stellar_teff)

        # figure out how much brighter this star is than its stellar_teff-equivalent dwarf
        ratio = luminosities/dwarf_luminosities
        return dwarf_absolutej - 2.5*np.log10(ratio)



    @property
    def distance(self):

        distance = self.standard['stellar_distance'] + 0.0

        # SUUUUUUPER KLUDGE
        for i in range(2):
            try:
                test = distance.data
                assert(len(test) == len(distance))
                distance = np.array(test)
            except:
                pass
        bad = ((distance > 0.1) == False) + (np.isfinite(distance) == False)

        kludge = (self.stellar_teff == 0.0) + (np.isfinite(self.stellar_teff) == False)
        self.standard['stellar_teff'][kludge] = 5780


        distance[bad] = 10**(1 + 0.2*(self.J - self.absoluteJ))[bad]

        if self.__class__.__name__ != 'TESS':
            probablybad = distance == 10.0
            distance[probablybad] = np.nan
        return distance

    @property
    def teq(self):
        try:
            return self.standard['teq']
        except KeyError:
            a_over_r = self.a_over_r
            bad = (np.isfinite(a_over_r) == False) | (a_over_r == 0.0)


            try:
                stellar_mass = self.stellar_mass
            except AttributeError:
                #KLUDGE!
                stellar_mass = self.stellar_radius

            period = self.period
            stellar_radius = self.stellar_radius
            otherestimate = (craftroom.units.G*(period*craftroom.units.day)**2*(stellar_mass*craftroom.units.Msun)/4/np.pi**2/(stellar_radius*craftroom.units.Rsun)**3)**(1./3.)

            a_over_r[bad] = otherestimate[bad]

            '''plt.ion()
            plt.figure(self.label)
            plt.plot(a_over_r, otherestimate, '.')
            for i in range(len(a_over_r)):
                try:
                    plt.text(a_over_r[i], otherestimate[i], self.name[i])
                except:
                    pass'''


            return self.stellar_teff/np.sqrt(2*self.a_over_r)

    @property
    def insolation(self):
        u = craftroom.units
        a_over_rearth = u.au/u.Rsun
        teqearth = u.Tsun/np.sqrt(2*u.au/u.Rsun)
        return (self.teq/teqearth)**4

        #        u = craftroom.units
        #        return 1.0/self.a_over_r**2*self.stellar_teff**4/(u.Rsun/u.au)**2/u.Tsun**4


    @property
    def duration(self):

        d = self.standard['transit_duration']

        try:

            # try to calculate it from the parameters
            bad = np.isfinite(d)
            estimate = self.period/np.pi/self.a_over_r*np.sqrt(1-self.b**2)
            d[bad] = estimate[bad]


            # try to calculate it from the parameters
            bad = np.isfinite(d)
            estimate = self.period/np.pi/self.a_over_r
            d[bad] = estimate[bad]
        except TypeError:
            # try to calculate it from the parameters
            if np.isfinite(d):
                return d
            else:
                if np.isfinite(self.b):
                    b = self.b
                else:
                    b = 0.0
                return self.period/np.pi/self.a_over_r*np.sqrt(1-b**2)


            # try to calculate it from the parameters
            bad = np.isfinite(d)
            estimate = self.period/np.pi/self.a_over_r
            d[bad] = estimate[bad]

        return d#self.period/np.pi/self.a_over_r

    @property
    def photons(self):
        return 10**(-0.4*self.J)

    @property
    def surfacegravity(self):
        try:
            self.planet_mass
        except AttributeError:
            return 1000.0
        planet_mass = self.planet_mass
        g = craftroom.units.G*self.planet_mass*craftroom.units.Mearth/(self.planet_radius*craftroom.units.Rearth)**2
        g[g <= 0] = 2000.0
        g[g =='???'] = 0.0
        return g

        # for comparing with mass-less planets!

    @property
    def mu(self):
        # for comparing planets of different compositions
        return 2.32

    @property
    def scaleheight(self):
        return craftroom.units.k_B*self.teq/self.mu/craftroom.units.mp/self.surfacegravity

    @property
    def escape_velocity(self):
        e_grav = craftroom.units.G*self.planet_mass*craftroom.units.Mearth/self.planet_radius/craftroom.units.Rearth
        return np.sqrt(2*e_grav)


    @property
    def escape_parameter(self):
        e_thermal = craftroom.units.k_B*self.teq
        e_grav = craftroom.units.G*self.planet_mass*craftroom.units.Mearth*craftroom.units.mp/self.planet_radius/craftroom.units.Rearth
        return e_grav/e_thermal

    @property
    def noisepertransit(self):
        return 1.0/np.sqrt(self.photons*self.duration)

    @property
    def noisepertime(self):
        return 1.0/np.sqrt(self.photons/self.a_over_r)

    @property
    def noise(self):
        return 1.0/np.sqrt(self.photons)


    @property
    def transmissionsignal(self):
        H = self.scaleheight #cm
        Rp = self.planet_radius*craftroom.units.Rearth # cm
        Rs = self.stellar_radius*craftroom.units.Rsun
        return 2*H*Rp/Rs**2

    @property
    def emissionsignal(self):
        return (self.planet_radius/self.stellar_radius)**2*self.teq/self.stellar_teff

    @property
    def reflectionsignal(self):
        return self.depth/self.a_over_r**2


    @property
    def cheopsnoisepertransit(self):
        #150 ppm/min for a 9th magnitude star (in V)
        cheopsnoiseperminute =  150.0/1e6/np.sqrt(10**(-0.4*(self.V - 9.0)))
        durationinminutes = self.duration*24.0*60.0
        return cheopsnoiseperminute/np.sqrt(durationinminutes)

    @property
    def depth(self):
        return (self.planet_radius*craftroom.units.Rearth/self.stellar_radius/craftroom.units.Rsun)**2

    def plot(self, xname, yname, names=True, xlog=True, ylog=True):
        '''Plot one parameter against another.'''
        plt.ion()
        x, y = self.standard[xname], self.standard[yname]
        try:
            self.ax.cla()
        except:
            self.figure = plt.figure('Exoplanet Population')
            self.ax = plt.subplot()
        self.ax.set_xlabel(xname)
        self.ax.set_ylabel(yname)
        self.ax.plot(x, y, marker='o', alpha=0.5, color='gray', linewidth=0)
        if False:
            for i in range(len(x)):
                self.ax.text(x[i], y[i], self.table['NAME'][i])
        if xlog:
            plt.xscale('log')
        if ylog:
            plt.yscale('log')

        plt.draw()

    def thumbtack(self, maxr=1000, dr=100, labels=False):
        '''Plot the planets as thumbtacks.'''
        def scale(d):
            return np.array(d)**1.5
        r = scale(self.distance)
        x, y = r*np.cos(self.ra*np.pi/180), r*np.sin(self.ra*np.pi/180)
        plt.ion()
        plt.figure('thumbtacks')

        ax = plt.subplot()
        ax.cla()
        ax.set_aspect('equal')
        theta = np.linspace(0,2*np.pi,1000)
        angle = -90*np.pi/180

        gridkw = dict(alpha=0.25,  color='green')
        for originalradius in np.arange(dr,maxr*2,dr):
            radii = scale(originalradius)

            ax.plot(radii*np.cos(theta), radii*np.sin(theta), linewidth=3, **gridkw)
            ax.text(radii*np.cos(angle), radii*np.sin(angle), '{0:.0f} pc'.format(originalradius), rotation=90+ angle*180/np.pi, va='bottom', ha='center', size=13, weight='extra bold', **gridkw)

        ax.plot(x, y, marker='o', alpha=0.5, color='gray', linewidth=0, markeredgewidth=0)
        close = (self.name == 'WASP-94A b').nonzero()[0]#(self.distance < maxr).nonzero()[0]
        if labels:
            for c in close:
                plt.text(x[c], y[c], self.name[c])
        ax.set_xlim(-scale(maxr), scale(maxr))
        ax.set_ylim(-scale(maxr), scale(maxr))




    def compare(self, x='teq', y='radius', area='depth', color='stellar_radius'):

        xplot = self.__dict__[x]
        yplot = self.__dict__[y]
        sizeplot = self.__dict__[size]
        colorplot = self.__dict__[color]

        maxarea = 1000
        area = self.__dict__[area]
        sizeplot = np.sqrt(area/np.nanmax(area)*maxarea)

        plt.scatter(xplot, yplot, linewidth=0, marker='o', markersize=sizeplot)
