# general class for exoplanet populations
from ..imports import *
import string

# the mamajek relation is needed for estimating distances to stars, for those without
mamajek = thistothat.Mamajek()
mamajek._pithy = True

exoplanet_columns = [
'name',
'ra', 'dec',
'stellar_distance',
'J',
'discoverer']

transit_columns = [
'period',
'transit_epoch',
'transit_duration',
'transit_depth',
'stellar_teff',
'stellar_mass',
'stellar_radius',
'planet_radius',
'planet_mass',
'a_over_r',
'b']

necessary_columns = exoplanet_columns + transit_columns

desired_columns = [
'planet_mass_upper',
'planet_mass_lower',
'planet_radius_upper',
'planet_radius_lower',
'stellar_distance_upper',
'stellar_distance_lower']

# these are keywords that can be set for
default_plotkw = dict(color='black',
                      alpha=1,
                      zorder=0,
                      ink=True)

# what keywords can we set for the population plotkw?
allowed_plotkw = list(default_plotkw.keys())
allowed_plotkw += ['s',
                   'c',
                   'marker',
                   'cmap',
                   'norm',
                   'vmin',
                   'vmax'
                   'linewidths',
                   'edgecolors',
                   'facecolors']

class Population(Talker):
    '''
    Create a population from a standardized table.
    '''

    def __init__(self, standard, label='unknown', **plotkw):
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
        self.standard = Table(standard)

        # store a label for this population
        self.label = label

        # keywords to use for plotting
        self.plotkw = plotkw

        # make sure it's searchable via planet name
        self.standard.add_index('name')


    def __getitem__(self, key):
        '''
        Create a subpopulation of planets by indexing, slicing, or masking.
        '''

        try:
            # if the key is an index/slice/mask, return it
            subset = self.standard[key]

            # if the key is a column, raise an error
            if type(key) in self.standard.colnames:
                raise IndexError(f'''
                You seem to be trying to access a column from this
                population via `pop[{key}]`. For clarity, all `[]`
                indexing is reserved for selecting subsets of the
                population.

                To access your particular column, please try either
                `pop.{key}` or `pop.standard[{key}]` to return a
                1D array of the entries in that column.
                ''')
        except KeyError:
            # use a string or a list of strings to index the population by name
            if type(key) == str:
                # remove spaces, to match the cleaned "name" index column
                key = key.replace(' ', '')
            elif type(key[0]) == str:
                # remove spaces, to match the cleaned "name" index column
                key = [k.replace(' ', '') for k in key]
            # pull out rows by planet name
            subset = self.standard.loc[key]

        # create a new population out of this subset
        return Population(standard=subset,
                          label=f'Subset of {self.label}',
                          **self.plotkw)

    def __getattr__(self, key):
        '''
        If an attribute/method isn't defined for a population,
        look for it as a column of the standardized table.

        For example, `population.stellar_radius` will try to
        access `population.standard['stellar_radius']`.

        Parameters
        ----------
        key : str
            The attribute we're trying to get.
        '''
        if key == 'label':
            raise RuntimeError('Yikes!')
        try:
            # extract the column from the standardized table
            return self.standard[key]
        except KeyError:
            # try to get a plotkw from this pop, from the default, then None
            try:
                assert(key in allowed_plotkw)
                return self.plotkw.get(key, default_plotkw[key])
            except KeyError:
                return None

    def __setattr__(self, key, value):
        '''
        Define what happens when we try to set an attribute via `pop.attr = x`.
        If the keyword is a pre-defined "plotting" keyword in `allowed_plotkw`,
        then we should save it in a special `plotkw` dictionary. Otherwise,
        the attribute should be set as normal.

        Parameters
        ----------
        key : str
            The attribute we're trying to set.
        value : anything
            The value we're trying to give that attribute.
        '''

        if key in allowed_plotkw:
            # store plotting keywords in a separate plotting dictionary
            self.plotkw[key] = value
        else:
            # otherwise, store attributes as normal for objects
            self.__dict__[key] = value

    def __repr__(self):
        '''
        How should this object appear as a repr/str?
        '''
        return f'<{self.label} | population of {self.n} planets>'

    def single(self, name):
        '''
        Create a subpopulation of a single planet.
        '''

        # create a subset of the standardized table
        subset = self.standard.loc[name]

        # create a new object, from this subset
        return Population(standard=subset, label=name, **self.plotkw)

    def validate_columns(self):
        '''
        Make sure this standardized table has all the necessary columns.
        Summarize the amount of good data in each.
        '''

        N = len(self.standard)
        for k in necessary_columns:
            try:
                n = sum(self.standard[k].mask == False)
            except AttributeError:
                try:
                    n = sum(np.isfinite(self.standard[k]))
                except TypeError:
                    n = sum(self.standard[k] != '')
            self.speak(f'{k:>25} | {n:4}/{N} rows = {n/N:4.0%} are not empty')

    def find(self, name):
        '''
        Return index of a particular planet in the population.

        ??? = maybe this could/should be replaced with some table cleverness?
        '''

        return np.array([clean(name) in clean(x) for x in self.name]).nonzero()[0]

    def correct(self, name, **kwargs):
        '''
        Correct the properties of a particular planet,
        modifying its values in the standardized table.

        Parameters
        ----------
        name : str
            The name of the planet to fix.
        **kwargs : dict
            Keyword arguments will go into modifying
            the properties of that planet.
        '''

        # find the entry to replace
        match = self.find(name)
        if len(match) != 1:
            self.speak(f'failed when trying to modify parameters for {name}')
            return

        # loop over the keys, modifying each
        self.speak(f'for planet "{name}"')
        for k, new in kwargs.items():
            old = self.standard[k][match]
            self.speak(f' {k} changed from {old} to {new}')
            self.standard[k][match] = new

    def removeRows(self, indices):

        raise NotImplementedError('''
        The `removeRows` method has been removed. Please use something like
        `population[0:42]` or `population[ok]` to use slices, indices, or masks
        to pull subsets out of this population.
        ''')

    @property
    def n(self):
        '''
        How many planets are in this population?
        '''
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


class PredefinedPopulation(Population):
    '''
    Population object keeps track of an exoplanet population.
    '''

    expiration = 0.00001
    def __init__(self, label='exoplanets', remake=False, **plotkw):
        '''
        Initialize a population, by trying the following steps:
                1) Load a standardized ascii table.
                2) Ingest a raw table, and standardize it.

        Parameters
        ----------
        label : str
            The name of this population, for use both in filenames
            and labeling points on plots.
        remake : bool
            Should we re-ingest this table from its raw ingredients?
        **plotkw : dict
            All other keywords are stored as plotting suggestions.
        '''

        # set the name for this population
        self.label = label

        try:
            # try to load the standardized table
            assert(remake == False)
            standard = self.load_standard()
        except (IOError,FileNotFoundError,AssertionError):
            # or create a new standardized table and save it
            standard = self.ingest_table(remake=remake)

        # initialize with a standard table
        Population.__init__(self,
                            standard=standard,
                            label=label,
                            **plotkw)

    @property
    def fileprefix(self):
        '''
        Define a fileprefix for this population, to be used
        for setting the filename of the standardized population.
        '''
        return clean(self.label)

    def ingest_table(self, **kwargs):
        '''
        Ingest a new population table of arbitrary format,
        and then standardize it, using the tools defined in
        inherited population classes.'''


        # load the raw table
        raw = self.load_raw()

        # trim elements from raw table as necessary
        trimmed = self.trim_raw(raw)

        # create a standardized table from the array
        standard = self.create_standard(trimmed)

        # save the standardized table
        self.save_standard(standard)

        return standard

    @property
    def standard_path(self):
        '''
        Define the filepath for the standardized table.
        '''
        return os.path.join(directories['data'],
                            f'standardized-{self.fileprefix}.txt')

    def load_standard(self):
        '''
        Load a standardized population table. Generally this
        will be from a file like ~/.exopop/standardized-*.txt

        Returns
        -------

        standard : astropy.table.Table
            A table of planet properties,
            with a minimal set of columns.
        '''

        # make sure this file is recent enough
        old = check_if_needs_updating(self.standard_path, self.expiration)
        assert(old == False)


        # keywords for reading a standardized table
        readkw = dict(delimiter='|',
                      fill_values=[('',np.nan), ('--', np.nan)])

        standard = ascii.read(self.standard_path, **readkw)
        self.speak(f'Loaded standardized table from {self.standard_path}')

        # ??? change this to do something more clever with tables
        # masked = np.ma.filled(standard, fill_value = np.nan)

        return standard

    def save_standard(self, standard):
        '''
        Save the standardized table out to a text file
        like ~/exopop/standardized-*.txt
        '''

        # save it as an ascii table for humans to read
        standard.write(self.standard_path,
                            format='ascii.fixed_width',
                            bookend=False,
                            delimiter='|',
                            overwrite=True )
        self.speak(f'Saved a standardized text table to {self.standard_path}')
