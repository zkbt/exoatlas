# general class for exoplanet populations
from ..imports import *
import string

exoplanet_columns = [
'name',
'ra', 'dec',
'stellar_distance',
'discoverer']

transit_columns = [
'period',
'semimajoraxis', #pl_orbsmax
'e', 'omega', 'inclination',
'transit_epoch',
'transit_duration',
'transit_depth',
'stellar_teff',
'stellar_mass',
'stellar_radius',
'planet_radius',
'planet_mass',
'transit_ar',
'transit_b']

necessary_columns = exoplanet_columns + transit_columns

desired_columns = [
'planet_mass_uncertainty_upper',
'planet_mass_uncertainty_lower',
'planet_radius_uncertainty_upper',
'planet_radius_uncertainty_lower',
'stellar_distance_uncertainty_upper',
'stellar_distance_uncertainty_lower']

# these are keywords that can be set for
default_plotkw = dict(color='black',
                      alpha=1,
                      zorder=0,
                      ink=True,
                      label_planets=False)

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
                          label=f'ExoplanetSubsets of {self.label}',
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
            try:
                # try to return the array of quantities (with units)
                return self.standard[key].quantity
            except TypeError:
                # columns without units don't have quantities
                return self.standard[key].data
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

    def uncertainty(self, key):
        '''
        Return an array of symmetric uncertainties on a column.

        Parameters
        ----------
        key : str
            The column for which we want errors.
        '''

        # first try for an `uncertainty_{key}` column
        try:
            return self.__getattr__(f'{key}_uncertainty')
        except (KeyError, AssertionError):
            # this can be removed after debugging
            self.speak(f'no symmetric uncertainties found for "{key}"')

        # then try for crudely averaging asymmetric uncertainties
        try:
            lower = self.__getattr__(f'{key}_uncertainty_lower')
            upper = self.__getattr__(f'{key}_uncertainty_upper')
            avg = 0.5*(np.abs(lower) + np.abs(upper))
            return avg
        except (KeyError, AssertionError):
            # this can be removed after debugging
            self.speak(f'no asymmetric uncertainties found for "{key}"')

        # then give up and return nans
        return np.nan*self.standard[key]

    def uncertainty_lowerupper(self, key):
        '''
        Return two arrays of lower and upper uncertainties on a column.

        Parameters
        ----------
        key : str
            The column for which we want errors.

        Returns
        -------
        lower : np.array
            The magnitude of the lower uncertainties (x_{-lower}^{+upper})
        upper : np.array
            The magnitude of the upper uncertainties (x_{-lower}^{+upper})
        '''

        # first try for actual asymmetric uncertainties
        try:
            lower = self.__getattr__(f'{key}_uncertainty_lower')
            upper = self.__getattr__(f'{key}_uncertainty_upper')
            return np.abs(lower), np.abs(upper)
        except (KeyError, AssertionError):
            # this can be removed after debugging
            self.speak(f'no asymmetric uncertainties found for "{key}"')

        # first try for an `uncertainty_{key}` column
        try:
            sym = self.__getattr__(f'{key}_uncertainty')
            return np.abs(sym), np.abs(sym)
        except (KeyError, AssertionError):
            # this can be removed after debugging
            self.speak(f'no symmetric uncertainties found for "{key}"')

        # then give up and return nans
        unc = np.nan*self.__getattr__(key)
        return unc, unc

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
        to create new sub-populations that extract subsets from this one.
        ''')

    @property
    def n(self):
        '''
        How many planets are in this population?
        '''
        return len(self.standard)

    @property
    def semimajoraxis(self):
        '''
        Have a safe way to calculate the semimajor axis of planets,
        that fills in gaps as necessary. Basic strategy:

            First from table.
            Then from NVK3L.
            Then from a/R*.

        '''

        # pull out the actual values from the table
        a = self.standard['semimajoraxis'].copy().quantity

        # try to replace bad ones with NVK3L
        bad = np.isfinite(a) == False
        self.speak(f'{sum(bad)}/{self.n} semimajoraxes are missing')

        # calculate from the period and the stellar mass
        P = self.period[bad]
        M = self.stellar_mass[bad]
        G = con.G
        a[bad] = ((G*M*P**2/4/np.pi**2)**(1/3)).to('AU')

        # replace those that are still bad with the a/R*
        stillbad = np.isfinite(a) == False
        self.speak(f'{sum(stillbad)}/{self.n} are still missing after NVK3L')
        # (pull from table to avoid potential for recursion)
        a_over_rs = self.standard['transit_ar'][stillbad].quantity
        rs = self.standard['stellar_radius'][stillbad].quantity
        a[stillbad] = a_over_rs*rs

        return a

    @property
    def a_over_rs(self):
        '''
        Have a safe way to calculate the scaled semimajor axis of planets,
        that fills in gaps as necessary. Basic strategy:

            First from table, mostly derived from transit.
            Then from the semimajor axis.
        '''

        # pull out the values from the table
        a_over_rs = self.standard['transit_ar'].copy()

        # try to replace bad ones with NVK3L
        bad = np.isfinite(a_over_rs) == False
        self.speak(f'{sum(bad)}/{self.n} values for a/R* are missing')

        a = self.semimajoraxis[bad]
        R = self.stellar_radius[bad]
        a_over_rs[bad] = a/R

        stillbad = np.isfinite(a_over_rs) == False
        self.speak(f'{sum(stillbad)}/{self.n} are still missing after a and R*')

        return a_over_rs

    @property
    def stellar_luminosity(self):
        T = self.stellar_teff
        R = self.stellar_radius
        sigma = con.sigma_sb
        return (4*np.pi*R**2*sigma*T**4).to(u.Lsun)

    @property
    def e(self):
        '''
        FIXME -- assumes are missing eccentricities are 0!
        '''

        # pull out the actual values from the table
        e = self.standard['e'].copy().quantity

        # try to replace bad ones with NVK3L
        bad = np.isfinite(e) == False
        self.speak(f'{sum(bad)}/{self.n} eccentricities are missing')
        self.speak(f'assuming they are all zero')
        e[bad] = 0

        return e

    @property
    def omega(self):
        '''
        (FIXME! we need better longitudes of periastron)
        '''

        # pull out the actual values from the table
        omega = self.standard['omega'].copy().quantity

        # try to replace bad ones with NVK3L
        bad = np.isfinite(omega) == False
        self.speak(f'{sum(bad)}/{self.n} longitudes of periastron are missing')
        e_zero = self.e == 0
        self.speak(f'{sum(e_zero)} have eccentricities assumed to be 0')
        omega[e_zero] = 0*u.deg

        return omega

    @property
    def b(self):
        '''
        Transit impact parameter.
        (FIXME! split this into transit and occultation)
        '''

        # pull out the actual values from the table
        b = self.standard['transit_b'].copy().quantity

        # try to replace bad ones with NVK3L
        bad = np.isfinite(b) == False
        self.speak(f'{sum(bad)}/{self.n} impact parameters are missing')

        # calculate from the period and the stellar mass
        a_over_rs = self.a_over_rs[bad]
        i = self.standard['inclination'][bad].quantity
        e = self.e[bad]
        omega = self.omega[bad]
        b[bad] = a_over_rs*np.cos(i)*((1-e**2)/(1+e*np.sin(omega)))

        # report those that are still bad
        stillbad = np.isfinite(b) == False
        self.speak(f'{sum(stillbad)}/{self.n} are still missing after using i')

        return b



    # the 1360 W/m^2 that Earth receives from the Sun
    earth_insolation = (1*u.Lsun/4/np.pi/u.AU**2).to(u.W/u.m**2)

    @property
    def insolation(self):
        '''
        The insolation the planet receives, relative to Earth.
        '''

        # calculate the average insolation the planet receives
        insolation = self.stellar_luminosity/4/np.pi/self.semimajoraxis**2
        return insolation.to(u.W/u.m**2)

    @property
    def teq(self):
        '''
        The equilibrium temperature of the planet.
        '''
        f = self.insolation
        sigma = con.sigma_sb
        A = 1
        return ((f*A/4/sigma)**(1/4)).to(u.K)

    @property
    def transit_depth(self):
        '''
        The depth of the transit
        (FIXME, clarify if this is 1.5-3.5 or what)
        '''

        # pull out the actual values from the table
        d = self.standard['transit_depth'].copy().quantity

        # try to replace bad ones with NVK3L
        bad = np.isfinite(d) == False
        self.speak(f'{sum(bad)}/{self.n} transit depths are missing')

        Rp = self.planet_radius[bad]
        Rs = self.stellar_radius[bad]

        d[bad] = (Rp/Rs).decompose()**2

        # report those that are still bad
        stillbad = np.isfinite(d) == False
        self.speak(f'{sum(stillbad)}/{self.n} are still missing after Rp/Rs')

        return d


    @property
    def transit_duration(self):
        '''
        The duration of the transit
        (FIXME, clarify if this is 1.5-3.5 or what)
        '''

        # pull out the actual values from the table
        d = self.standard['transit_duration'].copy().quantity

        # try to replace bad ones with NVK3L
        bad = np.isfinite(d) == False
        self.speak(f'{sum(bad)}/{self.n} transit durations are missing')


        P = self.period[bad]
        a_over_rs = self.a_over_rs[bad]
        b = self.b[bad]

        T0 = P/np.pi/a_over_rs
        T = T0*np.sqrt(1-b**2)

        e = self.e[bad]
        omega = self.omega[bad]
        factor = np.sqrt(1 - e**2)/(1 + e*np.sin(omega))

        d[bad] = (T*factor).to(u.day)

        # report those that are still bad
        stillbad = np.isfinite(d) == False
        self.speak(f'{sum(stillbad)}/{self.n} are still missing after P, a/R*, b')

        return d

    @property
    def photons(self):
        '''
        FIXME -- make this an actual flux?
        should it be a function for photons/s/m2/nm
        '''
        return 10**(-0.4*self.Jmag)

    @property
    def surface_gravity(self):
        '''
        (FIXME) -- make an assumption for planets without masses
        '''

        G = con.G
        M = self.planet_mass
        R = self.planet_radius

        g = (G*M/R**2).to('m/s**2')
        return g

    @property
    def mu(self):
        '''
        The mean molecular weight of an atmosphere.
        Here, assumed to be H/He-dominated at reasonable temperature.
        '''

        # for comparing planets of different compositions
        return 2.32

    @property
    def scale_height(self):
        '''
        The scale height of the atmosphere, at equilibrium temperature.
        '''
        k = con.k_B
        T = self.teq
        mu = self.mu
        m_p = con.m_p
        g = self.surface_gravity
        return k*T/mu/m_p/g

    @property
    def escape_velocity(self):
        '''
        The escape velocity of the planet.
        '''
        G = con.G
        M = self.planet_mass
        R = self.planet_radius
        return np.sqrt(2*G*M/R).to('m/s')


    @property
    def escape_parameter(self):
        '''
        The Jeans atmospheric escape parameter for atomic hydrogen,
        at the equilibrium temperature of the planet.
        '''
        k = con.k_B
        T = self.teq
        mu = 1
        m_p = con.m_p
        G = con.G
        M = self.planet_mass
        R = self.planet_radius

        e_thermal = k*T
        e_grav = G*M*m_p/R
        return e_grav/e_thermal

    @property
    def distance_modulus(self):
        mu = 5*np.log10(self.stellar_distance/(10*u.pc))
        return mu

    # PICK UP FROM HERE! THESE ARE ALL RELATIVE, SHOULD WE MAKE THEM ABSOLUTE?
    # (e.g. define a R=10-JWST-hour as m**2*s)
    @property
    def noisepertransit(self):
        return 1.0/np.sqrt(self.photons*self.transit_duration)

    @property
    def noisepertime(self):
        return 1.0/np.sqrt(self.photons/self.transit_ar)

    @property
    def noise(self):
        return 1.0/np.sqrt(self.photons)


    @property
    def transmissionsignal(self):
        H = self.scale_height #cm
        Rp = self.planet_radius*craftroom.units.Rearth # cm
        Rs = self.stellar_radius*craftroom.units.Rsun
        return 2*H*Rp/Rs**2

    @property
    def emissionsignal(self):
        return (self.planet_radius/self.stellar_radius)**2*self.teq/self.stellar_teff

    @property
    def reflectionsignal(self):
        return self.depth/self.transit_ar**2


    @property
    def cheopsnoisepertransit(self):
        #150 ppm/min for a 9th magnitude star (in V)
        cheopsnoiseperminute =  150.0/1e6/np.sqrt(10**(-0.4*(self.V - 9.0)))
        durationinminutes = self.transit_duration*24.0*60.0
        return cheopsnoiseperminute/np.sqrt(durationinminutes)


    def scatter(self, xname, yname, c=None, s=None, names=True, xlog=True, ylog=True, **kw):
        '''Plot one parameter against another.'''
        plt.ion()
        x, y = self.__getattr__(xname), self.__getattr__(yname)
        try:
            self.ax.cla()
        except:
            self.figure = plt.figure('Exoplanet Population')
            self.ax = plt.subplot()

        self.ax.set_xlabel(xname)
        self.ax.set_ylabel(yname)
        self.ax.scatter(x, y, c=c, s=s, **kw)
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
        r = scale(self.stellar_distance)
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
        close = (self.name == 'WASP-94A b').nonzero()[0]#(self.stellar_distance < maxr).nonzero()[0]
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
    def __init__(self, label='exoplanets', remake=False, skip_update=False, **plotkw):
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
        skip_update : bool
            Should we skip checking for updates in the existing data?
        **plotkw : dict
            All other keywords are stored as plotting suggestions.
        '''

        # set the name for this population
        self.label = label

        try:
            # try to load the standardized table
            assert(remake == False)
            standard = self.load_standard(skip_update=skip_update)
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

    def load_raw(self):
        raise NotImplementedError('''
        Yikes! The `.load_raw` method has not been defined
        for whatever object is trying to call it!
        ''')

    def trim_raw(self, raw):
        '''
        Trim bad/unnecessary rows out of a raw table of planet properties.
        '''

        # no trimming necessary
        trimmed = raw

        # for debugging, hang onto the trimmed table as a hidden attribute
        self._trimmed = trimmed

        # a trimmed table
        return self._trimmed


    def load_standard(self, skip_update=False):
        '''
        Load a standardized population table. Generally this
        will be from a file like ~/.exopop/standardized-*.txt

        Returns
        -------

        standard : astropy.table.Table
            A table of planet properties,
            with a minimal set of columns.
        skip_update : bool
            Should we skip checks to see if the data are too stale?
        '''

        # make sure this file is recent enough (unless we're skipping updates)
        if not skip_update:
            old = check_if_needs_updating(self.standard_path, self.expiration)
            assert(old == False)


        # keywords for reading a standardized table
        readkw = dict(format='ecsv', fill_values=[('',np.nan), ('--', np.nan)])

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
                            format='ascii.ecsv',
                            overwrite=True )
        self.speak(f'Saved a standardized text table to {self.standard_path}')
