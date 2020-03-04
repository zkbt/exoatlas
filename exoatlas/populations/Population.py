# general class for exoplanet populations
from ..imports import *
from ..telescopes import *

import string

exocolumns = [
'name',
'ra', 'dec',
'distance',
'discoverer']

transit_columns = [
'period',
'semimajoraxis',
'e', 'omega', 'inclination',
'transit_epoch',
'transit_duration',
'transit_depth',
'stellar_teff',
'stellar_mass',
'stellar_radius',
'radius',
'mass',
'transit_ar',
'transit_b']

necessary_columns = exocolumns + transit_columns

desired_columns = [
'mass_uncertainty_upper',
'mass_uncertainty_lower',
'radius_uncertainty_upper',
'radius_uncertainty_lower',
'distance_uncertainty_upper',
'distance_uncertainty_lower']

# these are keywords that can be set for
default_plotkw = dict(color='black',
                      alpha=1,
                      zorder=0,
                      ink=True,
                      exact=False,
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

    #kludge?
    _pithy = True
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
            label = f'Subset of {self.label}'

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
            if isinstance(key, str):
                # remove spaces, to match the cleaned "name" index column
                key = key.replace(' ', '')
                label = key
            elif isinstance(key[0], str):
                # remove spaces, to match the cleaned "name" index column
                key = [k.replace(' ', '') for k in key]
                label = '+'.join(key)

            # pull out rows by planet name
            subset = self.standard.loc[key]
        # create a new population out of this subset
        return Population(standard=subset,
                          label=label,
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
            # try to get a plotkw from this pop, from the plotting defaults, from None
            try:
                assert(key in allowed_plotkw)
                return self.plotkw.get(key, default_plotkw[key])
            except (AssertionError, KeyError):
                raise AttributeError(f"""
                Alas, there seems to be no way to find `.{key}`
                as an attribute or propetry of {self}.
                """) #AtlasError

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
        except (KeyError, AssertionError, AtlasError, AttributeError): # is including AttributeError a kludge?
            # this can be removed after debugging
            self.speak(f'no symmetric uncertainties found for "{key}"')

        # then try for crudely averaging asymmetric uncertainties
        try:
            lower = self.__getattr__(f'{key}_uncertainty_lower')
            upper = self.__getattr__(f'{key}_uncertainty_upper')
            avg = 0.5*(np.abs(lower) + np.abs(upper))
            return avg
        except (KeyError, AssertionError, AtlasError, AttributeError):
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
        except (KeyError, AssertionError, AttributeError):
            # this can be removed after debugging
            self.speak(f'no asymmetric uncertainties found for "{key}"')

        # first try for an `uncertainty_{key}` column
        try:
            sym = self.__getattr__(f'{key}_uncertainty')
            return np.abs(sym), np.abs(sym)
        except (KeyError, AssertionError, AttributeError):
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

    def __len__(self):
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
        The insolation the planet receives, in W/m^2.
        '''

        # calculate the average insolation the planet receives
        insolation = self.stellar_luminosity/4/np.pi/self.semimajoraxis**2
        return insolation.to(u.W/u.m**2)

    @property
    def relative_insolation(self):
        '''
        The insolation the planet receives, relative to Earth.
        '''
        return self.insolation/self.earth_insolation

    @property
    def log_relative_insolation(self):
        return np.log10(self.relative_insolation)

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

        Rp = self.radius[bad]
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
    def kludge_mass(self):
        '''
        Have a safe way to calculate the mass of planets,
        that fills in gaps as necessary. Basic strategy:

            First from table.
            Then from msini.
        '''

        # pull out the actual values from the table
        M = self.standard['mass'].copy().quantity

        # try to replace bad ones with NVK3L
        bad = np.isfinite(M) == False
        self.speak(f'{sum(bad)}/{self.n} masses are missing')

        # estimate from the msini
        try:
            M[bad] = self.msini[bad]
        except (KeyError, AssertionError, AtlasError, AttributeError):
            pass

        # replace those that are still bad with the a/R*
        stillbad = np.isfinite(M) == False
        self.speak(f'{sum(stillbad)}/{self.n} are still missing after msini')

        return M

    @property
    def surface_gravity(self):
        '''
        (FIXME) -- make an assumption for planets without masses
        '''

        G = con.G
        M = self.mass
        R = self.radius

        g = (G*M/R**2).to('m/s**2')
        return g

    @property
    def density(self):
        '''
        The density of the planet.
        '''
        mass = self.mass
        volume = 4/3*np.pi*(self.radius)**3
        return (mass/volume).to('g/cm**3')

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
        return (k*T/mu/m_p/g).to('km')

    @property
    def escape_velocity(self):
        '''
        The escape velocity of the planet.
        '''
        G = con.G
        M = self.mass
        R = self.radius
        return np.sqrt(2*G*M/R).to('km/s')


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
        M = self.mass
        R = self.radius

        e_thermal = k*T
        e_grav = G*M*m_p/R
        return e_grav/e_thermal

    @property
    def distance_modulus(self):
        '''
        The distance modulus to the system, in magnitudes.
        '''
        mu = 5*np.log10(self.distance/(10*u.pc))
        return mu

    # FIXME -- how do we make it possible to change the threshold for inclusion?
    @property
    def transmission_signal(self, threshold=2):
        '''
        What is the transit depth of 1 scale height of an H2-rich
        atmosphere transiting in front of the star.

        Parameters
        ----------
        threshold : float
            By how many sigma must the planet mass be detected?

        '''
        with np.errstate(invalid='ignore'):

            H = self.scale_height
            Rp = self.radius
            Rs = self.stellar_radius
            depth = (2*H*Rp/Rs**2).decompose()

            dlnm = self.uncertainty('mass')/self.mass
            bad = dlnm > 1/threshold
            depth[bad] = np.nan
            return depth


    @property
    def reflection_signal(self):
        '''
        What is the reflected light eclipse depth,
        for an albedo of 100%?
        '''
        return 0.25*(self.radius/self.semimajoraxis).decompose()**2


    def emission_signal(self, wavelength=5*u.micron):
        '''
        What is the thermal emission eclipse depth,
        assuming Planck spectra for both star and planet?

        This calculation assumes a Bond albedo of 0
        and that heat is uniformly distributed over the planet.

        Parameters
        ----------
        wavelength : astropy.unit.Quantity
            The wavelength at which it should be calculated.
        '''

        # create thermal emission sources for both star and planet
        import rainbowconnection as rc
        star = rc.Thermal(teff=self.stellar_teff, radius=self.stellar_radius)
        planet = rc.Thermal(teff=self.teq, radius=self.radius)

        # calculate the depth as the luminosity ratio
        depths = planet.spectrum(wavelength)/star.spectrum(wavelength)


        return depths

    def stellar_brightness(self, wavelength=5*u.micron):
        '''
        How many photons/s/m^2/micron do we receive from the star?

        This is calculated from the distance, radius, and
        stellar effective temperature of the stars.

        (It could be potentially be improved with PHOENIX
        model grids and/or cleverness with photometry.)

        Parameters
        ----------
        wavelength : astropy.unit.Quantity
            The wavelength at which it should be calculated.
        '''

        # import some tools for easy cartoon spectra
        import rainbowconnection as rc

        # create source with right temperature, size, distance
        teff, radius = self.stellar_teff, self.stellar_radius
        star = rc.Thermal(teff=teff,
                          radius=radius).at(self.distance)

        # calculate the energy flux
        flux_in_energy = star.spectrum(wavelength)

        # convert to photon flux
        photon_energy = con.h*con.c/wavelength/u.ph
        flux_in_photons = flux_in_energy/photon_energy

        # return the
        return flux_in_photons.to('ph s^-1 m^-2 micron^-1')


    def photons_in_one_HST_orbit(self, wavelength=1.4*u.micron,
                                       R=20):
        '''
        The number of photons that hit the Hubble Space
        Telescope mirror, over the course of 1 HST orbit,
        accounting for 50% loss due to Earth occultations.
        '''


        duty_cycle=0.5
        dt = 1*u.def_unit('orbit',
                          96*u.minute,
                          doc=f'''
                          This custom unit represents the
                          the time of one HST orbit.
                          ''')

        telescope_unit = define_telescope_unit_by_name('HST',
                                      wavelength=wavelength,
                                      R=R,
                                      dt=dt)

        collecting_power = telescope_unit*duty_cycle

        b = self.stellar_brightness(wavelength=wavelength)

        return (b*collecting_power).decompose()

    def photons_in_one_WFC3_orbit(self, **kw):
        N = self.photons_in_one_HST_orbit(**kw)
        throughput = 0.4
        duty_cycle = 0.8
        return N*throughput*duty_cycle

    def depth_uncertainty_in_one_WFC3_orbit(self, **kw):
        N = self.photons_in_one_WFC3_orbit(**kw)/u.ph
        sigma_D = 1/np.sqrt(N)
        oot = np.sqrt(2)
        return sigma_D*oot

    def radius_uncertainty_in_one_WFC3_orbit(self, **kw):
        Rs = self.stellar_radius
        D = self.transit_depth
        sigma_D = self.depth_uncertainty_in_one_WFC3_orbit(**kw)
        sigma_Rp = (Rs*sigma_D/2/np.sqrt(D)).to(u.Rearth)
        return sigma_Rp

    def radius_uncertainty_in_one_WFC3_transit(self, **kw):
        N_orbits = (self.transit_duration/(96*u.minute)).decompose()
        per_orbit = self.radius_uncertainty_in_one_WFC3_orbit(**kw)
        return per_orbit/np.sqrt(N_orbits)




    def stellar_brightness_JWST(self, wavelength=5*u.micron):
        '''
        The stellar brightness, converted to JWST units.

        Parameters
        ----------
        wavelength : astropy.unit.Quantity
            The wavelength at which it should be calculated.
        '''

        flux_in_photons = self.stellar_brightness(wavelength)
        unit = photon_unit/JWST_transit_unit(wavelength)
        return flux_in_photons.to(unit)

    def stellar_brightness_HST(self, wavelength=1.4*u.micron):
        '''
        The stellar brightness, converted to HST units.

        Parameters
        ----------
        wavelength : astropy.unit.Quantity
            The wavelength at which it should be calculated.
        '''

        flux_in_photons = self.stellar_brightness(wavelength)
        unit = photon_unit/HST_orbit_unit(wavelength)
        return flux_in_photons.to(unit)

    # PICK UP FROM HERE! THESE ARE ALL RELATIVE, SHOULD WE MAKE THEM ABSOLUTE?
    # (e.g. define a R=20-JWST-hour as m**2*s)
    @property
    def noisepertransit(self):
        raise NotImplementedError('!')
        return 1.0/np.sqrt(self.photons*self.transit_duration)

    @property
    def noisepertime(self):
        raise NotImplementedError('!')
        return 1.0/np.sqrt(self.photons/self.transit_ar)

    @property
    def noise(self):
        return 1.0/np.sqrt(self.photons)


    @property
    def cheopsnoisepertransit(self):
        raise NotImplementedError('!')
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
        will be from a file like ~/.exoatlas/standardized-*.txt

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
        like ~/exoatlas/standardized-*.txt
        '''

        # save it as an ascii table for humans to read
        standard.write(self.standard_path,
                            format='ascii.ecsv',
                            overwrite=True )
        self.speak(f'Saved a standardized text table to {self.standard_path}')
