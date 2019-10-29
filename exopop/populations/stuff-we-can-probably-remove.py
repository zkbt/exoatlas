

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
