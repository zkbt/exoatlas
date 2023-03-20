from ..imports import *
from .Exoplanets import *
from .curation.TransitingExoplanets import curate

__all__ = ["TransitingExoplanets"]


class TransitingExoplanets(Exoplanets):
    def __init__(self, **kw):
        Exoplanets.__init__(self, **kw)

        # set the label
        self.label = "Transiting Exoplanets"

        # remove non-transiting
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            ok_transit = self.detected_in_transit == 1
            has_no_radius = np.isnan(self.radius)

        # apply any planet-specific curation to the population
        curate(self)

        # tidy up this population
        self.remove_nontransiting()
        self.remove_bad_radii()
        self.tidy_eccentricities()

    def remove_nontransiting(self):
        """
        Remove non-transiting planets from the population.

        Returns
        -------
        (modifies the input population in place, but returns nothing)
        """
        ok = self.detected_in_transit == 1
        self.standard = self.standard[ok]

    def remove_bad_radii(self):
        """
        Remove planets with bad radii.

        Returns
        -------
        (modifies the input population in place, but returns nothing)
        """

        self.speak("Trying to remove bad radii.")

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            has_no_radius = np.isnan(self.radius)
            upper_limit_on_radius = 30 * u.Rearth
            ok = (self.radius < upper_limit_on_radius) | has_no_radius
            bad = ok == False

        N = sum(bad)
        s = ", ".join(self.name[bad])
        self.speak(
            f"""
        {N} planets had non-existent radii, or radii larger
        than {upper_limit_on_radius}. Let's assume those are bad.
        """
        )
        if N > 0:
            self.speak(f"The affected planets are {s}")

        self.standard = self.standard[ok]

    def tidy_eccentricities(self):
        """
        Often in the exoplanet archive, a planet with e!=0 is in fact an
        upper limit. We make the (probably a little too bold) assumption
        that if an eccentricity has actually been determined for a planet
        then it will also have finite omega. Therefore, we'll call everything
        that has an eccentricity but no omega an upper limit and force it
        to zero. FIXME -- we should handle this more responsibly in the future!

        Returns
        -------
        (modifies the input population in place, but returns nothing)
        """

        self.speak("Trying to tidy the eccentricities.")

        # figure out which e's are troublesome
        finite_e = self.e > 0
        no_omega = np.isfinite(self.omega) == False
        is_probably_upper_limit = finite_e & no_omega

        N = sum(is_probably_upper_limit)
        s = ", ".join(self.name[is_probably_upper_limit])
        self.speak(
            f"""
        {N} eccentricities were finite but had no omega.
        Let's assume those are upper limits and set e=0
        for all of them.
        """
        )

        if N > 0:
            self.speak(f"The affected planets are {s}")

        # force those ones to zero
        self.standard["eccentricity"][is_probably_upper_limit] = 0
