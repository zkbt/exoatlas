"""define a Transit class, to represent an individual planetary transit"""
from ..imports import *


class Transit(Talker):
    """a Transit object is one transit of one (associated) planet"""

    def __init__(
        self, planet, i, phasefrommidtransit=0.0, observatory=None, buffer=1.5
    ):
        """
        initialize:
            planet = a Planet object (which contains a period and epoch)
            i = the epoch of this transit
            phasefrommidtransit = set this to 0.5 for eclipses (for circular)
            plan = the Plan object in which this is embedded
        """

        # set up the talking
        Talker.__init__(self)

        # store lots of links to other objects
        self.planet = planet
        self.planet._pithy = True
        self.coord = SkyCoord(self.planet.ra, self.planet.dec)
        self.observatory = observatory

        # set the midtransit time of this particular transit
        self.midtransit = (
            planet.period * (i + phasefrommidtransit) + planet.transit_epoch
        )[0]
        self.i = i

        # how much padding should be on the other side of the transit
        self.buffer = buffer

    @property
    def duration(self):
        """the duration of the transit"""

        # this should be in days
        d = self.planet.transit_duration[0]

        if np.isfinite(d):
            return d
        else:
            return 5.0 * u.hour  # KLUDGE!

    @property
    def ingress(self):
        """when is ingress?"""
        return self.midtransit - self.duration * 0.5

    @property
    def egress(self):
        """when is egress?"""
        return self.midtransit + self.duration * 0.5

    @property
    def pretransit(self):
        """when should observations start?"""
        return self.ingress - self.duration * self.buffer

    @property
    def posttransit(self):
        """when should observations end?"""

        return self.egress + self.duration * self.buffer

    @property
    def withbuffer(self):
        """list of plot_date times from start to finish of observation"""
        return [self.pretransit, self.posttransit]

    @property
    def justtransit(self):
        """list of plot_date times from start to finish of transit"""
        return [self.ingress, self.egress]

    def plot(
        self, alpha=0.6, linewidth=5, marker=None, color="black", **kwargs
    ):  # y=0,
        """plot this transit within a night"""

        kwargs["marker"] = marker
        kwargs["linewidth"] = linewidth
        kwargs["alpha"] = alpha
        # kwargs['color'] = self.planet.color

        # plot the transits
        # vertical = [y,y]

        # figure out the airmass for this transit
        times = self.midtransit + np.linspace(-1, 1, 100) * self.duration * self.buffer
        astropy_times = Time(times, format="jd")
        altaz = self.observatory.altaz(self.coord, astropy_times)
        airmass = altaz.secz.value
        ok = altaz.alt > 0

        # why is this here?
        try:
            assert self.duration < 0.5 * u.day
        except AssertionError:
            return

        # pull out the intransit subset
        intransit = (times >= self.ingress) * (times <= self.egress)

        # plot the full airmass curve
        plt.plot(astropy_times.plot_date[ok], airmass[ok], color=color, **kwargs)

        # plot the airmass curve of just the transit
        plt.plot(
            astropy_times.plot_date[ok * intransit],
            airmass[ok * intransit],
            color=color,
            **kwargs,
        )

        # add a label to this transit
        y = np.min(airmass[ok * intransit]) - 0.05
        description = f"{self.planet.name[0]}\n"
        description += (
            rf"{self.planet.transit_depth[0]:.2%},{self.planet.radius[0].value:.1f}$R_\oplus$,{self.planet.stellar_radius[0].value:.1f}$R_*$"
            + "\n"
        )
        description += rf"{self.planet.relative_insolation[0]:.0f}$S_\oplus$,{self.planet.distance[0]:.0f}"
        plt.text(
            Time(self.midtransit, format="jd").plot_date,
            y,
            description,
            color=color,
            ha="center",
            va="bottom",
            weight="bold",
            alpha=kwargs["alpha"],
            fontsize=5,
        )

    def airmass(self, time):
        """airmass string, for a particular time"""
        astropy_time = Time(time, format="jd")
        return "{0:.1f}".format(
            self.observatory.altaz(self.coord, astropy_time).secz.value[0]
        )

    def sunalt(self, time):
        """sunaltitude, for a particular time"""
        astropy_time = Time(time, format="jd")

        return "{0:+3.0f}".format(self.observatory.sun(astropy_time).alt.deg)

    def simpletime(self, time):
        astropy_time = Time(time, format="jd")
        return "{0:.5f} = {1} UT (A={2},S={3})".format(
            astropy_time.jd,
            astropy_time.iso[0:16],
            self.airmass(astropy_time),
            self.sunalt(astropy_time),
        )

    def details(self):
        s = f"{self.planet.name[0]:<20} \n"  # [not accounting for heliocentric]
        s += " {0:3.1f}x duration before = {1}\n".format(
            self.buffer, self.simpletime(self.pretransit)
        )
        s += "              ingress = {0}\n".format(self.simpletime(self.ingress))
        s += "          mid-transit = {0}\n".format(self.simpletime(self.midtransit))
        s += "               egress = {0}\n".format(self.simpletime(self.egress))
        s += "  {0:.1f}x duration after = {1}\n".format(
            self.buffer, self.simpletime(self.posttransit)
        )

        print(s)
        return s
        # print('mid-transit = {0:.5f}  [{1} to {2} UT, secz = {3:.2f} to {4:.2f}, sun = {5:+3.0f} to {6:+3.0f}]'.format(

        #                    self.observatory.sun(t.pretransit).alt.deg,
        #                    self.observatory.sun(t.posttransit).alt.deg))
