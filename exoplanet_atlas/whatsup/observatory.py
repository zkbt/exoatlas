from ..imports import *
from pytz import timezone
from astropy.coordinates import EarthLocation, get_sun, AltAz
import matplotlib.dates as mdates

# create an empty dictionary
observatories = {}

# define a few observatories
#  parameters copied from skycalc.c source
observatories["CTIO"] = dict(
    name="Cerro Tololo",
    timezone="Chilean",
    standardzone=4.0 * u.hour,
    usedaylightsaving=-1,
    longitude=4.721 * u.hourangle,
    latitude=-30.165 * u.deg,
    elevsea=2215.0 * u.m,
    elev=2215.0 * u.m,  # /* for ocean horizon, not Andes! */
)

observatories["LCO"] = dict(
    name="Las Campanas Observatory",
    timezone="Chilean",
    standardzone=4.0 * u.hour,
    usedaylightsaving=-1,
    longitude=4.71333 * u.hourangle,
    latitude=-29.00833 * u.deg,
    elevsea=2282.0 * u.m,
    elev=2282.0 * u.m,  # /* for ocean horizon, not Andes! */
)

observatories["FLWO"] = dict(
    name="Mt. Hopkins, AZ",
    timezone="Mountain",
    standardzone=7.0 * u.hour,
    usedaylightsaving=0,
    longitude=7.39233 * u.hourangle,
    latitude=31.6883 * u.deg,
    elevsea=2608.0 * u.m,
    elev=500.0 * u.m,  # /* approximate elevation above horizon mtns */
)

observatories["APO"] = dict(
    name="Apache Point Observatory",
    timezone="Mountain",
    standardzone=7.0 * u.hour,
    usedaylightsaving=0,
    longitude=105 * u.deg,
    latitude=32 * u.deg,
    elevsea=2798.0 * u.m,
    elev=500.0 * u.m,  # /* approximate elevation above horizon mtns */
)


observatories["SBO"] = dict(
    name="Sommers-Bausch Observatory",
    timezone="Mountain",
    standardzone=7.0 * u.hour,
    usedaylightsaving=0,
    longitude=105.2635 * u.deg,
    latitude=40.00372 * u.deg,
    elevsea=1653.0 * u.m,
    elev=500.0 * u.m,  # /* approximate elevation above horizon mtns */
)


class Observatory(Talker):
    """
    The Observatory object, for keeping track of timing calculating alt-az.

    :param name:
        String listing short name (LCO, CTIO, FLWO) of obsrvatory.

    """

    def __repr__(self):
        return f"<{self.name}>"

    def __init__(self, abbreviation=None, plan=None):
        """
        Initialize an observatory.
        """

        Talker.__init__(self)

        # make sure an abbreviation is defined
        if abbreviation is None:
            self.speak("Pick an observatory from following options:")
            self.speak(str(observatories.keys()))
            abbreviation = self.input().strip()

        # populate the attributes of this observatory
        o = observatories[abbreviation]
        for k in o.keys():
            self.__dict__[k] = o[k]

        # set this observatory's location on Earth
        self.location = EarthLocation(
            lat=self.latitude, lon=-self.longitude, height=self.elevsea
        )

        # connect this observatory to a plan
        self.plan = plan

    def sun(self, times):
        """
        Get the alt/az of the Sun (at some times).
        """

        sunCoords = get_sun(times)
        altAzFrame = AltAz(obstime=times, location=self.location)
        return sunCoords.transform_to(altAzFrame)

    def altaz(self, coord, times):
        """
        Get the altaz of the star (at some times).
        """

        frame = AltAz(obstime=times, location=self.location)
        return coord.transform_to(frame)

    def plotSun(self, jd, ax=None, threshold=-12.0, color="black"):
        """
        Plot the alitude of the Sun, for some given times.
        """

        self.speak("plotting the sun")

        times = Time(jd, format="jd")
        sunAltAz = self.sun(times)

        dt = times[1:] - times[:-1]
        nudged = times[:-1] + dt / 2.0
        a = sunAltAz.alt.deg - threshold
        sunrises = ((a[1:] > 0) * (a[:-1] < 0)).nonzero()
        sunsets = ((a[1:] < 0) * (a[:-1] > 0)).nonzero()
        starts = nudged[sunrises]
        finishes = nudged[sunsets]
        maxlength = min(len(starts), len(finishes))
        if starts[0] > finishes[0]:
            starts, finishes = starts[: maxlength - 1], finishes[1:maxlength]
        else:
            starts, finishes = starts[:maxlength], finishes[:maxlength]

        assert (finishes > starts).all()
        for i in range(len(starts)):
            plt.axvspan(
                starts[i].plot_date,
                finishes[i].plot_date,
                color=color,
                alpha=0.8,
                zorder=10000,
            )

    def plot_empty_night(self, midnight, threshold=-12.0, ax=None):
        astropy_midnight = Time(midnight, format="jd")

        plt.figure(figsize=(8, 3), dpi=400)
        times = Time(
            np.arange(midnight.value - 1.0, midnight.value + 1.0, 5.0 / 60 / 24),
            format="jd",
        )

        sunAltAz = self.sun(times)

        dt = times[1:] - times[:-1]
        nudged = times[:-1] + dt / 2.0
        a = sunAltAz.alt.deg - threshold
        sunrises = ((a[1:] > 0) * (a[:-1] < 0)).nonzero()
        sunsets = ((a[1:] < 0) * (a[:-1] > 0)).nonzero()
        starts = nudged[sunrises]
        finishes = nudged[sunsets]
        maxlength = min(len(starts), len(finishes))
        if starts[0] > finishes[0]:
            starts, finishes = starts[: maxlength - 1], finishes[1:maxlength]
        else:
            starts, finishes = starts[:maxlength], finishes[:maxlength]

        assert (finishes > starts).all()
        for i in range(len(starts)):
            plt.axvspan(
                starts[i].plot_date,
                finishes[i].plot_date,
                color="black",
                alpha=0.8,
                zorder=10000,
            )

        # set the formatting on the time label
        fmt = plt.matplotlib.dates.DateFormatter("UT %m-%d|%H")
        hours = mdates.HourLocator(interval=1)  #
        ax = plt.gca()
        ax.xaxis.set_major_formatter(fmt)
        ax.xaxis.set_major_locator(hours)
        ax.tick_params(axis="x", labelsize=6.5)
        plt.xlim(astropy_midnight.plot_date - 0.5, astropy_midnight.plot_date + 0.5)
        plt.ylabel("Airmass")

        f = plt.gcf()
        f.autofmt_xdate(rotation=45)

    def plotAirmass(self, coord, times, **kwargs):
        """
        Plot the airmass of the object (at some times).
        """

        self.speak("plotting airmass")
        # calculate altaz
        altaz = self.altaz(coord, times)
        airmass = altaz.secz.value
        ok = altaz.alt > 0

        self.plan.ax["airmass"].plot(times.plot_date[ok], airmass[ok], **kwargs)
        self.plan.ax["airmass"].set_ylim(3, 1)
