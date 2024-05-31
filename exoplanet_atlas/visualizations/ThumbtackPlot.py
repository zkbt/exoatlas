from .imports import *
from .Panels import BubblePanel
from .TransitingExoplanets import NonKepler, Kepler, TESS
from .KOI import UnconfirmedKepler
import datetime

# the figure size
figsize = 7

# the aspect ratio of the figure
aspect = 768 / 1024.0


# the angle of the distance labels
angle = 30 * np.pi / 180  # 67.5


class ThumbtackPlot(BubblePanel):
    """Plot exoplanet populations on a (possibly animated) thumbtack plot."""

    def __init__(self, pops, lightyears=False, **kwargs):
        """initialize the thumbtack object"""
        self.lightyears = lightyears

        if self.lightyears:
            self.maxlabeldistance = 100.0 * craftroom.units.ly / craftroom.units.pc
        else:
            self.maxlabeldistance = 30.0
        self.planetfontsize = 6
        BubblePanel.__init__(self, pops=pops, **kwargs)
        self.title = "Thumbtacks"
        self.xlabel = ""
        self.ylabel = ""
        self.circlegrid = [3, 10, 30, 100, 300, 1000]
        if self.lightyears:
            self.circlegrid = (
                np.array([3, 10, 30, 100, 300, 1000, 3000])
                * craftroom.units.ly
                / craftroom.units.pc
            )

    def get_sizes(self):
        """by default, all points are the same size"""
        return 50.0 / 1000.0

    @property
    def r(self):
        """convert the system's distance into a radial coordinate"""
        return self.stretch(self.pop.distance)

    @property
    def theta(self):
        """convert RA to theta for plotting"""
        return self.toclock(self.pop.ra * 24.0 / 360.0)

    @property
    def x(self):
        """convert distance and RA to y"""
        return self.stretch(self.pop.distance) * np.cos(self.theta)

    @property
    def y(self):
        """convert distance and RA to x"""
        return self.stretch(self.pop.distance) * np.sin(self.theta)

    def setup(self):
        """setup the initial plotting window"""

        # set the figure's aspect ratio
        self.aspect = aspect

        # make the figure
        plt.figure(self.title, figsize=(figsize, figsize * self.aspect), dpi=200)

        # populate with an axis that fills entire frame
        gs = plt.matplotlib.gridspec.GridSpec(
            1, 1, left=0, right=1, bottom=0, top=1, hspace=0, wspace=0
        )
        self.ax = plt.subplot(gs[:, :])

        # make sure x and y directions are plotted on same scale
        self.ax.set_aspect("equal")

        # make the lines go away
        self.ax.get_xaxis().set_visible(False)
        self.ax.get_yaxis().set_visible(False)

        # add RA labels around the clock face
        hours = np.arange(0, 24, 3)
        theta = self.toclock(hours)
        months = ["September", "December", "March", "June"]
        clockr = 0.925
        for i in range(len(hours)):
            x = clockr * (np.cos(theta[i]) * 0.5) * self.aspect + 0.5
            y = clockr * np.sin(theta[i]) * 0.5 + 0.5
            self.ax.text(
                x,
                y,
                "{0}h".format(hours[i]),
                rotation=-90 + theta[i] * 180 / np.pi,
                transform=self.ax.transAxes,
                va="center",
                ha="center",
                size=10,
                alpha=0.5,
            )
            if (i % 2) == 0:
                self.ax.text(
                    x,
                    y,
                    "\n\n" + months[np.int(i / 2)],
                    size=6,
                    rotation=-90 + theta[i] * 180 / np.pi,
                    transform=self.ax.transAxes,
                    va="center",
                    ha="center",
                    alpha=0.35,
                )

    def circles(self, radii):
        """plot circles at particular radii"""

        # create grid of angles
        theta = np.linspace(0, 2 * np.pi, 1000)

        # graphics keywords for plotting the circle lines
        gridkw = dict(color="gray", alpha=0.5, zorder=-100)

        # store the plotted circles, so they can be erased elsewhere
        self.circlelabels = {}

        # nudge=1.15
        nudge = 0.0
        for originalradius in radii:
            # get the self.stretched radial coordinate
            r = self.stretch(originalradius)

            # plot the circle
            self.ax.plot(r * np.cos(theta), r * np.sin(theta), linewidth=4, **gridkw)

            # label the circles with their distances
            if self.lightyears:
                circletext = "{0:.0f} ly".format(
                    originalradius * craftroom.units.pc / craftroom.units.ly
                )
            else:
                circletext = "{0:.0f} pc".format(originalradius)
            self.circlelabels[originalradius] = self.ax.text(
                (nudge + r) * np.cos(angle),
                (nudge + r) * np.sin(angle),
                circletext,
                rotation=-90 + angle * 180 / np.pi,
                va="center",
                ha="center",
                size=13,
                weight="extra bold",
                **gridkw
            )

    def build(self, keys=None, distances=[10, 30, 100, 300, 1000]):
        try:
            del self.ax
            del self.signature
        except AttributeError:
            pass
        plt.cla()
        if keys is None:
            keys = self.pops.keys()
        for key in keys:
            self.plot(key)
            for z in distances:
                self.zoom(z)
                self.clearnames()
                for thiskey in keys:
                    self.namestars(thiskey)
                plt.draw()
                plt.savefig(self.fileprefix(key))
        return self

    def movie(
        self,
        step=1.02,
        bitrate=20000,
        highlight="",
        keys=None,
        maxdistance=1500.0,
        fileprefix="exoplanetszoom",
    ):
        metadata = dict(
            title="Exoplanets Zoom", artist="Zach Berta-Thompson (zkbt@mit.edu)"
        )
        self.writer = animation.FFMpegWriter(fps=15, metadata=metadata, bitrate=bitrate)

        try:
            del self.ax
            del self.signature

        except AttributeError:
            pass
        plt.cla()
        if keys is None:
            keys = self.pops.keys()
        for key in keys:
            self.plot(key)
            if highlight == "completeness":
                self.highlight(
                    (self.pop.radius > 2)
                    * (self.pop.radius < 4)
                    * (self.pop.period < 10),
                    "2-4 Earth radii\nP<10 days",
                )
            if highlight == "small":
                self.highlight((self.pop.radius < 2), "<2 Earth radii")
            if highlight == "habitable":
                self.highlight(
                    (self.pop.radius > 0.7)
                    * (self.pop.radius < 1.6)
                    * (self.pop.teq < 310)
                    * (self.pop.teq > 200),
                    "Potentially Habitable Planets",
                )

        f = plt.gcf()
        filename = "{}_{}{}.mp4".format(fileprefix, "+".join(keys), highlight)
        self.speak("writing movie to {0}".format(filename))
        z = 2.0

        with self.writer.saving(f, filename, 1024.0 / figsize):
            # loop over exposures
            while z < maxdistance:
                self.zoom(z)

                #                    self.clearnames()
                #                    if highlight != 'habitable':
                #                        self.namestars('nonkepler')
                #                        #KLUDGE! self.namestars('new')
                #                        #if key == 'kepler':
                #                        self.namestars('kepler')

                self.writer.grab_frame()
                self.speak("zoomed to {0}".format(z))
                z *= step
                # plt.draw()
            # plt.savefig(self.fileprefix(key))

    def zoom(self, distance):
        self.outer = distance
        height = self.stretch(self.outer) * 1.2
        width = height / self.aspect
        self.ax.set_ylim(-height, height)
        self.ax.set_xlim(-width, width)

        # handle adding and removing names
        self.clearnames()
        for thiskey in self.pops.keys():
            self.namestars(thiskey)

        nudge = 0.05 * self.stretch(self.outer)
        for d in self.circlelabels.keys():
            self.circlelabels[d].set_visible(d > 0.2 * distance)
            self.circlelabels[d].set_x((nudge + self.stretch(d)) * np.cos(angle))
            self.circlelabels[d].set_y((nudge + self.stretch(d)) * np.sin(angle))

    def label(self, key):
        return (
            self.title.replace(" ", "")
            + "_"
            + key.lower()
            + "_{0:04.0f}pc".format(self.outer)
            + ".pdf"
        )

    def namestars(self, key):
        if True:
            old = self.key + ""
            self.point_at(key)

            xlim = self.ax.get_xlim()
            ylim = self.ax.get_ylim()
            onplot = (
                (self.x > np.min(xlim))
                * (self.x < np.max(xlim))
                * (self.y > np.min(ylim))
                * (self.y < np.max(ylim))
            )
            nottooclose = self.pop.distance > self.outer * 0.2

            nottoofar = self.pop.distance <= self.maxlabeldistance
            # print np.sum(onplot*nottooclose*nottoofar)

            # KLUDGE!!!!!
            if key.lower() == "tess":
                tolabel = (nottooclose * onplot).nonzero()[0]
            else:
                tolabel = (nottooclose * onplot * nottoofar).nonzero()[0]
            # tolabel = self.pop.find('WASP94Ab')
            if np.size(tolabel) > 1:
                tolabel = tolabel[np.unique(self.x[tolabel], return_index=True)[1]]

            for c in tolabel:
                try:
                    assert self.pop.alpha == 0
                    downwardnudge = ""
                except (AttributeError, AssertionError):
                    downwardnudge = "\n\n"

                text = downwardnudge + r"{}".format(self.pop.name[c].replace(" ", ""))
                self.labeled[c] = plt.text(
                    self.x[c],
                    self.y[c],
                    text,
                    color=self.pop.color,
                    alpha=self.pop.alpha,
                    va="center",
                    ha="center",
                    weight="bold",
                    size=self.planetfontsize,
                    zorder=self.pop.zorder,
                )
            self.point_at(old)

    def clearnames(self):
        """Remove system names that have been added to the plot."""
        for k in self.labeled:
            try:
                self.labeled[k].remove()
            except ValueError:
                pass

    def plot(self, key, labels=False):
        self.point_at(key)
        try:
            self.ax
        except:
            self.setup()
            self.circles(self.circlegrid)
        kw = self.kw()
        kw["facecolors"] = self.pop.color
        kw["edgecolors"] = "none"
        try:
            self.signature
        except:
            date = datetime.datetime.now().strftime("%B %Y")
            self.signature = self.ax.text(
                0.02,
                0.02,
                "animation by Zach Berta-Thompson\ndata from NASA Exoplanet Archive ({})".format(
                    date
                ),
                transform=self.ax.transAxes,
                ha="left",
                va="bottom",
                alpha=0.5,
                size=6,
            )
        try:
            kw["alpha"] = self.pop.alpha
        except AttributeError:
            pass
        self.ax.scatter(self.x, self.y, **kw)
        self.ax.scatter(
            0, 0, marker="x", s=100, alpha=1, color="lightgray", linewidth=4, zorder=-90
        )

        # kludge!
        self.leg = plt.legend(
            loc="upper right",
            fontsize=10,
            framealpha=0.0,
            scatterpoints=1,
            markerscale=1.5,
            title="Transiting Exoplanets",
        )

    def highlight(self, indices, label="special!"):
        print(self.pop.standard[indices])
        kw = self.kw()
        kw["marker"] = "*"
        kw["edgecolors"] = self.pop.color.replace("black", "gray")
        kw["facecolors"] = "white"
        kw["s"] *= 8
        kw["alpha"] = 1
        kw["label"] = None
        kw["zorder"] = 20000
        handle = self.ax.scatter(self.x[indices], self.y[indices], **kw)
        # self.highlightleg = plt.legend(label, fontsize=10, framealpha=0, markerscale=1.5)
        self.speak("highlighting {0} points".format(len(self.x[indices])))
        t = self.pop.teq[indices]
        r = self.pop.radius[indices]
        # kic = self.pop.standard['kepid'][indices]
        for i in range(len(t)):
            self.ax.text(
                self.x[indices][i],
                self.y[indices][i],
                r"{0:.0f}K, {1:.1f}R$_\oplus$".format(t[i], r[i]) + "\n",
                size=8,
                ha="center",
                va="bottom",
                color=self.pop.color.replace("black", "gray"),
            )
            print(r[i], t[i])  # , self.pop.distance[indices][i]
        # (self.name == 'WASP-94A b').nonzero()[0]

    def toclock(self, hour):
        """convert an hour to angle (in radians) on the plot"""
        return (6 - hour) * 2 * np.pi / 24.0

    def stretch(self, d, style="equalarea"):
        """stretch radial coordinate for plotting"""
        if style == "equalarea":
            return np.array(d) ** 1.5
        assert False


def tessComparison():
    non = NonKepler()
    kep = Kepler()
    unc = UnconfirmedKepler()
    kep.standard = vstack([kep.standard, unc.standard])
    tes = TESS()
    non.color = "royalblue"
    non.zorder = -1
    kep.color = "black"
    kep.zorder = 0
    tes.color = "coral"
    tes.zorder = 1
    t = ThumbtackPlot(pops=dict(nonkepler=non, kepler=kep, tess=tes))
    return t
