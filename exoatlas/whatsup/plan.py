from ..imports import *
from .block import Block
from .observatory import Observatory
from .transit import Transit
from matplotlib import animation


# a list of colors to give transits (will repeat)
colors = [  'salmon',
            'firebrick',
            'hotpink',
            'mediumvioletred',
            'coral',
            'orangered',
            'darkorange',
            'plum',
            'fuchsia',
            'darkviolet',
            'purple',
            'indigo',
            'slateblue',
            'limegreen',
            'seagreen',
            'darkgreen',
            'olive',
            'teal',
            'steelblue',
            'royalblue',
            'sandybrown',
            'chocolate',
            'saddlebrown',
            'maroon']

class Plan(Talker):
    '''a Plan object to keep track of potentially observable transits'''

    def __init__(self,
                    pop=None,
                    start=None, finish=None,
                    observatory='SBO',
                    maxairmass=2.5,
                    maxsun=-12.0,
                    buffer=1.0,
                    directory='upcoming-transits/',
                    name='default'):

        '''initialize, setting the observatory, block, and population'''

        # initialize the talking setup
        Talker.__init__(self)

        # set the maximum airmass and the maximum sun altitude
        self.maxairmass = maxairmass
        self.maxsun = maxsun
        self.buffer = buffer

        # set up a directory to store results in
        self.directory = directory
        mkdir(self.directory)

        # set up our observatory
        self.observatory = Observatory(observatory, plan=self)

        # set up the time range this should cover
        self.block = Block(start=start, finish=finish, plan=self)

        # load the population of planets
        self.population = pop

        # tidy up that population
        self.tidy_population()

        # print out a summary of this plan
        self.describe_inputs()

    def describe_inputs(self):
        '''print text to terminal describing this observing plan'''

        self.speak('Observatory is {}'.format(self.observatory))
        self.speak('Time range spans {} to {}'.format(self.block.start, self.block.finish))
        self.speak('Population contains {} objects'.format(len(self.population)))

    def tidy_population(self):
        '''load a population of exoplanets, defaulting to the list
        of confirmed transiting planets from the NASA Exoplanet Archive'''

        # remove planets that are impossible to see from this latitude
        zenith = np.abs(self.population.dec - self.observatory.latitude)
        possible = zenith < 60.0*u.deg
        self.speak('{}/{} targets are visible from {} latitude'.format(
                            np.sum(possible), len(possible),
                            self.observatory.latitude))

        # trim the population
        self.population = self.population[possible]

        self.speak('the trimmed population contains {} objects'.format(
                        len(self.population.standard)))

        #KLUDGE
        bad_duration = np.isfinite(self.population.transit_duration) == False

        #print('These planets have bad durations!')
        #print(self.population[bad_duration].name)

        self.population = self.population[bad_duration == False]
        # add some colors to the population
        N = len(self.population)
        c = np.array(colors)[np.arange(N) % len(colors)]
        self.population.standard['observing_color'] = c
        self.population._pithy = True

    def find_planets_transiting_on_night(self, midnight):
        '''identify all the transits for all the planets'''

        now = midnight
        pop = self.population


        # pull out the period and epoch of all the planets
        P = pop.period.to('day')
        T0 = pop.transit_epoch.to('day')
        coords = SkyCoord(ra=pop.ra, dec=pop.dec)

        # find the closest exact transit epochs
        pop.standard['closest_epoch_number'] = np.round((now - T0)/P)
        pop.standard['closest_epoch_day'] = pop.closest_epoch_number*P + T0

        # ask which planets have transits happening tonight/today
        dt = now - pop.closest_epoch_day
        is_happening_tonight = np.abs(dt) < 0.5*u.day
        sum(is_happening_tonight), len(is_happening_tonight)

        # pull out only the times where a transit is happening within 0.5 days
        times_happening =  Time(pop.closest_epoch_day[is_happening_tonight], format='jd')
        buffer = self.buffer*pop.transit_duration[is_happening_tonight]
        coords_happening = coords[is_happening_tonight]

        # figure out if the targets are up
        min_alt = 30*u.deg
        up_happening = self.observatory.altaz(coords_happening, times_happening).alt > min_alt
        up_happening *= self.observatory.altaz(coords_happening, times_happening + buffer).alt > min_alt
        up_happening *= self.observatory.altaz(coords_happening, times_happening - buffer).alt > min_alt

        # figure out if the sun is down
        max_sun = -12*u.deg
        dark_happening = self.observatory.sun(times_happening).alt < max_sun
        dark_happening *= self.observatory.sun(times_happening+buffer).alt < max_sun
        dark_happening *= self.observatory.sun(times_happening-buffer).alt < max_sun

        # populate a mask of which transits we should both thinking about
        should_consider = True
        should_consider *= is_happening_tonight
        should_consider[is_happening_tonight] *= up_happening
        should_consider[is_happening_tonight] *= dark_happening

        planets_tonight = pop[should_consider]
        return planets_tonight

    def plot_transits_on_night(self, midnight):


        # plot black bars to mask out when the Sun is up
        self.observatory.plot_empty_night(midnight)

        planets = self.find_planets_transiting_on_night(midnight)
        for p in planets:
            t = Transit(p, p.closest_epoch_number, observatory=self.observatory)
            t.plot(color=p.observing_color[0])

        plt.ylim(2.3, 0.9)

    def plot_all_transits(self):
        '''
        Plot all observable transits that fall within this night.
        '''

        for midnight in tqdm(self.block.midnights):

            self.plot_transits_on_night(midnight)
            date = Time(midnight - self.observatory.standardzone -0.5*u.day,
                        format='jd').iso[:10]
            plt.title(f'Transits Observable from {self.observatory.name} on {date}')

            filename = os.path.join(self.directory,
                                    f'{date}.pdf')
            plt.savefig(filename)
            plt.clf()

    def print_transits(self, filename='upcoming_transits.txt'):
        '''print all the observable transits'''
        with open(filename, 'w') as f:
            for p in self.planets:
                p.speak(p.name)
                for t in p.transits:
                    s = t.details()
                    f.write(s + '\n')
                #(t.posttransit - t.pretransit).value*24))
                #self.observatory.sun(t.pretransit).alt.deg[0],
                #self.observatory.sun(t.posttransit).alt.deg[0]))


    def movieTransits(self, filename=None, fps=10, bitrate=10000):
        '''
        Make a movie of the whole block,
        showing all the transits for your population.
        '''

        self.speak('plotting transits')
        self.find_transits()
        plt.ioff()

        # set up the plotting window (one night at a time)
        self.figure = plt.figure('upcoming transits', figsize=(15, 6), dpi=200)
        self.gs = plt.matplotlib.gridspec.GridSpec(2,1, hspace=0, wspace=0, height_ratios=[0.3, 1.0], bottom=0.35)
        self.ax = {}

        # an ax to include the names of the planets
        self.ax['transits'] = plt.subplot(self.gs[0])
        plt.setp(self.ax['transits'].get_xticklabels(), visible=False)
        plt.setp(self.ax['transits'].get_yticklabels(), visible=False)
        self.ax['transits'].plot_date([],[])
        plt.setp(self.ax['transits'].get_xaxis(), visible=False)
        plt.setp(self.ax['transits'].get_yaxis(), visible=False)
        self.ax['transits'].set_title('Observability from {0}'.format(self.observatory.name, self.block.name))

        # an ax to include the airmass curves
        self.ax['airmass'] = plt.subplot(self.gs[1],
            sharex=self.ax['transits'])
        self.ax['airmass'].set_ylabel('Airmass')

        # plot all the individual transits for all the planets
        for i in range(len(self.planets)):
            self.planets[i].plotTransits(y=i)

        # nudge the limits on the transit labels
        # self.ax['transits'].set_xlim(self.block.start.plot_date, self.block.finish.plot_date)
        self.ax['transits'].set_ylim(-0.1*len(self.planets), 1.2*len(self.planets))

        # plot black bars to mask out when the Sun is up
        for k in self.ax.keys():
            self.observatory.plotSun(jd=self.block.times, ax=self.ax[k])

        # set the formatting on the time label
        fmt =  plt.matplotlib.dates.DateFormatter('%Y-%m-%d\n%H:%M UT')
        self.ax['airmass'].xaxis.set_major_formatter(fmt)
        f = plt.gcf()
        f.autofmt_xdate(rotation=90)
        plt.tight_layout()


        self.speak('making a movie of transits')

        # initialize the animator
        metadata = dict(title='Observability from {0}'.format(self.observatory.name),
                        artist='whatsup')
        #self.writer = animation.FFMpegWriter(fps=fps, metadata=metadata, bitrate=bitrate)
        self.writer = animation.PillowWriter(fps=fps, metadata=metadata, bitrate=bitrate)
        if filename is None:
            filename = self.directory + 'transits-{}planets.mp4'.format(len(self.planets))
        self.speak('trying to save movie to {0}'.format(filename))
        with self.writer.saving(self.figure, filename, self.figure.get_dpi()):
            # loop over exposures
            window = 7/24.0

            # approximate kludge; no daylight savings; need to incorporate time zones
            offset = self.observatory.standardzone.to('day').value

            for m in tqdm(self.block.midnights):
                xlim = self.ax['airmass'].get_xlim()
                center = Time(m, format='jd').plot_date + offset
                self.ax['airmass'].set_xlim((center - window), (center + window))

                nightof = Time(m - 0.5, format='jd').iso.split()[0]
                self.ax['transits'].set_title('Observability from {} on {}'.format(self.observatory.name, nightof))
                self.writer.grab_frame()
