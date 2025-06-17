from ...imports import *
import astroplan as ap
from astroplan.plots import plot_airmass

sbo = ap.Observer(
    longitude=-105.2630 * u.deg,
    latitude=40.00371 * u.deg,
    elevation=1653 * u.m,
    name="Sommers-Bausch Observatory",
    timezone="US/Mountain",
)


def altaz(self, where=sbo, when="tonight", **kw):
    """
    Altitude, Azimuth (astropy.coordinates.SkyCoord)

    This returns an `astropy.coordinates.SkyCoord` object
    representing the altitude and azimuth of the objects
    in the population from a particular observer's
    location on Earth (`where`) at a particular moment (`when`).

    Caveats:
    It's (currently) meaningless for Solar System objects.
    It doesn't propagate uncertainties.
    It doesn't include proper motions or parallax.

    Parameters
    ----------
    where : astroplan.Observer
        This must be an astroplan Observer object
        encoding the location (and timezone) of an
        Earth-bound observatory.
    when : str, astropy.time.Time
        This can be a Time object specifying the moment at
        which the altitude and azimuth should be calculated.
        Or, it can be the string "tonight" meaning to use
        the closest local midnight to the current time.
    """

    if when == "tonight":
        when = where.midnight(Time.now())
    coordinates = SkyCoord(ra=self.ra(), dec=self.dec())
    altaz = where.altaz(when, coordinates)
    return altaz


def airmass(self, where=sbo, when="tonight", **kw):
    """
    Airmass (unitless)

    Calculate the airmass of a target from a particular observer's
    location on Earth (`where`) at a particular moment (`when`).

    Parameters
    ----------
    where : astroplan.Observer
        This must be an astroplan Observer object
        encoding the location (and timezone) of an
        Earth-bound observatory.
    when : str, astropy.time.Time
        This can be a Time object specifying the moment at
        which the altitude and azimuth should be calculated.
        Or, it can be the string "tonight" meaning to use
        the closest local midnight to the current time.
    """

    secz = self.altaz(where=where, when=when).secz

    airmass = np.where(secz > 0, secz, np.nan)
    return airmass


def plot_airmass_for_transit(row, max_airmass=2.0, savefig=False):
    """
    Make an airmass plot for one transit.

    For each row of a table produced by `show_upcoming_transits`,
    this function produces a plot showing the target's airmass
    from a particular location, with the transit indicated.

    Parameters
    -----------
    max_airmass : float
        The worst airmass to show on the plot.
    savefig : bool
        Should plots be saved out automatically
        into a folder called "airmass-plots"?
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        plt.figure()

        # use astroplan to do most of the work
        plot_airmass(
            targets=row["target"],
            observer=row.meta["where"],
            time=row["midpoint"],
            brightness_shading=True,
            altitude_yaxis=True,
            max_airmass=max_airmass,
            style_kwargs=dict(color="black"),
        )

        # plot the transit as vertical lines
        plt.axvline(
            row["midpoint"].plot_date,
            color="black",
            linestyle="--",
            label="mid-transit",
        )
        ingress_kw = dict(color="black", linestyle=":")
        plt.axvline(row["ingress"].plot_date, **ingress_kw)
        plt.axvline(row["egress"].plot_date, **ingress_kw)

        # add title
        plt.title(
            f'{row["name"]}\n{row.meta["where"].name} | {row["midpoint"].iso} UTC'
        )

        # set xlimit to be exactly sunset to sunrise
        plt.xlim(
            row.meta["where"].sun_set_time(row["midpoint"], which="previous").plot_date,
            row.meta["where"].sun_rise_time(row["midpoint"], which="next").plot_date,
        )

    # save the figure, if desired
    if savefig:
        mkdir("airmass-plots")
        plt.savefig(
            f'airmass-plots/{clean(row["name"])}-{row["midpoint"].iso.split()[0]}.png',
            bbox_inches="tight",
        )


def show_upcoming_transits(
    self,
    where=sbo,
    when="now",
    window=3 * u.day,
    allow_partial_transits=False,
    min_altitude=30 * u.deg,
    max_altitude=90 * u.deg,
    visualize=True,
):
    """
    Calculate all observable transits for a population.

    Parameters
    ----------
    where : astroplan.Observer
        This must be an astroplan Observer object
        encoding the location (and timezone) of an
        Earth-bound observatory.
    when : str, astropy.time.Time
        Specify the start of the window to start
        considering to search for transits, either
        as an astropy Time or as the string "now",
        meaning starting tonight for whenver you
        call the function.
    window : astropy.units.Quantity
        How long of a window should we consider for
        searching for potentially observable transits?
    allow_partial_transits : bool
        If True, only the transit midpoint must
        meet the constraints to be observable.
        If False, ingress + midpoint + egress must
        all together meet the constraints.
    min_altitude : astropy.units.Quantity
        The lowest altitude to consider observable.
        (30 degrees corresponds to airmass 2)
    max_altitude : astropy.units.Quantity
        The highest altitude to consider observable.
        (some telescopes can't observe at zenith!)
    visualize : bool
        Should we automatically produce an airmass
        plot for each observable transit event?
    """

    # construct skeleton table and a place to store transit events
    initial_planning_table = self.create_table(
        [
            "name",
            "ra",
            "dec",
            "magnitude_gaia",
            "transit_depth",
            "stellar_radius",
            "distance",
        ]
    )
    tables_for_individual_planets = []

    # define the start of the observing window
    if when == "now":
        when = Time.now()
    start_of_observing_window = where.midnight(when) - 0.5 * u.day

    # (keep track of if we need to warn about bad data)
    need_to_warn = True

    # there are some ~1s timing warnings we can probably ignore
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        # lopp through planets
        for i in range(len(self)):

            # set the target location on the sky
            target = ap.FixedTarget(
                coord=SkyCoord(ra=self.ra()[i], dec=self.dec()[i]), name=self.name()[i]
            )
            try:
                # set the orbital parameters
                es = ap.EclipsingSystem(
                    primary_eclipse_time=Time(self.transit_midpoint()[i], format="jd"),
                    orbital_period=self.period()[i],
                    duration=self.transit_duration()[i],
                )

                # calculate enough upcoming transits to fill dt
                n_transits = int(np.ceil(window / self.period()[i]))

                # (n_transits,) Time array
                midtransit_times = es.next_primary_eclipse_time(
                    start_of_observing_window, n_eclipses=n_transits
                )
                # (n_transits, 2) Time array
                ingress_and_egress_times = es.next_primary_ingress_egress_time(
                    when, n_eclipses=n_transits
                )
            except ValueError:
                # warn us if a planet has bad planning data
                if need_to_warn:
                    print("We could not predict transits for these planets:\n")
                    need_to_warn = False
                print(
                    f"☹️ {self.name()[i]:>20} | T0={self.transit_midpoint()[i]:>20} | P={self.period()[i]:10.3f} | dt={self.transit_duration()[i]:10.3f} ☹️"
                )
                continue

            # decide whether ingress/egress must be observable too
            if allow_partial_transits:
                times_ingress_egress = None
            else:
                times_ingress_egress = ingress_and_egress_times

            # calculate which of all possible events are observable
            ok = ap.is_event_observable(
                constraints=[
                    ap.AtNightConstraint(max_solar_altitude=-12 * u.deg),
                    ap.AltitudeConstraint(min=min_altitude, max=max_altitude),
                ],
                observer=where,
                target=target,
                times=midtransit_times,
                times_ingress_egress=times_ingress_egress,
            ).flatten()

            # if there are transits, record them!
            N = sum(ok)
            if N > 0:
                table_of_transits_for_this_planet = Table(
                    dict(
                        name=[self.name()[i]] * N,
                        ingress=ingress_and_egress_times[ok, 0],
                        midpoint=midtransit_times[ok],
                        egress=ingress_and_egress_times[ok, 1],
                        target=[target] * N,
                    )
                )
                tables_for_individual_planets.append(table_of_transits_for_this_planet)

    if len(tables_for_individual_planets) == 0:
        print("📭 No transits were found. Sorry! 📭")
        return None

    # stack and join individual planets into one giant table with metadata
    table_of_all_transits = vstack(tables_for_individual_planets)
    transit_planning_table = join(table_of_all_transits, initial_planning_table)
    transit_planning_table.meta["where"] = where
    transit_planning_table.meta["when"] = when
    transit_planning_table.meta["window"] = window

    if visualize:
        for row in transit_planning_table:
            plot_airmass_for_transit(row, savefig=True)

    return transit_planning_table
