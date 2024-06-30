import exoatlas as ea
from astropy.table import Table, vstack, QTable
import numpy as np
import matplotlib.pyplot as plt
from astropy.time import Time
import astropy.units as u
import datetime as dt
from datetime import date
from datetime import datetime
from astropy.coordinates import SkyCoord
from astroquery.ipac.nexsci.nasa_exoplanet_archive import NasaExoplanetArchive
from astroplan import (
    FixedTarget,
    EclipsingSystem,
    AtNightConstraint,
    is_event_observable,
    AltitudeConstraint,
    LocalTimeConstraint,
)
from astroplan.plots import plot_airmass
from astroplan import Observer, TimeConstraint
import pytz
from pytz import timezone


def find_observable_transits(
    population, observer, obs_start_time, *args, **kwargs
):  # constraints=None, obs_start_time = datetime.now(), n_transits=1):
    """
    Finds n observable transits for a population of 1 object starting from a given observation start time.

    Parameters
    -----------
    population: an exoatlas Population object
        Only accepts population objects that have a length of 1 (only 1 target)
    observer: an astroplan Observer object
        Where the observations are taking place
    obs_start_time: an astropy Time object
        Whne the observations are beginning

    Args
    -----
    constraints: observing contraints from astroplan
        default: observations occur after civil twilight
    n_transits: int
        default: 1 transit from given obs_start_time

    Returns
    --------
    transit_table: astropy table
        an astropy table with the population name, ra, dec, transit period, transit midpoint, transit duration, transit ingress time,transit midpoint time
        ,and transit egress time

    """
    constraints = kwargs.get("constraints", AtNightConstraint.twilight_civil())
    # obs_start_time = kwargs.get('obs_start_time', Time(datetime.now(tz=pytz.utc)))
    n_transits = kwargs.get("n_transits", 1)

    name = population.name
    ra = population.ra.to_value("deg")
    dec = population.dec.to_value("deg")
    period = population.period.to_value("day")
    midpoint = population.transit_midpoint.to_value("day")
    duration = population.transit_duration.to_value("day")

    for i in range(len(population)):
        if len(population) > 1:
            print(
                "Population length cannot exceed 1. If transits for more than 1 object are desired, use function 'find_all_observable_transits'"
            )
            break
        elif np.isnan(ra) == True:
            print("Object ra not found. Please use a valid ra.")
            break
        elif np.isnan(dec) == True:
            print("Object dec not found. Please use a valid dec.")
            break
        elif np.isnan(period) == True:
            print("Object period not found. Please use a valid period.")
            break
        elif np.isnan(midpoint) == True:
            print(
                "Object transit midpoint not found. Please use a valid transit midpoint."
            )
            break
        elif np.isnan(duration) == True:
            print(
                "Object transit duration not found. Please use a valid transit duration."
            )
            break
        else:
            pass

        c = SkyCoord(ra, dec, unit="deg")
        target = FixedTarget(coord=c, name=name)

        midtransit_time = Time(midpoint, format="jd")
        transit_duration = duration * u.day
        period = period * u.day

        transit_table = QTable(
            names=(
                "name",
                "ra",
                "dec",
                "period",
                "transit_midpoint",
                "transit_duration",
                "ingress_time",
                "midpoint_time",
                "egress_time",
            ),
            dtype=(
                population.name.dtype,
                population.ra.dtype,
                population.dec.dtype,
                population.period.dtype,
                population.transit_midpoint.dtype,
                population.transit_duration.dtype,
                "U",
                "U",
                "U",
            ),
        )

        es = EclipsingSystem(
            primary_eclipse_time=midtransit_time,
            orbital_period=period,
            duration=transit_duration,
        )

        transit_times = es.next_primary_eclipse_time(
            obs_start_time, n_eclipses=n_transits
        )  # creating the array of mid-transit times
        ing_eg_times = es.next_primary_ingress_egress_time(
            obs_start_time, n_eclipses=n_transits
        )  # array of ingress/egress times

        in_iso = []
        eg_iso = []
        for j in range(len(transit_times)):
            in_iso = np.append(
                in_iso, Time(ing_eg_times[j, 0], format="iso", scale="utc")
            )
            eg_iso = np.append(
                eg_iso, Time(ing_eg_times[j, 1], format="iso", scale="utc")
            )

        in_eg_arr = np.vstack((in_iso, eg_iso)).T
        obs_arr = is_event_observable(
            constraints=constraints,
            observer=observer,
            target=target,
            times=transit_times,
            times_ingress_egress=in_eg_arr,
        )

        midpoint_times = []
        ingress = []
        egress = []
        for j in range(n_transits):
            if obs_arr[:, j] == True:
                midpoint_times = np.append(midpoint_times, transit_times[j])
                ingress = np.append(ingress, in_iso[j])
                egress = np.append(egress, eg_iso[j])

        for i in range(len(ingress)):
            transit_table.add_row(
                vals=(
                    population.name,
                    ra,
                    dec,
                    period,
                    midpoint,
                    duration,
                    ingress[i].value,
                    midpoint_times[i].value,
                    egress[i].value,
                )
            )

    return transit_table


def find_all_observable_transits(population, observer, obs_start_time, *args, **kwargs):
    """
    Finds n observable transits for a population of objects starting from a given observation start time.

    Parameters
    -----------
    population: an exoatlas Population object
    observer: an astroplan Observer object
        Where the observations are taking place
    obs_start_time: an astropy Time object
        Whne the observations are beginning

    Args
    -----
    constraints: observing contraints from astroplan
        default: observations occur after civil twilight
    n_transits: int
        default: 1 transit from given obs_start_time

    Returns
    --------
    transit_table: astropy table
        an astropy table with the population name, ra, dec, transit period, transit midpoint, transit duration, transit ingress time,transit midpoint time
        ,and transit egress time

    """
    constraints = kwargs.get("constraints", AtNightConstraint.twilight_civil())
    # obs_start_time = kwargs.get('obs_start_time', datetime.now(tz=pytz.utc))
    n_transits = kwargs.get("n_transits", 1)

    transit_table = QTable(
        names=(
            "name",
            "ra",
            "dec",
            "period",
            "transit_midpoint",
            "transit_duration",
            "ingress_time",
            "midpoint_time",
            "egress_time",
        ),
        dtype=(
            population.name.dtype,
            population.ra.dtype,
            population.dec.dtype,
            population.period.dtype,
            population.transit_midpoint.dtype,
            population.transit_duration.dtype,
            "U",
            "U",
            "U",
        ),
    )

    for i in range(len(population)):
        row = find_observable_transits(
            population=population[i],
            observer=observer,
            constraints=constraints,
            obs_start_time=obs_start_time,
            n_transits=n_transits[i],
        )

        transit_table = vstack([transit_table, row])
    return transit_table


def make_airmass_plot(transit_table, observer, savefig=False):
    """
    creates airmass plots from a given transit table

    Parameters
    -----------
    transit_table: an astropy table
        Result of find_observable transits or find_all_observable transits
    observer: an astroplan Observer object
        Where the observations are taking place


    Returns
    --------
    airmass plots:
        Airmass plots for each transit event in the given transit table

    """
    for i in range(len(transit_table)):
        name = transit_table["name"][i]
        ra = transit_table["ra"][i]
        dec = transit_table["dec"][i]
        midpoint_times = transit_table["midpoint_time"][i]
        ingress = transit_table["ingress_time"][i]
        egress = transit_table["egress_time"][i]

        c = SkyCoord(ra, dec, unit="deg")
        target = FixedTarget(coord=c)

        plot_airmass(
            targets=target,
            observer=observer,
            time=Time(midpoint_times),
            brightness_shading=True,
            altitude_yaxis=True,
            use_local_tz=True,
        )
        plt.axvline(
            Time(midpoint_times).plot_date,
            color="orange",
            linestyle="--",
            label="mid-transit",
        )
        plt.axvline(
            Time(ingress).plot_date, color="green", linestyle="--", label="ingress"
        )
        plt.axvline(
            Time(egress).plot_date, color="green", linestyle="--", label="egress"
        )
        plt.title("%s Transit midpoint: %s" % (name, Time(midpoint_times).value))
        plt.legend()
        if savefig == True:
            plt.savefig("%s_%d.png" % (name[i], j), bbox_inches="tight")
        plt.show()
