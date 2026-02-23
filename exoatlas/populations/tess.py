from ..imports import *


def get_tess_sector_information(pop):
    """
    Get TESS sectors, cameras, CCDs, and pixels for a population.

    Parameters
    ----------
    pop : Population
        The population to test for observability with TESS.

    Returns
    -------
    table : Table
        An astropy table in which the first column `id`
        corresponds to the index of the population.
        Stars may appear 0 times (= not observed) or
        one or more times (= one row for each sector).
    """

    from tess_stars2px import tess_stars2px_function_entry

    ra = pop.ra().to_value("deg")
    dec = pop.dec().to_value("deg")
    id = np.arange(len(pop))
    (
        outID,
        outEclipLong,
        outEclipLat,
        outSec,
        outCam,
        outCcd,
        outColPix,
        outRowPix,
        scinfo,
    ) = tess_stars2px_function_entry(id, ra, dec)

    results = Table(
        dict(
            id=outID,
            ecliptic_longitude=outEclipLong,
            ecliptic_latitude=outEclipLat,
            tess_sector=outSec,
            tess_camera=outCam,
            tess_ccd=outCcd,
            tess_col=outColPix,
            tess_row=outRowPix,
        )
    )
    return results


def attach_tess_sectors(pop):
    """
    Create one population with a `.tess_sector` column
    that contains an array of TESS sectors containing a star.

    Parameters
    ----------
    pop : Population
        The population to test for observability with TESS.

    Returns
    -------
    pop_with_tess : Population
        A new Population, with TESS sectors attached.
    """

    # calculate the TESS sectors + positions of each target
    results = get_tess_sector_information(pop)

    # each star gets multiple sectors listed
    sectors = []
    for i in np.arange(len(pop)):
        is_this_star = results["id"] == i
        if np.sum(is_this_star) > 0:
            sectors_for_this_star = np.array(results["tess_sector"][is_this_star])
        else:
            sectors_for_this_star = np.array([])
        sectors.append(sectors_for_this_star)
    population_with_sectors = pop[:]
    population_with_sectors.add_column(name="tess_sector", data=sectors)
    return population_with_sectors


def split_into_tess_sectors(pop):
    """
    Create a dictionary of populations split by
    TESS sector.

    Parameters
    ----------
    pop : Population
        The population to test for observability with TESS.

    Returns
    -------
    pops : dict
        A dictionary of new Populations, split out
        according to their TESS sectors.
    """

    # calculate the TESS sectors + positions of each target
    results = get_tess_sector_information(pop)

    # one population per sector
    populations_per_sector = {}
    for sector in np.unique(results["tess_sector"]):
        s = int(sector)
        in_sector = results["tess_sector"] == s
        results_this_sector = results[in_sector]
        population_in_this_sector = pop[results_this_sector["id"]]
        population_in_this_sector.color = None
        population_in_this_sector.label = f"TESS S{s}"
        population_in_this_sector.s = 100
        for k in results_this_sector.colnames:
            population_in_this_sector.add_column(name=k, data=results_this_sector[k])
        populations_per_sector[s] = population_in_this_sector
    return populations_per_sector
