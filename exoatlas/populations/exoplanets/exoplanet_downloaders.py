"""
Define tool for downloading data from the NASA Exoplanet Archive.
"""

from ..downloader import *


class ExoplanetArchiveDownloader(Downloader):
    """
    ... to download Planetary Systems tables from the NASA Exoplanet Archive
    """

    # class-specific URL to access the archive API
    base = "https://exoplanetarchive.ipac.caltech.edu/TAP/sync?"

    # class-specific keywords to aid reading raw downloaded files into tables
    read_kw = dict(delimiter=",", fill_values=[("", np.nan), ("--", np.nan)])

    # class-specific names of supported tables
    supported_tables = ["ps", "pscomppars"]

    def __init__(
        self, which_table="ps", format="csv", select="select+*", planet="", **kw
    ):
        """
        Initialize an Exoplanet Archive downloader.

        Parameters
        ----------
        which_table : str
            The name of the table to be downloaded.
            - 'ps' is Planetary Systems with one line per reference
            - 'pscomppars' is Planetary Systems Composite Parameters with one line per planet
        format : str
            The format of the table to download.
        select : str
            The "select" string for the TAP query.
            This could be modified to select only a subset of columns.
        planet : str
            The name of planet to query, if desiring only a single planet.
            This is mostly used for testing the downloader with small datasets.
        **kw : dict
            Other keywords (expiration, timeout) are passed to the Downloader.
        """

        # keep track of what should be downloaded
        self.select = select
        assert which_table in self.supported_tables
        self.which_table = which_table
        self.planet = planet
        self.format = format
        Downloader.__init__(self, **kw)

    @property
    def url(self):
        """
        Define the download URL to acces the online table.
        """

        if self.planet is None:
            where_string = ""
        else:
            where_string = f"+where+pl_name+like+'{self.planet}'".replace(" ", "%20")

        # create the URL by stitching together keyword pairs
        url = (
            f"{self.base}"
            f"query={self.select}+"
            f"from+{self.which_table}"
            f"{where_string}"
            f"&format={self.format}"
        )

        return url

    @property
    def path(self):
        """
        Define a file path where the local copy of this file should be stored.
        """
        if self.planet == "":
            planet_string = ""
        else:
            planet_string = "-" + clean(self.planet)
        return os.path.join(
            directories["data"], f"nea-{self.which_table}{planet_string}.txt"
        )


planetary_systems_downloader = ExoplanetArchiveDownloader("ps")
composite_planetary_systems_downloader = ExoplanetArchiveDownloader("pscomppars")
