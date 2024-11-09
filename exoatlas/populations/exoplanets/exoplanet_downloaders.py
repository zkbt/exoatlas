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
    read_kw = dict(
        format="ascii", delimiter=",", fill_values=[("", np.nan), ("--", np.nan)]
    )

    # class-specific names of supported tables
    supported_tables = ["ps", "pscomppars"]

    def __init__(self, which_table="ps", format="csv", N=np.inf, planet="", **kw):
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
        N : int
            The maximum number of rows to return.
            This is mostly used for testing the downloader with small datasets.
        planet : str
            The name of planet to query, if desiring only a single planet.
            This is mostly used for testing the downloader with small datasets.
        **kw : dict
            Other keywords (expiration, timeout) are passed to the Downloader.
        """

        # keep track of what should be downloaded
        assert which_table in self.supported_tables
        self.which_table = which_table
        self.format = format
        self.N = N
        self.planet = planet
        Downloader.__init__(self, **kw)

    @property
    def top_string(self):
        # limit to a certain number of rows
        if (self.N == None) or (np.isfinite(self.N) == False):
            return ""
        else:
            return f"top+{self.N}+"

    @property
    def where_string(self):
        if (self.planet is "") or (self.planet is None):
            return ""
        else:
            return f"where+pl_name+like+'{self.planet}'".replace(" ", "%20")

    @property
    def from_string(self):
        return f"from+{self.which_table}+"

    @property
    def url(self):
        """
        Define the download URL to acces the online table.
        """

        # create the URL by stitching together keyword pairs
        url = (
            f"{self.base}"
            f"query=select+{self.top_string}*+"
            f"{self.from_string}"
            f"{self.where_string}"
            f"&format={self.format}"
        )

        return url

    @property
    def path(self):
        """
        Define a file path where the local copy of this file should be stored.
        """
        label = f"{self.top_string}{self.from_string}{self.where_string}".strip("+")
        return os.path.join(directories["data"], f"nasa-exoplanet-archive-{label}.txt")
