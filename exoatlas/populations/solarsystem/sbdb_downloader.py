"""
Define tool for downloading data from the NASA Exoplanet Archive.
"""

from ..downloader import *
import requests, json


class SBDBDownloader(Downloader):
    """
    ... to download SBDB data from JPL Solar System Dynamics.
    """

    # class-specific URL to access the archive API
    base = "https://ssd-api.jpl.nasa.gov/sbdb_query.api"

    def __init__(self, minimum_diameter=100 * u.km, **kw):
        """
        Initialize a Solar System small body downloader.

        Parameters
        ----------
        **kw : dict
            Other keywords (expiration, timeout) are passed to the Downloader.
        """

        # keep track of what should be downloaded
        self.minimum_diameter = minimum_diameter
        Downloader.__init__(self, **kw)

    def download_fresh(self):
        """
        Download a brand new table from the JPL/SBDB Query API.
        """

        print(
            f"""
        Attempting to freshly download data from
        {self.base}
        
        This may take a long time, and we're dreadully 
        sorry that we haven't coded in a clever way to 
        tell you how long. Archive tables might be up to
        a few 100MB, so however long it takes your computer
        to download that is probably a reasonable-ish amount 
        of time to wait before getting stressed out. 
        """
        )

        self._parameters = {}
        self._parameters["fields"] = (
            "full_name,name,class,e,a,q,i,om,w,ma,tp,per,H,G,diameter,GM,density,rot_per,albedo"
        )
        self._parameters["sb-cdata"] = (
            '{"AND":["diameter|GT|' + f"{self.minimum_diameter.to_value('km')}" + '"]}'
        )

        # download the file
        self._downloaded = requests.get(
            self.base, params=self._parameters, timeout=self.timeout.to_value("second")
        )
        self.url = self._downloaded.url
        print("Download successful. Processing into a table.")

        # process into an astropy Table
        d = json.loads(self._downloaded.text)
        self._raw_data = d
        # create a dictionary to organize data
        organized = {}
        for i, k in enumerate(d["fields"]):
            this_column = [x[i] for x in d["data"]]
            if k in ["full_name", "name", "class"]:
                organized[k] = np.array(this_column).astype(str)
            else:
                organized[k] = np.array(this_column).astype(float)

        self.table = Table(organized)
        self.table.meta["url"] = self.url
        self.table.meta["minimum_diameter"] = self.minimum_diameter
        self.table.write(self.path, overwrite=True)

        print(f"Download successful! Saved file to {self.path}")

    @property
    def path(self):
        """
        Define a file path where the local copy of this file should be stored.
        """
        label = f"{self.minimum_diameter}km"
        return os.path.join(directories["data"], f"jpl-ssd-sbdb-{label}.hdf5")
