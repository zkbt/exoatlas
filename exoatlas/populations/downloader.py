"""
Define a downloader tool to somewhat politely use a local file and/or 
download a big one from online, if it's older than some particular date.
"""

from ..imports import *


class Downloader(Talker):
    """
    ...to aid downloading large archive files, storing them locally,
    and redownloading them only when/if absolutely necessary.
    """

    # class-specific keywords to help reading files
    read_kw = dict()

    def __init__(
        self,
        expiration=1.0 * u.day,
        timeout=20 * u.minute,
    ):
        """
        Initialize a Downloader.

        Parameters
        ----------
        expiration : astropy.units.quantity.Quantity
            After how long should we consider the raw downloaded table expired,
            and therefore download a new one?
        timeout : astropy.units.quantity.Quantity
            What's the longest we should wait before giving up on the download?
        """

        self.expiration = expiration
        self.timeout = timeout

    def get(self, remake=False):
        """
        Get the table, downloading it from online if necessary.
        If the table is older than the downloader's threshold,
        then ask the user whether or not they want to download.

        Parameters
        ----------
        remake : bool
            Whether we should (definitely) re-download the data.
                If True, download data no matter what.
                If None, download data only if (there is none) or (it's expired and the user agrees).
                If False, download data only if there is none.
        """

        # if file doesn't exist, download it

        if remake == True:
            should_we_download = True
        if remake == None:
            should_we_download = check_if_needs_updating(self.path, self.expiration)
        elif remake == False:
            should_we_download = os.path.exists(self.path) == False

        # either download a fresh file, or load a local one
        if should_we_download:
            self.download_fresh()
        else:
            self.speak(f"Loading local file from {self.path}")

        # read the actual file
        table = ascii.read(self.path, **self.read_kw)

        # possibly add some new metadata to the table
        self.add_metadata(table)

        return table

    def download_fresh(self):
        """
        Download a brand new table from the Exoplanet Archive.
        """

        self.speak(
            f"""
        Attempting to freshly download data from
        {self.url}
        
        This may take a long time, and we're dreadully 
        sorry that we haven't coded in a clever way to 
        tell you how long. Archive tables might be up to
        a few 100MB, so however long it takes your computer
        to download that is probably a reasonable-ish amount 
        of time to wait before getting stressed out. 
        """
        )

        # download the file from the URL (20-minute timeout)
        temporary_path = download_file(
            self.url, cache=False, timeout=self.timeout.to_value("second")
        )

        # copy the file to its new location
        shutil.copyfile(temporary_path, self.path)

        self.speak(f"Download successful! Saved file to {self.path}")

    def add_metadata(self, table):
        pass
