# tools for accessing tables from the NASA Exoplanet Archive,
# either by downloading them downloading or loading local tables


# to-do:
# [] implement filtering by the "where" keyword to the archive

from ..imports import *

# how many days old can local data be, before we try to download it
maximum_age = 1.0

class Downloader(Talker):
    # anything special to know about reading this file format?
    readkw = dict(delimiter='|')

    def time_from_modified(self):
        '''
        How long ago was this file last modified?
        '''
        try:
            dt = Time.now().unix - os.path.getmtime(self.path)
            return dt/60/60/24
        except FileNotFoundError:
            return np.inf

    def get(self, remake=False, ask=True):
        '''
        Get the table, downloading it from online if necessary.
        If the table is older than some particular threshold,
        then ask the user whether or not they want to download.
        '''

        # how long ago was the data updated?
        dt = self.time_from_modified()

        # if file doesn't exist, download it
        if dt == np.inf:
            remake = True
        # if file is old, ask if we should remake it
        elif dt > maximum_age:
            self.speak(f'{self.path} is {dt:.3f} days old.')
            if ask:
                answer = 'y' in self.input('Should we download a new one? [y/N]').lower()
                remake = remake | answer

        # either download a fresh file, or load a local one
        if remake:
            self.download_fresh()
        else:
            self.speak(f'Loading local file from {self.path}')

        # read the actual file
        return ascii.read(self.path, **self.readkw)


    def download_fresh(self):
        '''
        Download a brand new table from the Exoplanet Archive.

        Parameters
        ----------
        table : str
            Which table should we download from the archive?

        Return
        '''

        self.speak(f'Attempting to freshly download data from \n{self.url}')

        # download the file from the URL (10-minute timeout)
        temporary_path = download_file(self.url, cache=False, timeout=600)

        # copy the file to its new location
        shutil.copyfile(temporary_path, self.path)

        self.speak(f'Download successful! Saved file to {self.path}')


class ExoplanetArchiveDownloader(Downloader):

    # define the base of all URLs for to access the archive API
    base = 'http://exoplanetarchive.ipac.caltech.edu/cgi-bin/nstedAPI/nph-nstedAPI?'

    # what format do we want for the downloaded table?
    format='bar-delimited'

    # what columns do we want? (probably easiest to go with everything)
    select = '*'

    # what tables are currently support through this mode of access?
    supported_tables = ['exoplanets', 'compositepars', ]

    def __init__(self, table='exoplanets'):#, where='*'):
        self.table = table
        #self.where = where

    @property
    def url(self):
        '''
        Define the download URL to acces the online table.
        '''

        # create the URL by stitching together keyword pairs
        url = (f'{self.base}'
               f'table={self.table}'
               f'&select={self.select}'
               f'&format={self.format}')

        return url

    @property
    def path(self):
        '''
        Where should the local copy of this file be stored?
        '''
        return os.path.join(directories['data'],
                            f'nea-{self.table}.txt')


class ExoFOPDownloader(Downloader):

    def __init__(self):
        self.url = 'https://exofop.ipac.caltech.edu/tess/download_toi.php?sort=toi&output=pipe'

    @property
    def path(self):
        '''
        Where should the local copy of this file be stored?
        '''
        return os.path.join(directories['data'], 'TOI.txt')


all_exoplanets = ExoplanetArchiveDownloader('exoplanets')
composite_exoplanets = ExoplanetArchiveDownloader('compositepars')
toi_exofop = ExoFOPDownloader()
