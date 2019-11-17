# tools for accessing tables from the NASA Exoplanet Archive,
# either by downloading them downloading or loading local tables


# to-do:
# [] implement filtering by the "where" keyword to the archive

from ..imports import *

class Downloader(Talker):
    expiration = np.inf

    # anything special to know about reading this file format?
    readkw =  dict(delimiter='|',
                   fill_values=[('',np.nan), ('--', np.nan)])



    def get(self, remake=False, skip_update=False):
        '''
        Get the table, downloading it from online if necessary.
        If the table is older than some particular threshold,
        then ask the user whether or not they want to download.

        Parameters
        ----------
        remake : bool
            Should we definitely redownload the table?
        skip_update : bool
            Should we skip checking if the table's out of date?
        '''


        # if file doesn't exist, download it
        if not skip_update:
            remake = remake | check_if_needs_updating(self.path,
                                                      self.expiration)

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
    supported_tables = ['exoplanets', 'compositepars']

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




exoplanets = ExoplanetArchiveDownloader('exoplanets')
composite_exoplanets = ExoplanetArchiveDownloader('compositepars')


class MergedExoplanetArchiveDownloader(ExoplanetArchiveDownloader):
    def __init__(self):
        self.table = 'merged'

    def download_fresh(self):
        '''
        Download a brand new merged table from the Exoplanet Archive.
        This grabs both the `exoplanets` and `compositepars` tables,
        and merges them together into one massive table with lots
        of columns for options.
        '''

        self.speak(f'Creating a merged exoplanet table from the NASA Exoplanet Archive.')

        # load the individual tables
        e = exoplanets.get()
        c = composite_exoplanets.get()

        # join the two tables together on the planets' names
        self.speak('Joining the two Exoplanet Archive tables together.')
        self.speak(' (This may take a while. The tables are big!)')
        c.rename_column('fpl_name', 'pl_name')
        j = join(e, c, keys='pl_name', table_names=['exoplanets', 'composite'])

        # tidy up some of the column names
        for k in j.colnames:
            if '_exoplanets' in k:
                j.rename_column(k, k.replace('_exoplanets', ''))

        # write the merged table out to a file
        self.speak(f'Merge successful! Saving file to {self.path}.')
        self.speak(' (This may take a while. The table is big!)')
        j.write(self.path,
                format='ascii.fixed_width',
                bookend=False,
                delimiter='|',
                overwrite=True)

        self.speak('File saved.')

merged_exoplanets = MergedExoplanetArchiveDownloader()


class ExoFOPDownloader(Downloader):
    expiration = 0.0
    def __init__(self):
        self.url = 'https://exofop.ipac.caltech.edu/tess/download_toi.php?sort=toi&output=pipe'

    @property
    def path(self):
        '''
        Where should the local copy of this file be stored?
        '''
        return os.path.join(directories['data'], 'TOI-exofop.txt')

toi_exofop = ExoFOPDownloader()


class MergedTOIDownloader(ExoFOPDownloader):
    '''
    Download the TOIs from the exoatlas table, but also search the MAST
    archive to pull out extra parameters for each star from the TIC.
    '''

    @property
    def path(self):
        '''
        Where should the local copy of this file be stored?
        '''
        return os.path.join(directories['data'], 'TOI-merged.txt')

    def download_fresh(self):
        '''
        Download the TOIs from ExoFOP, then search for their entries in
        the TIC catalog on the MAST to download more data.
        '''

        self.speak(f'Creating a merged TOI table from ExoFOP and MAST.')

        # download the table of TOIs from the ExoFOP
        self.speak('Downloading TOIs from ExoFOP.')
        t = toi_exofop.get(remake=True)

        # download the TIC entries for these stars
        self.speak(f'Searching for {len(t)} stars in the TIC on the MAST.')

        # import Catalogs only when we need it
        # (otherwise, we'll need the internet to ever run exoatlas)
        from astroquery.mast import Catalogs
        tic_table = Catalogs.query_criteria(catalog="Tic", ID=np.unique(t['TIC ID']))

        # preface all the columns with TIC so we can keep them straight
        for k in tic_table.colnames:
            tic_table[k].name = f'TIC {k}'

        # make sure the indices are integers so we can join on them
        tic_table['TIC ID'] = np.array(tic_table['TIC ID']).astype(np.int)

        # join the two tables together
        self.speak('Joining the TOI table with data from the TIC.')
        withtic = join(t, tic_table, 'TIC ID', join_type='left')

        # write the merged table out to a file
        self.speak(f'Merge successful! Saving file to {self.path}.')
        withtic.write(self.path,
                        format='ascii.fixed_width',
                        bookend=False,
                        delimiter='|',
                        overwrite=True)

        self.speak('File saved.')


toi_merged = MergedTOIDownloader()
