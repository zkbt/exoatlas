from .downloaders import *

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
    Download the TOIs from the exopop table, but also search the MAST
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
