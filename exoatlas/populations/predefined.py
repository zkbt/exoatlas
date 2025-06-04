from .population_core import *


class PredefinedPopulation(Population):
    """
    Predefined Populations are derived from archival data tables,
    representing the basic building blocks for many subsamples.
    """

    label = "?"
    _expiration = 10 * u.day

    def __init__(
        self, remake=False, standard=None, label=None, verbose=False, **plotkw
    ):
        """
        Initialize a predefined population.

        This will first try to load a standardarized local file. If that doesn't work,
        it will download raw data  (or load a recent local cached copy) and then
        carefully create a standard table, which can be saved for future use.
        Behind the scenes, this initialization is effectively something like:

        #   ._ingest_standardized_data
        #   ...or...
        #   .._ingest_raw_data (= ._download_raw_data + ._create_standardized_data)

        Parameters
        ----------
        remake : bool
            Whether we re-ingest this table from its raw ingredients.
        standard : astropy.table.QTable
            A standardized table with which to initialize this population;
            this is mostly used for creating subsets of existing predefined
            populations without having to reinitialize a cumbersome full
            object each time.
        label : str
            A custom label to apply to this population.
        verbose : bool
            Say where the data were loaded from?
        **plotkw : dict
            All other keywords are stored as plotting suggestions.
        """

        self.verbose = verbose
        if standard is None:
            try:
                # try to load the standardized table
                assert remake == False
                standard = self._ingest_standardized_data()
            except (IOError, FileNotFoundError, AssertionError):
                # or create a new standardized table and save it
                standard = self._ingest_raw_data(remake=remake)

        assert type(standard) in [QTable, Table, Row]

        # initialize with a standard table
        Population.__init__(
            self,
            standard=standard,
            label=label or self.label,
            **plotkw,
        )

    def _ingest_raw_data(self, remake=None):
        """
        Ingest a new population table of arbitrary format,
        and then standardize it, using the tools defined in
        inherited population classes.

        Parameters
        ----------
        remake : bool
            Whether we should (definitely) re-download the data.
                If True, download data no matter what.
                If None, download data only if (there is none) or (it's expired and the user agrees).
                If False, download data only if there is none.
        """

        # load the raw table
        raw = self._download_raw_data(remake=remake)

        # create a standardized table from the array
        standard = self._create_standardized(raw)

        # save the standardized table as an ascii table for humans to read
        standard.write(
            self._standardized_data_path, format="ascii.ecsv", overwrite=True
        )
        print(f"Saved a standardized text table to {self._standardized_data_path}")

        return standard

    def _download_raw_data(self, remake=False):
        """
        Load raw data, possibly downloading them if needed.

        Parameters
        ----------
        remake : bool
            Should the raw tables be re-made/re-downloaded,
            even if recent local ones exist?

        Returns
        -------
        raw : astropy.table.QTable
            A raw, untrimmed, unstandardized table.
        """

        #
        self._raw = self._downloader.get_table(remake=remake)

        # return both, so they can both be turned into standardized tables
        return self._raw

    def _ingest_standardized_data(self):
        """
        Load a standardized population table.

        This will *try* to ingest a standardized population table from
        a local file, probably one created by `_create_standardized_table`.
        It can/will fail with either an IOError or FileNotFoundError if
        that file doesn't exist, likely prompting the `__init__()` that
        probably called this to run `._create_standardized_data()` to
        create and save a new standard table.

        Returns
        -------
        standard : astropy.table.QTable
            A table of planet properties,
            with a minimal set of columns.
        """

        # keywords for reading a standardized table
        read_kw = dict(format="ecsv", fill_values=[("", np.nan), ("--", np.nan)])

        standard = ascii.read(self._standardized_data_path, **read_kw)
        if self.verbose:
            print(f"Loaded standardized table from {self._standardized_data_path}")
        return standard
