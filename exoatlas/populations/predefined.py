from .population import *


class PredefinedPopulation(Population):
    """
    Predefined Populations are derived from archival data tables,
    representing the basic building blocks for many subsamples.
    """

    label = "?"
    expiration = 10 * u.day

    def __init__(self, remake=False, standard=None, label=None, **plotkw):
        """
        Initialize a predefined population.

        This will first try to load a standardarized local file. If that doesn't work,
        it will download raw data  (or load a recent local cached copy) and then
        carefully create a standard table, which can be saved for future use.
        Behind the scenes, this initialization is effectively something like:

        #   .ingest_standardized_data
        #   ...or...
        #   ..ingest_raw_data (= .download_raw_data + .create_standardized_data)

        Parameters
        ----------
        remake : bool
            Whether we re-ingest this table from its raw ingredients.
        standard : astropy.table.Table
            A standardized table with which to initialize this population;
            this is mostly used for creating subsets of existing predefined
            populations without having to reinitialize a cumbersome full
            object each time.
        label : str
            A custom label to apply to this population.
        **plotkw : dict
            All other keywords are stored as plotting suggestions.
        """

        if standard is None:
            try:
                # try to load the standardized table
                assert remake == False
                standard = self.ingest_standardized_data()
            except (IOError, FileNotFoundError, AssertionError):
                # or create a new standardized table and save it
                standard = self.ingest_raw_data(remake=remake)

        assert (type(standard) == Table) or (type(standard) == Row)

        # initialize with a standard table
        Population.__init__(
            self, standard=standard, label=label or self.label, **plotkw
        )

    def ingest_raw_data(self, remake=None):
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
        raw = self.download_raw_data(remake=remake)

        # create a standardized table from the array
        standard = self.create_standardardized(raw)

        # save the standardized table
        self.save_standardized_data(standard)

        return standard

    def download_raw_data(self):
        raise NotImplementedError(
            """
        Yikes! The `.download_raw_data` method has not been defined
        for whatever object is trying to call it!
        """
        )

    def ingest_standardized_data(self):
        """
        Load a standardized population table.

        This will *try* to ingest a standardized population table from
        a local file, probably one created by `create_standardized_table`.
        It can/will fail with either an IOError or FileNotFoundError if
        that file doesn't exist, likely prompting the `__init__()` that
        probably called this to run `.create_standardized_data()` to
        create and save a new standard table.

        Returns
        -------
        standard : astropy.table.Table
            A table of planet properties,
            with a minimal set of columns.
        """

        # keywords for reading a standardized table
        read_kw = dict(format="ecsv", fill_values=[("", np.nan), ("--", np.nan)])

        standard = ascii.read(self.standardized_data_path, **read_kw)
        self.speak(f"Loaded standardized table from {self.standardized_data_path}")

        return standard
