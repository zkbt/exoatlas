from ...imports import *
from ..maps import ErrorMap, BubbleMap
import arviz as az


class Shoreline:
    """
    This `Shoreline` object handles all model
    calculations related to the cosmic shoreline
    from Berta-Thompson et al. (2025).
    """

    var_names = ["log_f_0", "p", "q", "ln_w"]

    def __init__(self, posterior=None, **kw):
        """
        Initialize a Shoreline object.

        Parameters
        ----------
        posterior : arviz.InferenceData, string, None
            The samples from the posterior probability distribution.
            If
        """
        if isinstance(posterior, az.InferenceData):
            self.posterior = posterior
        else:
            # download from zenodo and load it in with arviz
            raise NotImplementedError

        self.summary = az.summary(self.posterior, kind="all", stat_focus="median")

    def __repr__(self):
        return f"<🏝️{self.var_names}>"

    def best_parameters(self):
        """
        Return one "best" set of parameters.
        """
        best_parameters = self.summary["median"]
        return best_parameters

    def sampled_parameters(self, N_samples=100):
        """
        Return samples from the parameters.
        """
        df = self.posterior.to_dataframe(var_names=self.var_names)
        sampled_parameters = df[:: int(len(df) / N_samples)]
        return sampled_parameters

    def log_f_shoreline(self, log_f_0=0.0, p=0.0, q=0.0, log_v=0, log_L=0):
        """
        Calculate the location of the cosmic shoreline.

        If the parameters are scalars, the independent variables
        can be arrays (for example, many data points at once).

        If the independent variables are scalar, the parameters
        can be arrays (for example, many posterior samples).

        Parameters
        ----------
        log_f_0 : float, Quantity, array
            (model parameter), log10(flux normalization)
        p : float, Quantity, array
            (model parameter), escape velocity slope
        q : float, Quantity, array
            (model parameter), stellar luminosity slope
        log_v : float, Quantity, array
            (independent variable), log10(escape velocity relative to Earth)
        log_L : float, Quantity, array
            (independent variable), log10(stellar luminosity relative to Sun)
        """
        return log_f_0 + p * log_v + q * log_L

    def probability_of_atmosphere(
        self, log_f_0=1.0, p=4.0, q=0.0, ln_w=0, log_v=0, log_L=0, log_f=0
    ):
        """
        Calculate the probability a planet has an atmosphere.

        If the parameters are scalars, the independent variables
        can be arrays (for example, many data points at once).

        If the independent variables are scalar, the parameters
        can be arrays (for example, many posterior samples).

        Parameters
        ----------
        log_f_0 : float, Quantity, array
            (model parameter), log10(flux normalization)
        p : float, Quantity, array
            (model parameter), escape velocity slope
        q : float, Quantity, array
            (model parameter), stellar luminosity slope
        ln_w : float, Quantity, array
            (model parameter), ln(intrinsic width of the shoreline)
        log_v : float, Quantity, array
            (independent variable), log10(escape velocity relative to Earth)
        log_L : float, Quantity, array
            (independent variable), log10(stellar luminosity relative to Sun)
        log_L : float, Quantity, array
            (independent variable), log10(bolometric flux relative to Earth)
        """
        distance_from_shoreline = log_f - self.log_f_shoreline(
            log_f_0=log_f_0, p=p, q=q, log_v=log_v, log_L=log_L
        )
        width_of_shoreline = np.exp(ln_w)
        return 1 / (1 + np.exp(distance_from_shoreline / width_of_shoreline))


class ShorelineErrorMap(ErrorMap, Shoreline):
    """
    Define a cosmic shoreline map, for looking at
    some slice of the 3D cosmic shoreline.

    It contains the powers both of an ErrorMap
    and of a Shoreline with model posteriors.
    """

    def __init__(self, **kw):
        super().__init__(**kw)
        Shoreline.__init__(self, **kw)


class ShorelineBubbleMap(BubbleMap, Shoreline):
    """
    Define a cosmic shoreline map, for looking at
    some slice of the 3D cosmic shoreline.

    It contains the powers both of an BubbleMap
    and of a Shoreline with model posteriors.
    """

    def __init__(self, **kw):
        super().__init__(**kw)
        Shoreline.__init__(self, **kw)
