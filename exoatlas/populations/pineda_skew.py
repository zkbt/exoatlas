from ..imports import *

from scipy.stats.mstats import mquantiles
from scipy.optimize import fsolve
import scipy.stats as STS

from scipy.signal import medfilt
from astropy.table import Table


def confidenceInterval(values, interval=0.682689, ThSig=False):
    "For determining the median and central confidence interval of input values array"
    if ThSig:
        interval = 0.997300203937  # 3 sigma
    qL = (1.0 - interval) / 2.0
    qU = 1.0 - qL
    out = mquantiles(values, prob=[qL, 0.5, qU])
    return out


from scipy.stats import skewnorm, norm

gaussian_central_1sigma = np.diff(norm(loc=0, scale=1).cdf([-1, 1]))[0]


def plot_skewnormal(mu=0, sigma=1, alpha=0):
    fi, ax = plt.subplots(2, 1, sharex=True, constrained_layout=True, figsize=(7, 3))
    p = skewnorm(a=alpha, loc=mu, scale=sigma)
    x = np.linspace(mu - 5 * sigma, mu + 5 * sigma, 1000)
    skew_median = p.median()
    skew_lower, skew_upper = p.interval(gaussian_central_1sigma)
    skew_mode = x[np.argmax(p.pdf(x))]

    median_description = rf"""
    median = {skew_median:.2f}  
    median - {gaussian_central_1sigma/2:.2%} = {skew_median - skew_lower:.2f}  
    {gaussian_central_1sigma/2:.2%} - median = {skew_upper - skew_median:.2f}  
    upper/lower = {(skew_upper - skew_median)/(skew_median - skew_lower):.2f}  
    """

    mode_description = rf"""
    mode = {skew_mode:.2f}  
    mode - {gaussian_central_1sigma/2:.2%} = {skew_mode - skew_lower:.2f}  
    {gaussian_central_1sigma/2:.2%} - mode = {skew_upper - skew_mode:.2f}  
    upper/lower = {(skew_upper - skew_mode)/(skew_mode - skew_lower):.2f}  
    """

    plt.sca(ax[0])
    plt.plot(x, p.pdf(x))
    plt.text(
        0,
        1,
        median_description,
        ha="left",
        va="top",
        transform=ax[0].transAxes,
        color="royalblue",
    )
    plt.text(
        1,
        1,
        mode_description,
        ha="right",
        va="top",
        transform=ax[0].transAxes,
        color="purple",
    )

    plt.ylabel("p(x)")
    plt.title(f"\nmu = {mu:.2f}, sigma = {sigma:.2f}, alpha = {alpha:.2f}")
    plt.sca(ax[1])
    plt.plot(x, p.cdf(x))
    plt.ylabel(r"$\int_0^{x} p(x') dx$")

    for a in ax:
        plt.sca(a)
        plt.axvline(skew_median, color="royalblue")
        plt.axvline(skew_mode, color="purple")

        plt.axvspan(skew_lower, skew_upper, alpha=0.3, color="gray")


def estSkewNorm(xin, conf=(0.16, 0.5, 0.84), Guess=None, Mode="Med", Check=True):
    """
    Used for estimating the parameters of a Skew Normal distribution that matches input mode/median and given confidence interval.

    Numerically solves system of equations to give the 3 output parameters of corresponding skew normal. May require multiple guess to result in right solution. Output should give practically exact values, see `Check'.

    Returns the skew normal parameters (mu, sigma, alpha) --- see wikipedia on skew normal

    WARNING: Code is provided with minimal explanation, and verfication --- use with that in mind.

    ZKBT borrowed this amazing code from Sebastian Pineda's
    https://github.com/jspineda/stellarprop/
    on 13 November 2024. If you want to understand it more and/or use it, please see Appendix B
    of https://ui.adsabs.harvard.edu/abs/2021ApJ...918...40P/abstract, and cite that paper!

    >---

    Parameters
    ----------
    xin : array_like, tuple, or list, with 3 entries
        entries correspond to the data value matching the corresponding index in the `conf' key word. By default first entry indicates the lower bound for the central 68% confidence interval, second entry corresponds to median/mode value depending on keyword `Mode', and the third entry is the upper bound for the central 68% confidence interval.

    conf : array_like, tuple, or list, with 3 entries
        indicates the values that the skew normal cumulative probability distribution should give with input `xin'. By default, set to median and central 68% confidence interval. If Mode is `Peak' the median equation is replaced by one corresponding to peak of distribution.

    Guess : array_like, tuple or list, with 3 entries
        allows for user to input starting point for numerical equation solve. Default values are a general guess. If output does not converge, use Guess to change starting point for computation. May require iteration for adequete solution. Use `Check' to verify. If there is difficult, input parameters may not be well suited for approximation with skew normal.

    Mode : str one of ['Peak','Med2','Med','SF']
        Defines to set of equations used in defining the skew normal distribution. If 'Peak' system sets second entry to be the mode of skew normal instead of median. All others are for setting the median, but with slightly different numerical implementations. 'Peak' and 'Med' are the recommended modes.

    Check : boolean
        By default it is True. Used as verification on output solution, gives printed diagnostics as check. Outputs for converged solutions should be exact if fully successful.


    Returns
    -------

    out : array_like with 3 entries
        gives the (mu, sigma, alpha) parameters that define the skew normal distribution matching inputs

    Notes
    -----
    Printed warnings also given from scipy from fsolve to diagnose progress of numerical solution to system of equations

    Examples
    --------

    ## Note that here we use data from https://github.com/jspineda/stellarprop for illustration ; see also Pineda et al. 2021b

    >> trace = np.genfromtxt('../resources/massradius/fractional_01/chains.csv',delimiter=',',names=True)
    >> scatlb, scatmid, scatub = confidenceInterval(trace['scatter'])  # the scatter psoterior distribution is asymetric, these are typically the reported values in literature
    >> print([scatlb,scatmid,scatub])
    [0.02787424918238516, 0.0320051813038165, 0.03692976181631807]
    >> params = estSkewNorm( [scatlb, scatmid, scatub])
    Mode at [0.03121118]
    Median at 0.032005181304171265
    Result gives centeral 68.0% confidence interval: (0.027874249182851436, 0.03692976181636316)
    Peak at [0.03121118] - [0.00333693]  +  [0.00571858]
    >> print(params)
    [0.02771848 0.0065575  1.95731243]

    ## Note that Check outputs reported numerically match nearly exactly to inputs, these would be kinda off if iteration needed
    ## In this example alpha ~2, indicating positive skewness, peak (mode) is at 0.031, a little less than median at 0.032   -- see appendix of Pineda et al. 2021b


    """

    xl, x0, xu = xin
    cl, c0, cu = conf
    if Guess is not None:
        p0 = Guess
    else:
        p0 = (x0, (xu - xl) / 2.0, ((xu - x0) - (x0 - xl)) / ((xu - xl) / 2.0) * 3)

    ## if block used to toggle set of equations to solve using scipy fsolve
    if Mode == "Peak":
        # print("Setting Peak of Distribution")
        # KLUDGE!!!!!
        mu_guess = x0
        sigma_guess = np.sqrt((xu - x0) ** 2 + (x0 - xl) ** 2)
        alpha_guess = np.log((xu - x0) / (x0 - xl)) * 3
        p0 = mu_guess, sigma_guess, alpha_guess

        def eq_sys(p):
            mu, sigma, alpha = p
            t = (x0 - mu) / sigma
            return (
                STS.skewnorm.cdf(xl, alpha, mu, sigma) - conf[0],
                STS.skewnorm.cdf(xu, alpha, mu, sigma) - conf[2],
                alpha * STS.norm.pdf(alpha * t) - STS.norm.cdf(alpha * t) * t,
            )

    elif Mode == "Med2":
        # print("Setting Median of Distribution")

        def eq_sys(p):
            mu, sigma, alpha = p
            return np.power(STS.skewnorm.cdf(xin, alpha, mu, sigma) - np.array(conf), 2)

    elif Mode == "SF":
        # print("Setting Median of Distribution")

        def eq_sys(p):
            mu, sigma, alpha = p
            return STS.skewnorm.isf(1 - np.array(conf), alpha, mu, sigma) - np.array(
                xin
            )

    elif Mode == "Med":
        # print("Setting Median of Distribution")

        def eq_sys(p):
            mu, sigma, alpha = p
            return (
                STS.skewnorm.cdf(xl, alpha, mu, sigma) - cl,
                STS.skewnorm.cdf(x0, alpha, mu, sigma) - c0,
                STS.skewnorm.cdf(xu, alpha, mu, sigma) - cu,
            )

    out, infodict, ier, message = fsolve(eq_sys, p0, factor=0.1, full_output=True)
    ok = bool(ier)
    mu, sigma, alpha = out

    if Check:
        ff = lambda a: STS.norm.pdf(out[2] * a) * out[2] - a * STS.norm.cdf(a * out[2])
        tm = fsolve(ff, 0.2 * out[2])
        xm = tm * out[1] + out[0]
        print("Mode at {}".format(xm))
        print("Median at {}".format(STS.skewnorm.median(out[2], out[0], out[1])))
        print(
            "Result gives centeral {0}% confidence interval:".format(
                (conf[2] - conf[0]) * 100
            ),
            STS.skewnorm.interval(conf[2] - conf[0], out[2], out[0], out[1]),
        )
        print("Peak at {0} - {1}  +  {2} ".format(xm, xm - xl, xu - xm))

        # ZKBT adds a plot
        plot_skewnormal(mu, sigma, alpha)

    return mu, sigma, alpha, ok


def make_skewnormal_parameters_to_interpolate(N=1000):
    """
    Make a table to convert from sigma_upper/sigma_lower to skewnormal parameters.

    This function is/was used to generate the table in
    exoatlas/populations/data/skewnormal_parameters_for_interpolating_mode.ecsv
    that gets loaded whenever exoatlas is imported, to be
    able to quickly create skewnormal distributions for
    quantities with asymmetric uncertainties.

    """
    center = 0
    sigma_lower = 1
    sigmas_upper = np.logspace(-2, 2, N)
    tables = {}
    mylabels = {"Med": "median", "Peak": "mode"}
    for mode in ["Peak", "Med"]:
        iterations = 2

        x0 = center
        xu = center + sigmas_upper
        xl = center - sigma_lower
        mu_guess = x0 * np.ones_like(sigmas_upper)
        sigma_guess = np.sqrt((xu - x0) ** 2 + (x0 - xl) ** 2)
        if mode == "Peak":
            alpha_guess = np.log((xu - x0) / (x0 - xl)) * 3  # 100
        elif mode == "Med":
            alpha_guess = ((xu - x0) - (x0 - xl)) / (xu - xl) * 100
        guesses = zip(mu_guess, sigma_guess, alpha_guess)

        for i in range(iterations):
            parameters = np.array(
                [
                    estSkewNorm(
                        [center - sigma_lower, center, center + sigma_upper],
                        Mode=mode,
                        Guess=guess,
                        Check=False,
                    )
                    for sigma_upper, guess in zip(sigmas_upper, guesses)
                ]
            ).T
            ok = parameters[3].astype(bool)
            parameters = parameters[:, ok]

            mu_skew = parameters[0]
            sigma_skew = parameters[1]
            alpha_skew = parameters[2]

            filter_size = 11

            plt.figure()
            plt.plot(
                sigmas_upper,
                mu_skew,
                label=r"skew $\mu$",
                color="red",
                alpha=0.3,
                linewidth=3,
            )
            plt.plot(
                sigmas_upper,
                sigma_skew,
                label=r"skew $\sigma$",
                color="green",
                alpha=0.3,
                linewidth=3,
            )
            plt.plot(
                sigmas_upper,
                alpha_skew,
                label=r"skew $\alpha$",
                color="blue",
                alpha=0.3,
                linewidth=3,
            )
            plt.xscale("log")
            plt.legend(frameon=False)
            plt.xlabel(r"$\sigma_{upper}$")
            plt.title(
                f"mode = {mode}; "
                + r"$\mu_{original} = $"
                + f"{center}; "
                + r"$\sigma_{lower} = $"
                + f"{sigma_lower} | "
                + f"iteration #{i}"
            )

            plt.plot(sigmas_upper, mu_guess, color="red", linestyle="--")
            plt.plot(sigmas_upper, sigma_guess, color="green", linestyle="--")
            plt.plot(sigmas_upper, alpha_guess, color="blue", linestyle="--")
            plt.ylim(-30, 30)

            # update guess
            mu_guess = medfilt(mu_skew, filter_size)
            sigma_guess = medfilt(sigma_skew, filter_size)
            alpha_guess = medfilt(alpha_skew, filter_size)
            guesses = zip(mu_guess, sigma_guess, alpha_guess)

        t = Table(
            dict(
                sigma_upper=sigmas_upper,
                mu_skew=mu_guess,
                sigma_skew=sigma_guess,
                alpha_skew=alpha_guess,
            )
        )
        t.meta["mu"] = 0
        t.meta["sigma_lower"] = 1
        filename = f"skewnormal_parameters_for_interpolating_{mylabels[mode]}.ecsv"
        t.write(
            filename,
            overwrite=True,
        )
        print(f"saved parameter interpolating table to {filename}")
        tables[mylabels[mode]] = t
    return tables


interpolation_filename = os.path.join(
    code_directory, "populations/data/skewnormal_parameters_for_interpolating_mode.ecsv"
)
skew_table = ascii.read(interpolation_filename)


def make_skew_samples_from_lowerupper(
    mu=0, sigma_lower=1, sigma_upper=1, N_samples=1000
):
    """
    Generate a skew-normal distribution from a mode and asymmetric error bars.

    Parameters
    ----------
    mu : float, Quantity, and/or array-like
        The mode of the distribution(s) = the peak of the distribution.
        For nearly symmetric distributions the mode and mean and median are
        all very close, but for more skewed distributions this is making a
        (slightly arbitrary) strong choice about distribution shapes.

    sigma_lower : float, Quantity, and/or array-like
        The magnitude of the lower uncertainty = the distance from the mode
        mu down to the 15.87% confidence interval, which would be exactly
        1 sigma below the mode/median/mean of symmetric Gaussian.

    sigma_lower : float, Quantity, and/or array-like
        The magnitude of the lower uncertainty = the distance from the mode
        mu down to the 15.87% confidence interval, which would be exactly
        1 sigma below the mode/median/mean of symmetric Gaussian.

    N_samples : int
        The number of samples to generate.

    Returns
    -------


    """

    assert np.max([len(np.shape(x)) for x in [mu, sigma_lower, sigma_upper]]) <= 1

    N_data = np.max([np.size(mu), np.size(sigma_lower), np.size(sigma_upper)])

    # find coefficients for mu=0, sigma_lower=1, sigma_upper/sigma_lower
    mu_skew = np.interp(
        sigma_upper / sigma_lower, skew_table["sigma_upper"], skew_table["mu_skew"]
    )
    sigma_skew = np.interp(
        sigma_upper / sigma_lower, skew_table["sigma_upper"], skew_table["sigma_skew"]
    )
    alpha_skew = np.interp(
        sigma_upper / sigma_lower, skew_table["sigma_upper"], skew_table["alpha_skew"]
    )

    # populate empty array with nans, to skip asking skewnormal to handle them
    ok = np.isfinite(mu_skew * sigma_skew * alpha_skew)
    samples = np.ones((N_samples, N_data)) * mu[np.newaxis, :] * np.nan

    # create zero-centered, unnormalized skewnormal distribution + samples
    p = skewnorm(loc=mu_skew[ok], scale=sigma_skew[ok], a=alpha_skew[ok])
    samples[:, ok] = p.rvs(size=(N_samples, sum(ok))) * sigma_lower[ok] + mu[ok]

    return samples.T
