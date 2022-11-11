"""
Tools for working with the habitable zone definitions
in Kopparapu et al. (2013).
"""


coefficients = {}
coefficients["recent-venus"] = [1.7753, 1.4316e-4, 2.9875e-9, -7.5702e-12, -1.1635e-15]
coefficients["runaway-greenhouse"] = [
    1.0512,
    1.3242e-4,
    1.5418e-8,
    -7.9895e-12,
    -1.8328e-15,
]
coefficients["moist-greenhouse"] = [
    1.0140,
    8.1774e-5,
    1.7063e-9,
    -4.3241e-12,
    -6.6462e-16,
]
coefficients["maximum-greenhouse"] = [
    0.348,
    5.8924e-5,
    1.6558e-9,
    -3.0045e-12,
    -5.2983e-16,
]
coefficients["early-mars"] = [0.3179, 5.4513e-5, 1.5313e-9, -2.7786e-12, -4.8997e-16]


def make_hz(which="runaway-greenhouse"):
    """
    Create a function f(Teff) describing one of
    the habitable zone boundaries presented in
    Table 3 of Kopparapu et al. (2013).
    """

    # pull out the coefficients
    try:
        S_o, a, b, c, d = coefficients[which]
    except KeyError:
        error_message = f"""
        Alas! There seems to be no way to determine
        coefficients for your choice of `{which}`.
        The available options are:
        """
        for k in which:
            error_message += k + "\n"
        raise ValueError(error_message)

    # define the function
    def f(Teff):
        f"""
        The `{which}` HZ boundary from Kopparapu et al. (2013),
        as a function of stellar effective temperature.
        """

        T = Teff - 5780
        return S_o + a * T + b * T**2 + c * T**3 + d * T**4

    # return that function definition
    return f
