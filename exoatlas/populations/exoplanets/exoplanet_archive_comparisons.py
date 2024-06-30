from ...imports import *
from .exoplanet_downloaders import *
from IPython.display import display

from astropy.table import join

kws = {
    "ps": dict(color="orchid"),
    "default": dict(color="black"),
    "pscp": dict(color="royalblue"),
}


class NASAExoplanetArchiveComparison:

    # define some variables to pay attention to in plots
    some_key_parameters = dict(
        stellar_basics=["st_mass", "st_rad", "st_dens", "st_lum", "st_teff"],
        stellar_fancy=["st_logg", "st_vsin", "st_met", "st_rotp", "st_age"],
        planet_basics=["pl_rade", "pl_bmasse", "pl_dens", "pl_orbper", "pl_orbsmax"],
        planet_fancy=[
            "pl_insol",
            "pl_eqt",
            "pl_imppar",
            "pl_ratror",
            "pl_ratdor",
            "pl_trueobliq",
            "pl_rvamp",
            "pl_orblper",
            "pl_orbtper",
            "pl_orbincl",
        ],
        planet_transit=[
            "pl_trandur",
            "pl_tranmid",
            "pl_trandep",
        ],
        stellar_positions=[
            "ra",
            "dec",
            "sy_dist",
            "sy_plx",
            "sy_pmra",
            "sy_pmdec",
            "st_radv",
        ],
    )

    def __init__(
        self, plot_directory="nasa-exoplanet-archive-table-comparisons", **downloadkw
    ):

        # empty dictionary of
        self.tables = {}

        # download both tables from NASA Exoplanet Archive
        self.tables["pscp"] = ExoplanetArchiveDownloader("pscomppars").get(**downloadkw)
        self.tables["ps"] = ExoplanetArchiveDownloader("ps").get(**downloadkw)

        # trim down to just default solutions
        is_default = self.tables["ps"]["default_flag"] == 1
        self.tables["default"] = self.tables["ps"][is_default]

        # create joined tables to be able to compare
        self.joined = {}
        for k in ["ps", "default"]:
            self.joined[f"pscp+{k}"] = join(
                self.tables["pscp"],
                self.tables[k],
                keys="pl_name",
                table_names=["pscp", k],
            )

        # create a directory to store plots
        self.plot_directory = plot_directory
        try:
            os.mkdir(plot_directory)
        except FileExistsError:
            pass

    def plot_parameter_comparison(self, k="pl_rade"):
        """
        Visualize how a parameter compares across the complete Planetary Systems,
        the default Planetary Systems, and Planetary Systems Composite Parameters.

        Parameters
        ----------
        k : str
            Which parameter to visualize.
        """
        plotkw = dict(marker=".", linewidth=0, alpha=0.5)

        fi, ax = plt.subplots(
            2,
            2,
            sharex="col",
            sharey="row",
            figsize=(8, 5),
            dpi=600,
            constrained_layout=True,
            gridspec_kw=dict(width_ratios=[1, 0.5]),
        )

        # plot one vs the other
        plt.sca(ax[0, 0])
        plt.title(k)
        for t in ["ps", "default"]:
            plt.loglog(
                self.joined[f"pscp+{t}"][f"{k}_pscp"],
                self.joined[f"pscp+{t}"][f"{k}_{t}"],
                label=t,
                **plotkw,
                **kws[t],
            )
        plt.ylabel(f'"{k}" from Planetary Systems')
        plt.legend(frameon=False, loc="upper left")

        # plot ratio, as points and as sideways histogram
        for t in ["ps", "default"]:
            ratio = (
                self.joined[f"pscp+{t}"][f"{k}_{t}"]
                / self.joined[f"pscp+{t}"][f"{k}_pscp"]
            )
            plt.sca(ax[1, 1])
            plt.hist(
                ratio, orientation="horizontal", bins=np.logspace(-3, 3, 300), **kws[t]
            )
            plt.xscale("log")

            plt.sca(ax[1, 0])
            plt.loglog(
                self.joined[f"pscp+{t}"][f"{k}_pscp"],
                ratio,
                label="all",
                **plotkw,
                **kws[t],
            )
        plt.xlabel(f'"{k}" from Planetary Systems Composite Parameters')
        plt.ylabel(f"ps/pscp")
        plt.ylim(0.5e-3, 2e3)
        plt.axhline(0.1, color="gray", linestyle="--")

        # display a text summary
        plt.sca(ax[0, 1])
        plt.axis("off")
        s = f"{k}:\n\n"
        for t in ["ps", "default", "pscp"]:
            this = self.tables[t][k]
            try:
                N = np.sum(this.mask == False)
            except AttributeError:
                N = np.sum(np.isfinite(this))
            s += f"{t} has data for {N/len(this):>5.1%} ({N}/{len(this)})\n"
        s += "\n"
        for t in ["ps", "default"]:
            this_joined = self.joined[f"pscp+{t}"]
            is_same = this_joined[f"{k}_{t}"] == this_joined[f"{k}_pscp"]
            N = np.sum(is_same)
            s += f"{t} == pscp  for {N/len(this_joined):>5.1%} ({N}/{len(this_joined)})\n"
        plt.text(0, 1, s, transform=ax[0, 1].transAxes, va="top", ha="left", fontsize=8)
        return s

    def do_completeness_check(self):
        """
        Check how complete each table is for some key parameters.
        """
        for group in self.some_key_parameters:
            # print(group)
            for k in self.some_key_parameters[group]:
                s = self.plot_parameter_comparison(k)
                # print(s)
                plt.savefig(
                    os.path.join(self.plot_directory, f"check-completeness-{k}.png")
                )
            # print

    def do_units_check(self, egregious_threshold=0.05):
        """
        Check how values compare when provided in Earth or Jupiter units.
        """
        # loop over radius and mass
        for p, q in zip(["rad", "bmass"], ["R", "M"]):

            plt.figure(figsize=(8, 5), constrained_layout=True)
            plt.title(f"'{p}' in different units")

            # loop over tables
            for k in ["ps", "default", "pscp"]:
                from_earth = self.tables[k][f"pl_{p}e"] * u.Unit(f"{q}earth")
                from_jupiter = self.tables[k][f"pl_{p}j"] * u.Unit(f"{q}jup")
                ratio = (from_earth / from_jupiter).decompose()
                plt.semilogx(
                    from_earth,
                    ratio,
                    marker=".",
                    linewidth=0,
                    alpha=0.25,
                    label=k,
                    **kws[k],
                )
                plt.ylabel(f"(pl_{p}e/pl_{p}j)" + r"$\times (M_{Jupiter}/M_{Earth})$")
                plt.ylim(0.75, 1.25)

                if k == "pscp":
                    is_egregious = np.abs(ratio - 1) > egregious_threshold
                    bad_planets = self.tables[k][is_egregious]["pl_name"]
                    print(
                        f"""
                    The parameter '{p}' is off by more than {egregious_threshold:%} for the {len(bad_planets)} planets 
                    {bad_planets}
                    simply because of weird rounding and/or Earth-Jupiter unit conversions.
                    """
                    )
                    print(list(self.tables[k][is_egregious]["pl_name"]))
            plt.legend(frameon=False)
            plt.savefig(os.path.join(self.plot_directory, f"check-units-{p}.png"))

            plt.show()

    def do_mass_radius_check(self):
        """
        Check that theoretical masses/radii get removed.
        """
        kw = dict(
            marker=".",
            linewidth=0,
            alpha=0.25,
        )
        plt.figure(figsize=(8, 5), constrained_layout=True)
        for k in ["ps", "default", "pscp"]:
            ok_mass_error = (self.tables[k]["pl_bmasseerr1"].mask == False) * (
                self.tables[k]["pl_bmasseerr2"].mask == False
            )
            ok_radius_error = (self.tables[k]["pl_radeerr1"].mask == False) * (
                self.tables[k]["pl_radeerr2"].mask == False
            )
            t = self.tables[k][ok_mass_error * ok_radius_error]
            plt.loglog(t["pl_bmasse"], t["pl_rade"], label=k, **kws[k], **kw)

        t = self.tables[k][ok_mass_error * ok_radius_error == False]
        plt.loglog(
            t["pl_bmasse"],
            t["pl_rade"],
            color="red",
            label="no mass|radius errors",
            **kw,
        )

        plt.axis("scaled")
        plt.legend(frameon=False)
        plt.xlabel("Planet Mass ($M_\oplus$)")
        plt.ylabel("Planet Radius ($R_\oplus$)")
        plt.savefig(os.path.join(self.plot_directory, f"check-mass-radius.png"))

    def do_luminosity_check(self):
        kw = dict(
            marker=".",
            linewidth=0,
            alpha=0.25,
        )
        fi, ax = plt.subplots(
            2, 1, figsize=(8, 5), constrained_layout=True, sharex=True
        )
        for k in ["ps", "default", "pscp"]:
            t = self.tables[k]  # [ok_mass_error * ok_radius_error]
            ok_lum = (t["st_lumerr1"].mask == False) * (t["st_lumerr2"].mask == False)
            ok_rad = (t["st_raderr1"].mask == False) * (t["st_raderr2"].mask == False)
            ok_teff = (t["st_tefferr1"].mask == False) * (
                t["st_tefferr2"].mask == False
            )
            ok = ok_lum * ok_rad * ok_lum
            T = t["st_teff"] * u.K
            L = 10 ** t["st_lum"] * u.Lsun
            R = t["st_rad"] * u.Rsun
            plt.sca(ax[0])
            x = (4 * np.pi * R**2 * con.sigma_sb * T**4).to(u.Lsun)
            y = L.to(u.Lsun)
            plt.loglog(x[ok], y[ok], label=k, **kws[k], **kw)
            plt.sca(ax[1])
            plt.semilogx(x[ok], (y / x)[ok], label=k, **kws[k], **kw)
            plt.ylim(0, 2)

            bad = ok == False
            plt.sca(ax[0])
            plt.loglog(
                x[bad],
                y[bad],
                label=f"{k} (no lum/teff/radius uncertainty)",
                color="red",
                **kws[k],
                **kw,
            )
            plt.plot([1e-6, 1e4], [1e-6, 1e4], linestyle="--", color="gray")
            plt.ylabel("L")
            plt.legend(frameon=False, bbox_to_anchor=(1, 1), loc="upper left")

            plt.sca(ax[1])
            plt.semilogx(x[bad], (y / x)[bad], color="red", **kws[k], **kw)
            plt.axhline(1, color="gray", linestyle="--")
            plt.ylabel("L / [$4\pi R^2 \sigma T_{eff}^4$]")

        plt.sca(ax[1])
        plt.xlabel("$4\pi R^2 \sigma T_{eff}^4$")
        plt.savefig(
            os.path.join(self.plot_directory, f"check-luminosity-radius-teff.png")
        )
