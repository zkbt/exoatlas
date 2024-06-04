from ..imports import *
from .downloaders import *

from astropy.table import join


class ExploreNEATables:
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
        self.tables["pscp"] = composite_planetary_systems_downloader.get(**downloadkw)
        self.tables["ps"] = planetary_systems_downloader.get(**downloadkw)

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
        kws = {"ps": dict(color="orchid"), "default": dict(color="black")}
        fi, ax = plt.subplots(
            2,
            2,
            sharex="col",
            sharey="row",
            figsize=(8, 8),
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

    def summarize_some_key_parameters(self):
        for group in self.some_key_parameters:
            print(group)
            for k in self.some_key_parameters[group]:
                s = self.plot_parameter_comparison(k)
                print(s)
                plt.savefig(os.path.join(self.plot_directory, f"pscp-{k}.png"))
            print
