from ..imports import *
from .downloaders import *

from astropy.table import join


class ExploreNEATables:
    def __init__(self, **downloadkw):

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

        # plot fractional difference
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
        plt.axhline(0.1, color="gray", linestyle="--")

        if True:

            plt.sca(ax[0, 1])
            plt.axis("off")
            s = f"{k}:\n\n"
            for t in ["ps", "default", "pscp"]:

                this = self.tables[t][k]
                N = np.sum(this.mask == False)
                s += f"{t} has data for {N/len(this):>5.1%} ({N}/{len(this)})\n"

            s += "\n"
            for t in ["ps", "default"]:
                this_joined = self.joined[f"pscp+{t}"]
                is_same = this_joined[f"{k}_{t}"] == this_joined[f"{k}_pscp"]
                N = np.sum(is_same)
                s += f"{t} == pscp  for {N/len(this_joined):>5.1%} ({N}/{len(this_joined)})\n"
            plt.text(
                0, 1, s, transform=ax[0, 1].transAxes, va="top", ha="left", fontsize=8
            )
