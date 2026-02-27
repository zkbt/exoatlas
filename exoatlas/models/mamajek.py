from .relation import *
import pandas as pd

__all__ = ["Mamajek"]


class Mamajek(Relation):
    """Eric Mamajek's relations for stars, linking spectral type, effective temperature, bolometric corrections, colors, abosolute magnitudes, luminosities, and mass,"""

    def __init__(self, filename="data/mamajek_dwarfproperties.txt"):
        Relation.__init__(self, filename)
        # self.table = self.table[np.argsort(self.table['Teff'])]

        # add some details
        self.name = "Mamajek"
        self.bibcode = "2013ApJS..208...9P"

        self.table = self.table.rename(
            columns={"R_Rsun": "radius", "Msun": "mass"}
        ).drop(columns="SpT2")
        self.table["logM"] = np.log10(self.table["mass"])
        self.table["logR"] = np.log10(self.table["radius"])

        # describe the columns in the relation
        self.descriptions = {
            "SpT": "(string) Spectral type of the star.",
            "Teff": "Stellar Effective Temperature, in K.",
            "logT": "log10(Teff), where Teff is in K.",
            "BCv": "Bolometric Correction from V band, in magnitudes.",
            "Mv": "Absolute V magnitude for a dwarf star.",
            "logL": "log10(Luminosity), in solar luminosities.",
            "mass": "Stellar mass, for dwarf stars, tentatively, in solar masses.",
            "M_G": "Absolute Gaia G magnitude.",
            "Mv": "Absolute V magnitude.",
            "M_J": "Absolute J magnitude.",
            "M_Ks": "Absolute Ks magnitude",
            "Mbol": "Absolute bolometric magnitude.",
            "radius": "Stellar Radius, in solar radii",
            "logM": "log10(stellar mass)",
            "logR": "log10(stellar radius)",
        }

        # create some additional color combinations that might be nice
        self.table["V-J"] = self.table["V-Ks"] - self.table["H-Ks"] - self.table["J-H"]

        # all the colors are the same text, so set them automatically
        for k in self.possible:
            if "-" in k:
                self.descriptions[k] = "{} color".format(k)

    def load(self, **kwargs):
        """
        Load the data table, and store it internally to this object.
        """

        # figure out the path to the data file, relative to this package
        # path = os.path.dirname(__file__) + '/'+ self.filename
        path = files(__name__) / self.filename

        # print("loading data from {0}".format(path))

        # load as an astropy table
        self.table = pd.read_csv(path, sep="\s+", comment="#")

        # give an update
        # print("   ...success!")


def test():
    stars = Mamajek()
    stars.plot()
    return stars
