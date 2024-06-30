from exoatlas.ThumbtackPlot import ThumbtackPlot
from exoatlas.transiting_exoplanets import GoodMass, BadMass, Kepler, NonKepler
from exoatlas.KOI import UnconfirmedKepler
from exoatlas.TESS import TESS

pops = dict(
    kepler=Kepler(), candidates=UnconfirmedKepler(), nonkepler=NonKepler(), new=TESS()
)


t = ThumbtackPlot(pops=pops, lightyears=False)
pops["kepler"].label = "Kepler"
pops["candidates"].label = None
pops["candidates"].color = "royalblue"
pops["kepler"].color = "royalblue"
pops["nonkepler"].color = "black"

pops["nonkepler"].zorder = 0
pops["kepler"].zorder = 1
pops["new"].zorder = 2
pops["new"].label = "Predicted TESS"
pops["new"].standard["name"] = ""
pops["new"].color = "darkorange"


t.movie(keys=["candidates", "kepler", "nonkepler", "new"], maxdistance=1300.0)
