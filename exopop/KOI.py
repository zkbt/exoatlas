# exoplanet population of all "confirmed" exoplanets from exoplanet archive
from imports import *
from Population import Population
from curation.KOI import correct


url ='http://exoplanetarchive.ipac.caltech.edu/cgi-bin/nstedAPI/nph-nstedAPI?table=cumulative&format=bar-delimited&select=*'
initial_filename = directories['data'] + 'exoplanetArchiveCumulativeCandidates.psv'

def downloadLatest():
    print 'downloading the latest list of confirmed exoplanets from the Exoplanet Archive'
    urllib.urlretrieve(url, initial_filename)
