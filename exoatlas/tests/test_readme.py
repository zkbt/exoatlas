def test_readme():
    '''
    Test the code that appears on the GitHub README.
    '''

    # import some population definitions and plotting tool
    from exoatlas import TransitingExoplanets, SolarSystem
    from exoatlas.visualizations import PlanetGallery

    # create a dictionary of populations
    exo = TransitingExoplanets()
    solar = SolarSystem()

    # use a default visualization to summarize these populations
    PlanetGallery().build([solar, exo])