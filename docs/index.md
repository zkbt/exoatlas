# 🌌 exoatlas

Welcome! `exoatlas` is a Python package for interacting with basic properties of exoplanets and Solar System planets. It provides a friendly way to download useful data from online archives like the [NASA Exoplanet Archive](https://exoplanetarchive.ipac.caltech.edu) and [JPL Solar System Dynamics](https://ssd.jpl.nasa.gov), easy access to planet properties, and some pre-packaged visualizations for summarizing and comparing populations of planets.

## 📦 How to Install `exoatlas`
For most users, simply installing via 

```pip install exoatlas``` 

should do the trick. If you want more detail, either because running that command seems scary *or* because you want to do something fancier, then please check out the [installation instructions](installation.ipynb).

## 📖 How to Use this Documentation
Let's be honest, you're probably going to skip over most of the documentation to look for your own very specific use. Good for you, that sounds like the right way to go! 

The documentation is targeted at someone with basic knowledge in astronomy, planetary science, and Python. Each page is fairly self-contained, so feel free to start wherever you like. Here's a little map help guide your search more efficiently. 

| If you're here because you want to...                 | ...then try jumping to:  | 
| ----------------------------------------------------- | ---------------------- | 
| get a quick sense of `exoatlas`                       | [Quickstart](quickstart.ipynb) |
| learn how `exoatlas` defines populations              | [Populations](user/populations.ipynb) |
| retrieve planet properties with `exoatlas`            | [Creating](user/creating.ipynb) | 
| do calculations and define subsets in `exoatlas`      | [Filtering](user/filtering.ipynb) |
| paint pretty population plots with `exoatlas`         | [Visualizing](user/visualizing.ipynb) |
| use `exoatlas` to estimate signal-to-noise ratios     | [Observing](user/observing.ipynb) | 
| plan telescope observations of `exoatlas` planets     | [Planning](user/planning.ipynb) |
| understand how `exoatlas` propagates uncertainties    | [Uncertainties]( user/uncertainties.ipynb) |
| use `exoatlas`'s [BTWM26](https://ui.adsabs.harvard.edu/abs/2025arXiv250702136B/abstract) cosmic shoreline calculator   | [Shoreline](user/shoreline.ipynb) |
| dive into the nitty gritty of `exoatlas` definitions  | [API Reference](autoapi/) |


## 🗂️ How to Cite `exoatlas`
A manuscript is currently under review at JOSS (as of late spring 2026). Once that is complete, we will link to the SciEx entry where you can access bibliographic data for the code/paper in multiple formats. In the meantime, [this arxiv posting](https://scixplorer.org/abs/2025arXiv250702210B/abstract) is a good reference.

## 📞 How to Get Help

If you run into trouble with `exoatlas` or are curious about a new feature, and you can't obviously find a solution here in the documentation, please consider [submitting an Issue](developer/github.ipynb) on GitHub. If you're wondering about something, someone else might be too; be kind and help us all out by asking your question!