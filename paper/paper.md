---
title: 'exoatlas: friendly Python code for exoplanet populations'
tags:
  - Python
  - astronomy
  - exoplanets
  - Solar System
authors:
  - name: Zach K. Berta-Thompson
    orcid: 0000-0002-3321-4924
    affiliation: 1
  - name: ???
    affiliation: ???
affiliations:
 - name: Department of Astrophysical and Planetary Sciences, University of Colorado Boulder, USA
   index: 1
   ror: 02ttsq026

date: 16 June 2025?
bibliography: paper.bib

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
aas-doi: ??? <- update this with the DOI from AAS once you know it.
aas-journal: The Astrophysical Journal Letters <- The name of the AAS journal.
---

# Summary

Planets are complicated. Understanding how they work requires context, comparing individual objects to relevant broader populations. Exoplanets are easier to imagine when seen next to their closest Solar System archetypes, and the planets in the Solar System look less lonely when viewed as one of thousands in the growing community of known exoplanets in the Milky Way. The `exoatlas` toolkit provides a friendly Python interface for retrieving and working with populations of planets, aiming to simplify the process of placing worlds in context. 

# Statement of need

We crafted `exoatlas` to meet a need among both researchers and educators for an intuitive Python tool to access planet populations. Particularly in working with students and junior scientists, for whom easy avenues for exploration and play would have particular benefit, we encountered significant barriers performing the following tasks:

- 🌎 retrieving basic properties for exoplanets + Solar System objects 
- 🧮 calculating derived planet quantities with propagated uncertainties
- 🗺️ comparing individual exoplanets to relevant comparison samples
- 🔭 planning future telescope observations of known exoplanet systems
- 🧑‍🏫 making beautiful planet plots for up-to-date data-rich lectures

Online planetary data archives have merged incredible curatorial efforts with powerful tools for data access and visualization, including the NASA Exoplanet Archive [@exo-archive], exoplanet.eu [@???], exo.MAST [@???], and Open Exoplanet Catalog [@open-exo] for exoplanets and the JPL Solar System Dynamics [@???], IAU Minor Planet Center [@???], the NSSDCA Planetary Facts sheets [@???] for Solar System objects. `exoatlas` does not intend to replace any of these important archival efforts (it pulls nearly all its data from the NASA Exoplanet Archive and JPL Solar System Dynamics datasets); rather, it aims to provide an approachable interface for exploratory analysis and illuminating visualizations to help the community make better use of these resources. 

# Mapping populations with `exoatlas` 

The user interface for `exoatlas` centers around the `Population` class, with each `Population` object containing a standardized table of planet properties and methods for interacting with that table. `exoatlas` makes extensive use of `astropy` [@astropy] to be as familiar as possible for modern astronomers, and is thoroughly documented with an astronomer audience in mind at [zkbt.github.io/exoatlas/](https://zkbt.github.io/exoatlas/). 

🌎 To retrieve planet data, `exoatlas` provides classes that inherit from `Population` and automatically access archive data for exoplanets as well as Solar System major planets, minor planets, and moons. Exoplanet data comes from the NASA Exoplanet Archive [@exo-archive] API, and Solar System data comes from the JPL Solar System Dynamics [@???] API or small reformatted data tables included in the package repository itself. Whatever the original data source, all `Population` objects act similarly and have uniform nomenclature for accessing data columns. Planet quantities all have physical units attached with `astropy.units` to facilitate unit conversions and minimize conceptual errors. New `Population` objects can also be created from `astropy.table` tables, enabling custom datasets to be included.

🧮 To calculate derived quantities for planets, a set of default methods are included within the core `Population` definition, or users may attach their own new calculation methods. If data for some quantity is missing, calculations can be used to  swap in alternate estimates; for example, a planet's semimajor axis $a$ will attempt first to pull from the original data table, and then second to calculate $a$ from the planet's period $P$ and the star's mass $M_\star$ assuming Newton's Version of Kepler's Third Law $P^2 = 4\pi^2 a^3/GM_\star$, and then third to calculate from a transit-derived ratio $a/R_*$. Uncertainties on all derived quantities can be numerically propagated using the `astropy.uncertainty` framework, where distributions of samples are generated for each original table quantity, carried through calculations, and then used to estimate confidence intervals.

🗺️ To extract sets of planets meeting particular criteria, `Population` objects can be indexed, sliced, and masked to generate new smaller `Population` objects. Coupling familiar array or table operations into the creation of subpopulations enables sophisticated filtering of datasets based on any combination of original archival table quantities, derived quantities, and/or quantity uncertainties. 

🔭 To plan telescope observations of known exoplanet systems, `exoatlas` can be used estimate the signal-to-noise (S/N) ratio acheivable for exoplanet observables with any particular telescope and to determine the visibility of targets from any particular observatory at any particular time. As derived quantities, all S/N estimates can include propagated uncertainties, allowing filtering to minimize biasing target samples toward systems with larger uncertainties. Using `astroplan` [@???] and `astropy.coordinates` [@???] to predict local altitude and azimuths, `exoatlas` can predict upcoming opportunities to observe transits from ground-based telescopes (when the star is above the horizon, the Sun is down, and the planet is passing in front of the star) for all elements in a `Population`. Including basic observation planning tools enables a workflow where targets can be filtered both by the estimated detectability of a signal and by whether a telescope can actually point at the star.

🧑‍🏫 To make illustrations of planet data, `exoatlas` includes a visual language to build up population comparisons. The core elements of this visual language are the `Plottable` (a quantity that should be represented with certain scaling, limits, and labels), the `Map` (a panel expressing plottable quantities with position, size, or color), and the `Gallery` (a collection of maps with linked axes and/or datasets). Preset visualizations built from this language can be used for quick contextual reference for a planet's fundamental properties, as in Figure \autoref{fig:exoatlas}, which was generated with only six lines of code:
```python
from exoatlas import *
from exoatlas.visualizations import *

e = TransitingExoplanets()
s = SolarSystem()
h = e["HD209458b"]

PlanetGallery().build([e, s, h])
```
![Example `exoatlas` visualization placing the first discovered transiting exoplanet HD209458b in context with other transiting exoplanets and the eight major Solar System planets. Errorbars use a color intensity that scales inversely with quantity uncertainties, to avoid giving undue visual weight to the least precise data. \label{fig:exoatlas}](joss-exoatlas-example.png)

# Research and learning
`exoatlas` was designed to support both the researcher who wants to contextualize planet populations in papers, proposals, or talks,  and the educator who wants help connect their students more immediately to real data from within the Solar System and beyond. The documentation aims to speak to an audience of both novice and expert users, hoping to combat quiet imposter thoughts by explicitly encouraging users to submit Issues for bugs and enhancements.

# Acknowledgements

We acknowledge the long commitment federally-funded archives have made to preserving and sharing data with the scientific community, and the heroic efforts of the people who build, maintain, and continually improve those archives.  

# References