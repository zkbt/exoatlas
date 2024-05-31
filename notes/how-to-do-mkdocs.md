# how to use mkdocs

christinahedges's very useful tutorial [here](https://christinahedges.github.io/astronomy_workflow/notebooks/3.0-building/mkdocs.html) goes through using `mkdocs` to make some automatic documentation.

To get started, install the developer installation. The `'.[develop]'` will get all the dependencies needed to generate the docs (= something like `mkdocs mkdocs-material mkdocstrings 'pytkdocs[numpy-style]' mkdocs-jupyter`).
```bash
git clone https://github.com/zkbt/exoplanet_atlas.git
cd exoplanet_atlas
pip install -e '.[develop]'
```

Then, we set up the docs by going into the base directory for this repository, and running
`mkdocs new .`
which made a `docs/` directory and a `mkdocs.yml`

Then, we copied Christina Hedges' template into the `docs/index.md` file and made appropriate changes for names.

Then, we edited the `docs/api.md` file to point to the objects and methods we want to explain. Docstrings for functions should follow the [`numpy` style conventions](https://numpydoc.readthedocs.io/en/latest/format.html). Options for rendering automatically are best described [here](https://mkdocstrings.github.io/python/usage/), although I don't fully understand everything going on in the YAML.

Then, we used the `mkdocs-jupyter` plugin to be able use jupyter notebooks as the source for writing docs, following the examples on their pages. We added some notebooks to the `docs/` diretory and pointed to them in the `mkdocs.yml` file.

Then, we ran `mkdocs serve`, and woah, a live version of the docs appeared at http://127.0.0.1:8000/. It was particularly cool (and way better than `sphinx` that I could make a change to any of the files and simply reload the page to see them update live into the docs). Hooray!

Then, we ran `mkdocs gh-deploy`, and double woah, it deployed a pretty version of the docs up at zkbt.github.io/chromatic! For the sake of not making the deployment `gh-pages` branch annoyingly large, add the `--no-history` option to erase the repository each time: `mkdocs gh-deploy --no-history`.
