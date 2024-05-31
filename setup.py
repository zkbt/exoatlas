"""
This setup.py file sets up our package to be installable on any computer,
so that folks can `import rainbowconnection` from within any directory.
Thanks to this file, you can...

...tell python to look for code in the current directory (which you
can continue to edit), by typing *one* of the following commands:

`pip install -e .`
or
`python setup.py develop`

...move a copy of this code to your site-packages directory, where python will
be able to find it (but you won't be able to keep editing it), by typing *one*
of the following commands:

`pip install .`
or
`python setup.py install`

...upload the entire package to the Python Package Index, so that other folks
will be able to install your package via the simple `pip install rainbow-connection`, by
running the following command:

`python setup.py release`

The template for this setup.py came was pieced together with help from
barentsen, christinahedges, timothydmorton, and dfm. Check them out on github
for more neat tricks!

[`python-packaging`](https://python-packaging.readthedocs.io/) is useful too!
"""

# import our basic setup ingredients
from setuptools import setup, find_packages
import os, sys

# running `python setup.py release` from the command line will post to PyPI
if "release" in sys.argv[-1]:
    os.system("python setup.py sdist")
    # uncomment the next line to test out on test.pypi.com/project/tess-zap
    # os.system("twine upload --repository-url https://test.pypi.org/legacy/ dist/*")
    os.system("twine upload dist/*")
    os.system("rm -rf dist/exoplanet_atlas*")
    sys.exit()

# a little kludge to get the version number from __version__
exec(open("exoplanet_atlas/version.py").read())

# run the setup function
setup(
    # the name folks can use to search for this with pip
    name="exoplanet_atlas",
    # this package will only be installed if the current version doesn't exist
    version=__version__,
    # what's a short description of the package?
    description="Tools for working with transiting exoplanet populations.",
    # what's a more detailed description?
    long_description="For detailed usage, please read the documentation at https://zkbt.github.io/exoplanet_atlas/",
    # who's the main author?
    author="Zach Berta-Thompson",
    # what's the main author's email?
    author_email="zach.bertathompson@colorado.edu",
    # what's the URL for the repository?
    url="https://github.com/zkbt/exoplanet_atlas",
    # this figures out what subdirectories to include
    packages=find_packages(),
    # are there directories of data that should be accessible when installed?
    include_package_data=True,
    # where are those data directories?
    package_data={},
    # any scripts will be copied into your $PATH, to run from the command line
    scripts=[],
    # some descriptions about this package (for searchability?)
    classifiers=[
        "Intended Audience :: Education",
        "Intended Audience :: Science/Research",
        "Programming Language :: Python",
        "Topic :: Scientific/Engineering :: Astronomy",
    ],
    # what other packages are required. these must be pip-installable
    install_requires=[
        "numpy>=1.13",
        "matplotlib>=2.0",
        "astropy>=3.2.3",
        "astroquery>=0.3.9",
        "astroplan>=0.10",
        "rainbow-connection>=0.0.1",
        "PyYAML",
        "tqdm",
        "pytz",
    ],
    extras_require={
        "develop": [
            "pytest",
            "black",
            "black[jupyter]",
            "jupyter",
            "ipython",
            "mkdocs",
            "mkdocs-material",
            "mkdocstrings",
            "mkdocstrings-python",
            "pytkdocs[numpy-style]",
            "mkdocs-jupyter",
            "mkdocs-exclude",
            "twine",
            "pre-commit",
        ]
    },
    # what version of Python is required?
    python_requires=">=3.9",
    # (I think just leave this set to False)
    zip_safe=False,
    # under what license is this code released?
    license="MIT",
)
