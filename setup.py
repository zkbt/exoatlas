from setuptools import setup, find_packages
import os,sys

def readme():
    with open('README.md') as f:
        return f.read()

# a little kludge to get the version number from __version__
exec(open('exopop/version.py').read())

setup(name = "exopop",
    version = __version__,
    description = "Tools for compiling and plotting populations of transiting exoplanets.",
    long_description = readme(),
    author = "Zachory K. Berta-Thompson",
    author_email = "zkbt@mit.edu",
    url = "https://github.com/zkbt/exopop",
    packages = find_packages(),
    package_data = {'':[]},
    scripts = [],
    classifiers=[
      'Intended Audience :: Science/Research',
      'Programming Language :: Python',
      'Topic :: Scientific/Engineering :: Astronomy'
      ],
    install_requires=['numpy>=1.13',
                      'matplotlib',
                      'astropy>=3.2',
                      'astroquery>=0.3.9',
                      'rainbow-connection>=0.0.1',
                      'tqdm'],
    python_requires='>3',
    zip_safe=False
)
