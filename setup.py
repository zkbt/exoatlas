from setuptools import setup, find_packages
import os,sys

def readme():
    with open('README.md') as f:
        return f.read()

# Hackishly inject a constant into builtins to enable importing of the
# package before the library is built.
import sys
if sys.version_info[0] < 3:
    import __builtin__ as builtins
else:
    import builtins
builtins.__EXOPOP_SETUP__ = True
import exopop
version = exopop.__version__

setup(name = "exopop",
    version = version,
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
    install_requires=['numpy>=1.9'],
    zip_safe=False
)
