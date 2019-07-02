#! /usr/bin/env python
from setuptools import setup, find_packages


setup(
    name="permamodel",
    version="0.1.2",
    author="Elchin Jafarov and Scott Stewart",
    author_email="james.stewart@colorado.edu",
    description="Permamodel",
    long_description=open("README.md").read(),
    packages=find_packages(),
    # install_requires=('numpy', 'gdal', 'pyproj'),
    install_requires=(
        "affine",
        "bmipy",
        "netCDF4",
        "numpy",
        "python-dateutil",
        "pyyaml",
        "scipy",
    ),
    package_data={"": ["examples/*.cfg", "examples/*.dat", "data/*"]},
)
