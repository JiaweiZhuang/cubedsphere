#!/usr/bin/env python

from setuptools import setup, find_packages

MAJOR = 0
MINOR = 1
MICRO = 1
VERSION = "{}.{}.{}".format(MAJOR, MINOR, MICRO)

DESCRIPTION = "Cubed-Sphere data processing with xarray"
DISTNAME = "cubedsphere"
AUTHOR = "Jiawei Zhuang"
AUTHOR_EMAIL = "jiaweizhuang@g.harvard.edu"
URL = "https://github.com/jiaweizhuang/cubedsphere"
LICENSE = "MIT"
DOWNLOAD_URL = URL + "/archive/v{}.tar.gz".format(VERSION)

# TODO: Add classifiers for development status, science area
CLASSIFIERS = [
    "License :: OSI Approved :: MIT License",
    "Operating System :: OS Independent",
    "Programming Language :: Python",
    "Programming Language :: Python :: 3.4",
    "Programming Language :: Python :: 3.5",
    "Programming Language :: Python :: 3.6",
    "Topic :: Scientific/Engineering",
]

setup(
    name=DISTNAME,
    version=VERSION,
    description=DESCRIPTION,
    url=URL,
    download_url=DOWNLOAD_URL,
    author=AUTHOR,
    author_email=AUTHOR_EMAIL,
    license=LICENSE,

    packages=find_packages(),
    package_data={},
    scripts=[],
    install_requires=[
        'xarray', 'dask', 'numpy', 'matplotlib', 'cartopy'
    ],

    classifiers=CLASSIFIERS
)
