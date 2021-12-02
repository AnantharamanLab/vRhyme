#! /usr/bin/env python3
# vRhyme
# Author: Kristopher Kieft
# University of Wisconsin-Madison

import setuptools

def get_descript():
    with open("README.md", "r") as dfile:
        d = dfile.read()
    return d

def get_version():
    with open("VERSION", 'r') as vfile:
        v = vfile.readline().strip()
    return v

def do_setup(v, d):
    setuptools.setup(
        name="vRhyme",
        version=v,
        author="Kristopher Kieft",
        author_email="kieft@wisc.edu",
        description="vRhyme: binning virus genomes from metagenomes",
        long_description_content_type="text/markdown",
        url="https://github.com/AnantharamanLab/vRhyme",
        classifiers=[
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License (GPLv3)"
        ],
        install_requires=['scikit-learn>=0.23.0', 'numpy>=1.17.0', 'numba>=0.50.0', 'pandas>=1.0.0', 'pysam>=0.15'],
        keywords=["bioinformatics", "metagenomics", "virus", "phage", "binning", "machine learning"],
        python_requires=">=3.6",
        scripts=['vRhyme'],
        packages=setuptools.find_packages(),
        license='GPLv3',
        long_description=d
    )


if __name__ == "__main__":
    do_setup(get_version(), get_descript())