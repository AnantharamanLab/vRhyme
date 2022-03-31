#! /usr/bin/env python3
# vRhyme
# Author: Kristopher Kieft
# University of Wisconsin-Madison

import setuptools
import subprocess
import os

def get_descript():
    with open("README.md", "r") as dfile:
        d = dfile.read()
    return d

def get_version():
    with open("VERSION", 'r') as vfile:
        v = vfile.readline().strip()
    return v

def do_setup(v, d):
    try:
        subprocess.check_output("which gzip", shell=True)
    except Exception:
        print("\nError: gzip cannot be found. gzip is necessary for unzipping one of the files. Exiting.")
        exit(1)
    subprocess.run(f'gunzip models/vRhyme_machine_model_ET.sav.gz 2> /dev/null', shell=True)

    aux = [f'vRhyme/aux/{f}' for f in os.listdir('vRhyme/aux') if f.endswith('.py')]
    threaders = [f'vRhyme/scripts/{f}' for f in os.listdir('vRhyme/scripts') if f.endswith('.py')]
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
        install_requires=['scikit-learn>=0.23.0', 'numpy>=1.17.0', 'numba>=0.50.0', 'pandas>=1.0.0', 'pysam>=0.15', 'networkx>=2.0'],
        keywords=["bioinformatics", "metagenomics", "virus", "phage", "binning", "machine learning"],
        python_requires=">=3.6",
        scripts=['vRhyme/vRhyme'] + aux + threaders,
        package_data={'vRhyme': ['models/vRhyme_machine_model_ET.sav', 'models/vRhyme_machine_model_NN.sav', 'LICENSE', 'VERSION', 'README.md', 'README.pdf']},
        packages=setuptools.find_packages(),
        license='GPLv3',
        long_description=d
    )


if __name__ == "__main__":
    do_setup(get_version(), get_descript())