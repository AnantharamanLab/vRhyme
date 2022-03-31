#! /usr/bin/env python3
# vRhyme
# Author: Kristopher Kieft
# University of Wisconsin-Madison

import sys
import subprocess
import pkg_resources
import os

def helper(arg):
    if arg == '-h' or arg == '--help':
        h = '''
        This script is used to test the setup of vRhyme dependencies. 
        There are no input arguments, simply run the script. 
        All test results are printed to the standard out.
        "Success" statements mean the dependency should be good.

'''
        sys.stdout.write(h)
        exit()

def test_packages():
    sys.stdout.write('\n  Python Dependencies\n  -------------------\n')
    packages = {str(p).split(' ')[0]:str(p).split(' ')[1] for p in pkg_resources.working_set}
    packages_names = set(packages.keys())
    test = {'scikit-learn': 0.23, 'numpy': 1.17, 'numba': 0.50, 'pandas': 1.0, 'pysam': 0.15, 'networkx': 2.0}
    for p,v in test.items():
        if p in packages_names:
            ver = packages[p]
            if ver.count('.') >= 2:
                check = float('.'.join(ver.split('.')[:2]))
            else:
                check = float(ver)
            if check < v:
                sys.stdout.write(f'\033[1m> {p}: Update Required (v{ver})\033[0m\n')
            else:
                sys.stdout.write(f'  {p}: Success (v{ver})\n')
        else:
            sys.stdout.write(f'\033[1m> {p}: Not Found!\033[0m\n')
    sys.stdout.write('\n')


def test_software():
    sys.stdout.write('\n  Program Dependencies\n  --------------------\n')
    software = {
        'mmseqs': 'Not Found! Required',
        'samtools': 'Not Found! Usually Required',
        'prodigal': 'Not Found! Optional',
        'mash': 'Not Found! Optional',
        'nucmer': 'Not Found! Optional',
        'bowtie2': 'Not Found! Optional',
        'bwa': 'Not Found! Optional' 
    }
    for s,v in software.items():
        try:
            subprocess.check_output(f"which {s}", shell=True)
            sys.stdout.write(f'  {s}: Success\n')
        except Exception:
            sys.stdout.write(f'\033[1m> {s}: {v}\033[0m\n')
    sys.stdout.write('\n')  


def test_model():
    sys.stdout.write('\n  Machine Learning Models\n  -----------------------\n')

    base_scripts = os.path.dirname(os.path.abspath(__file__)).replace('/aux','')
    base_bin = os.path.dirname(os.path.abspath(__file__)).replace('/bin','')
    if os.path.exists(f'{base_scripts}/models/vRhyme_machine_model_ET.sav.gz'):
        try:
            subprocess.check_output("which gzip", shell=True)
        except Exception:
            sys.stdout.write(f'  unzip ET model: Failed (gzip not found)\n')
        subprocess.run(f'gunzip {base_scripts}/models/vRhyme_machine_model_ET.sav.gz', shell=True)
    for ml in ['NN', 'ET']:
        if os.path.exists(f'{base_scripts}/models/vRhyme_machine_model_{ml}.sav'):
            sys.stdout.write(f'  {ml} model: Success\n')
        elif os.path.exists(f'{base_bin}/lib/python{sys.version_info.major}.{sys.version_info.minor}/site-packages/vRhyme/models/vRhyme_machine_model_{ml}.sav'):
            sys.stdout.write(f'  {ml} model: Success\n')
        else:
            sys.stdout.write(f'  {ml} model: Not Found!\n')
    sys.stdout.write('\n')
    
if __name__ == '__main__':
    try:
        helper(sys.argv[1])
    except IndexError:
        pass
    test_packages()
    test_software()
    test_model()