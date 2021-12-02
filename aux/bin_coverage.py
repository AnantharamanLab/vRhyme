#! /usr/bin/env python3
# vRhyme
# Author: Kristopher Kieft
# University of Wisconsin-Madison


import argparse
import os
import sys
import numpy as np


class BinCoverage:

    def __init__(self, folder, covtable, memfile, outfile, min_cov, present, round_int, skip_std):
        self.outfile = outfile
        self.folder = folder
        self.min = min_cov
        self.present = present
        self.round = round_int
        self.skip = skip_std

        if self.folder:
            self.get_files()
        else:
            self.covtable = covtable
            self.memfile = memfile

        self.get_cov()
        self.get_mem()
        self.get_avg()


    def get_files(self):
        if self.folder[-1] != '/':
            self.folder == '/'
        
        files = os.listdir(self.folder)
        mem = [f for f in files if f.endswith('membership.tsv')]
        self.memfile = f'{self.folder}{mem[0]}'
        self.covtable = f'{self.folder}vRhyme_coverage_files/vRhyme_coverage_values.tsv'
        
        if not os.path.exists(self.memfile):
            sys.stderr.write("\nError: could not identify membership file from -v folder. Exiting.\n\n")
            exit()
        if not os.path.exists(self.covtable):
            sys.stderr.write("\nError: could not identify coverage file from -v folder. Exiting.\n\n")
            exit()


    def get_cov(self):
        self.coverages = {} # {scaffold:[averages]}
        with open(self.covtable, 'r') as cov:
            self.header = cov.readline().strip('\n').split('\t')[1:]
            for line in cov:
                line = line.strip('\n').split('\t')
                if not line: continue
                self.coverages[line[0]] = [float(i) for i in line[1::2]]

    def get_mem(self):
        self.bins = {} # {bin:[members]}
        with open(self.memfile, 'r') as mem:
            next(mem)
            for line in mem:
                line = line.strip('\n').split('\t')
                if not line: continue
                self.bins.setdefault(line[1], []).append(line[0])
    
    def get_avg(self):
        with open(self.outfile, 'w') as out:
            if self.skip:
                self.header = '\t'.join(self.header[::2])
            else:
                self.header = '\t'.join(self.header)
            if self.present:
                out.write(f'bin\tpresent\t{self.header}\n')
            else:
                out.write(f'bin\t{self.header}\n')

            for key,vals in self.bins.items():
                arr = np.array([self.coverages[v] for v in vals])
                if self.round:
                    a = [round(i,0) for i in np.average(arr, axis=0)]
                    if not self.skip:
                        s = [round(i,0) for i in np.std(arr, axis=0)]
                else:
                    a = [i for i in np.average(arr, axis=0)]
                    if not self.skip:
                        s = [i for i in np.std(arr, axis=0)]
                if self.skip:
                    z = '\t'.join([str(i) for i in a])
                else:
                    z = '\t'.join([str(i) for pair in zip(a, s) for i in pair])
                if self.present:
                    pres = sum([i>=self.min for i in a])
                    out.write(f'{key}\t{pres}\t{z}\n')
                else:
                    out.write(f'{key}\t{z}\n')


if __name__ == '__main__':
 
    descript = """
    Calculate the average coverage of a bin per sample.
    Simply takes the average of all scaffolds' averages.
    Use -v alone or specify -c and -m. 
"""
    vRhyme = argparse.ArgumentParser(description=descript, formatter_class=argparse.RawTextHelpFormatter, usage=argparse.SUPPRESS)
    vRhyme.add_argument('-v', metavar='', type=str, nargs=1, default=[''], help="input vRhyme output folder (skip -c/-m)")
    vRhyme.add_argument('-c', metavar='', type=str, nargs=1, default=[''], help="input vRhyme format coverage file (skip -v, requires -m)")
    vRhyme.add_argument('-m', metavar='', type=str, nargs=1, default=[''], help="input vRhyme format membership file (skip -v, requires -c)")
    vRhyme.add_argument('-o', metavar='', type=str, nargs=1, required=True, help='output bin coverage table')
    vRhyme.add_argument('--min', metavar='', type=str, nargs=1, default=['1.0'], help='mimimum coverage to consider as present (only with --present) [1.0]')
    vRhyme.add_argument('--present', action='store_true', help='add column in output for the number of samples a bin is present')
    vRhyme.add_argument('--round', action='store_true', help='round all values to whole integers')
    vRhyme.add_argument('--skip_std', action='store_true', help='do not calculate standard deviations')
    #
    args = vRhyme.parse_args()
    #
    #
    if os.path.exists(args.o[0]):
        sys.stderr.write("\nError: output table (-o) already exists. Exiting.\n\n")
        exit()
    
    m = float(args.min[0])
    BinCoverage(args.v[0], args.c[0], args.m[0], args.o[0], m, args.present, args.round, args.skip_std)

#
#
#
