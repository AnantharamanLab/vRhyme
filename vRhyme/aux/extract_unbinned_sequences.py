#! /usr/bin/env python3
# vRhyme
# Author: Kristopher Kieft
# University of Wisconsin-Madison


import argparse
import os
import subprocess
import sys


class Unbinned:
    def __init__(self, infile, fasta, output):
        self.infile = infile
        self.fasta = fasta
        self.output = output
        self.bins = []

        self.get_bins()
        self.get_unbinned()

    def get_bins(self):
        with open(self.infile, 'r') as binsfile:
            next(binsfile)
            for line in binsfile:
                try:
                    line = line.rstrip("\n").split("\t")
                    self.bins.append(line[0])
                except Exception:
                    if line == '':
                        pass
                    else:
                        sys.stderr.write("\nError: The membership file could not be parsed. Is it in vRhyme format? Exiting.\n\n")
                        exit()
        #                
        self.bins = set(self.bins)
        if self.bins == set():
            sys.stderr.write("\nNo binned sequences were identified. Verify that the input file (-i) format is correct. Exiting.\n\n")
            exit()
        

    def get_unbinned(self):
        with open(self.fasta, 'r') as infasta, open(self.output, 'w') as outfile:
            seq = ''
            for line in infasta:
                if line.startswith(">"):
                    if seq != '':
                        if header not in self.bins:
                            outfile.write(f'>{header}\n{seq}\n')
                        seq = ''
                    header = line[1:].rstrip("\n")
                else:
                    seq += line.strip("\n")

            # last one
            if seq != '':
                if header not in self.bins:
                    outfile.write(f'>{header}\n{seq}\n')



if __name__ == '__main__':

    descript = """
    Extract unbinned sequences from input sequence file.
    Write all unbinned sequences to separate file.
    """
    vRhyme = argparse.ArgumentParser(description=descript, formatter_class=argparse.RawTextHelpFormatter)
    vRhyme.add_argument('-i', metavar='', type=str, nargs=1, required=True, help="input bins membership file")
    vRhyme.add_argument('-o', metavar='', type=str, nargs=1, required=True, help='name of output file to deposit unbinned sequences (fasta format)')
    vRhyme.add_argument('-f', metavar='', type=str, nargs=1, required=True, help='input scaffold sequences file')
    #
    args = vRhyme.parse_args()
    infile = str(args.i[0])
    output = str(args.o[0])
    fasta = str(args.f[0])
    #
    #
    if os.path.exists(output):
        sys.stderr.write("\nError: The output file already exists. Exiting.\n\n")
        exit()

    Unbinned(infile, fasta, output)

#
#
#
