#! /usr/bin/env python3
# vRhyme
# Author: Kristopher Kieft
# University of Wisconsin-Madison


import argparse
import os
import subprocess
import sys


class Main:
    def __init__(self, infile, output, fasta, proteins, genes, ex_fasta, ex_proteins, ex_genes, prefix):
        self.infile = infile
        self.output = output
        self.fasta = fasta
        self.proteins = proteins
        self.genes = genes
        self.ex_fasta = ex_fasta
        self.ex_proteins = ex_proteins
        self.ex_genes = ex_genes
        self.prefix = prefix
        self.bins = {}

        try:
            temp = self.infile.rsplit('/',1)[1]
            self.base = temp.rsplit('.',1)[0]
        except Exception:
            self.base = self.infile.rsplit('.',1)[0]

        self.get_bins()
        if self.fasta != '':
            self.get_fasta()
        if self.proteins != '':
            self.get_proteins()
        if self.genes != '':
            self.get_genes()


    def get_bins(self):
        with open(self.infile, 'r') as binsfile:
            next(binsfile)
            for line in binsfile:
                line = line.rstrip("\n").split("\t")
                self.bins.update({line[0]:line[1]})

        if self.bins == {}:
            sys.stderr.write("\nNo bins were identified. Verify that the input file (-i) format is tab-separated. Exiting.\n\n")
            exit()


    def get_fasta(self):
        with open(self.fasta, 'r') as infasta:
            seq = ''
            for line in infasta:
                if line.startswith(">"):
                    if seq != '':
                        try:
                            b = self.bins[header]
                            p = self.prefix.replace("#",b)
                            with open(f'{self.output}{self.base}.bin_{b}.{self.ex_fasta}', 'a') as outfasta:
                                outfasta.write(f">{p}{header}\n{seq}\n")
                        except KeyError:
                            pass # unbinned
                        seq = ''
                    header = line[1:].rstrip("\n")
                else:
                    seq += line.strip("\n")

            try:
                b = self.bins[header]
                p = self.prefix.replace("#",b)
                with open(f'{self.output}{self.base}.bin_{b}.{self.ex_fasta}', 'a') as outfasta:
                    outfasta.write(f">{p}{header}\n{seq}\n")
            except KeyError:
                pass # unbinned

    def get_proteins(self):
        with open(self.proteins, 'r') as inproteins:
            seq = ''
            for line in inproteins:
                if line.startswith(">"):
                    if seq != '':
                        try:
                            b = self.bins[header_base.rsplit('_',1)[0]]
                            p = self.prefix.replace("#",b)
                            with open(f'{self.output}{self.base}.bin_{b}.{self.ex_proteins}', 'a') as outfasta:
                                outfasta.write(f">{p}{header}\n{seq}\n")
                        except KeyError:
                            pass # unbinned
                        seq = ''
                    header = line[1:].strip("\n")
                    try:
                        header_base = header.split(" # ",1)[0]
                    except KeyError:
                        header_base = header
                else:
                    seq += line.strip("\n")

            try:
                b = self.bins[header_base.rsplit('_',1)[0]]
                p = self.prefix.replace("#",b)
                with open(f'{self.output}{self.base}.bin_{b}.{self.ex_proteins}', 'a') as outfasta:
                    outfasta.write(f">{p}{header}\n{seq}\n")
            except KeyError:
                pass # unbinned


    def get_genes(self):
        with open(self.genes, 'r') as ingenes:
            seq = ''
            for line in ingenes:
                if line.startswith(">"):
                    if seq != '':
                        try:
                            b = self.bins[header_base.rsplit('_',1)[0]]
                            p = self.prefix.replace("#",b)
                            with open(f'{self.output}{self.base}.bin_{b}.{self.ex_genes}', 'a') as outfasta:
                                outfasta.write(f">{p}{header}\n{seq}\n")
                        except KeyError:
                            pass # unbinned
                        seq = ''
                    header = line[1:].strip("\n")
                    try:
                        header_base = header.split(" # ",1)[0]
                    except KeyError:
                        header_base = header
                else:
                    seq += line.strip("\n")

            try:
                b = self.bins[header_base.rsplit('_',1)[0]]
                p = self.prefix.replace("#",b)
                with open(f'{self.output}{self.base}.bin_{b}.{self.ex_genes}', 'a') as outfasta:
                    outfasta.write(f">{p}{header}\n{seq}\n")
            except KeyError:
                pass # unbinned


if __name__ == '__main__':

    descript = """
    Write sequences from vRhyme alternative bins to fasta files.
    Write scaffold sequences, proteins and/or genes.
    The output files will be named according to -i input.
    """
    vRhyme = argparse.ArgumentParser(description=descript, formatter_class=argparse.RawTextHelpFormatter)
    vRhyme.add_argument('-i', metavar='', type=str, nargs=1, required=True, help="input alternative bins membership file")
    vRhyme.add_argument('-o', metavar='', type=str, nargs=1, required=True, help='folder to deposit output files')
    vRhyme.add_argument('-f', metavar='', type=str, nargs=1, default = [''], help='input scaffold sequences file (optional)')
    vRhyme.add_argument('-p', metavar='', type=str, nargs=1, default = [''], help='input scaffold proteins file (optional)')
    vRhyme.add_argument('-g', metavar='', type=str, nargs=1, default = [''], help='input scaffold genes file (optional)')
    vRhyme.add_argument('-x', metavar='', type=str, nargs=1, default = ['bin_#__'], help="prefix to append to binned scaffold names, use the '#' symbol where the bin number will be noted [bin_#__]")
    #
    args = vRhyme.parse_args()
    infile = str(args.i[0])
    output = str(args.o[0])
    fasta = str(args.f[0])
    proteins = str(args.p[0])
    genes = str(args.g[0])
    prefix = str(args.x[0])

    if '#' not in prefix:
        sys.stderr.write("\nError: Use the '#' symbol in prefix -x to denote where to identify the bin number. Exiting.\n\n")
        exit()

    if fasta == '' and proteins == '' and genes == '':
        sys.stderr.write("\nError: At least one of -f/-p/-g, or any combination, must be supplied. Exiting.\n\n")
        exit()
    if fasta != '':
        ex_fasta = fasta.rsplit('.',1)[1]
    else:
        ex_fasta = 'ex_fasta'
    if proteins != '':
        ex_proteins = proteins.rsplit('.',1)[1]
    else:
        ex_proteins = 'ex_proteins'
    if genes != '':
        ex_genes = genes.rsplit('.',1)[1]
    else:
        ex_genes = 'ex_genes'
    

    if len(set([ex_fasta, ex_proteins, ex_genes])) != 3:
        sys.stderr.write("\nError: The file extensions of -f/-p/-g cannot be the same. Exiting.\n\n")
        exit()

    if output[-1] != '/':
        output += '/'
    #
    #
    if os.path.exists(output):
        sys.stderr.write("\nError: The output folder already exists. Exiting.\n\n")
        exit()

    subprocess.run('mkdir ' + output, shell=True)


    Main(infile, output, fasta, proteins, genes, ex_fasta, ex_proteins, ex_genes, prefix)
#
#
#
#
