#! /usr/bin/env python3
# vRhyme
# Author: Kristopher Kieft
# University of Wisconsin-Madison


import argparse
import os
import sys


class CircularBins:
    def __init__(self, bins, ext, output):
        self.bins = bins
        self.ext = ext
        self.output = output
        self.main()
    
    def main(self):
        files = os.listdir(self.bins)
        files = [f for f in files if f.endswith(self.ext)]

        with open(f'{self.output}/{self.output}.tsv', 'w') as out, open(f'{self.output}/{self.output}.fasta', 'w') as f:
            out.write('bin\tscaffold\ttype\tmismatches\tlength\trepeat\n')
            for file in files:
                base = file.rsplit('.',1)[0]
                for result in SeqTools(f'{self.bins}{file}'):
                    if result.trType:
                        out.write(f'{base}\t{result.header}\t{result.trType}\t{result.mismatch}\t{result.trLen}\t{result.trSeq}\n')
                        f.write(f'>{result.header}\n{result.seq}\n')

class CircularFasta:
    def __init__(self, fasta, output):
        self.fasta = fasta
        self.output = output
        self.main()
    
    def main(self):
        with open(f'{self.output}/{self.output}.tsv', 'w') as out, open(f'{self.output}/{self.output}.fasta', 'w') as f:
            out.write('scaffold\ttype\tmismatches\tlength\trepeat\n')
            for result in SeqTools(self.fasta):
                if result.trType:
                    out.write(f'{result.header}\t{result.trType}\t{result.mismatch}\t{result.trLen}\t{result.trSeq}\n')
                    f.write(f'>{result.header}\n{result.seq}\n')


class SeqTools: # generator
    def __init__(self, fasta):
        self.fasta = fasta
        # entry point is __iter__
    
    def looper(self):
        self.trLen = 0
        self.mismatch = 0
        self.info()
        self.dtr()
        if not self.trType:
            self.revcomp()
            self.itr()

    def dtr(self):
        self.end = self.seq[-self.maxRep:] # last _maxRep_ nucleotides of sequence
        try:
            self.idx = self.end.index(self.seed)
            self.end_sub = self.end[self.idx+length:] # skip the seed
            self.getRemaining()
            self.trType = 'DTR'
        except ValueError:
            # error, no index
            self.trType = None
    
    def itr(self):
        try:
            self.idx = self.c_seq.index(self.seed)
            self.end_sub = self.c_seq[self.idx+length:] # skip the seed
            self.getRemaining()
            self.trType = 'ITR'
        except ValueError:
            # error, no index
            self.trType = None
    
    def getRemaining(self):
        self.start_sub = self.start[length:]
        # start_sub may be longer than end_sub -> trimmed by zip
        z = zip(self.start_sub, self.end_sub)
        iHolder = [] # extend nucleotides
        mHolder = [] # holds mismatches
        for i,j in z:
            if i == j:
                iHolder.append(i)
                mHolder.append(False)
            else:
                self.mismatch += 1
                iHolder.append(i)
                mHolder.append(True)
                if self.mismatch >= mis:
                    break
        
        if self.mismatch > 0:
            while True:
                if any(mHolder[-mis_2:]): # no mismatches at end, check mis_2 at end
                    if mHolder[-1]:
                        self.mismatch -= 1 # popped holder is a mismatch
                    iHolder.pop()
                    mHolder.pop()
                else:
                    break
        
        self.trSeq = self.seed + ''.join(iHolder)
        self.trLen = len(self.trSeq)

    def info(self):
        self.l = len(self.seq)
        if self.l > maxx_2:
            self.maxRep = maxx
        else:
            self.maxRep = int(self.l/2)
        self.start = self.seq[:self.maxRep] # first _maxRep_ nucleotides of sequence
        self.seed = self.start[:length] # first _length_ nucleotides of start

    def revcomp(self):
            complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'a': 't', 'c': 'g', 'g': 'c', 't': 'a', '\n': '\n'}
            self.c_seq = ''.join([complement.get(base, base) for base in self.seq[-self.maxRep:]][::-1]) # last _maxRep_ nucleotides of reverse complement sequence

    def __iter__(self):
        with open(self.fasta, 'r') as f:
            for line in f:
                if line[0] == '>':
                    try: 
                        self.seq = ''.join(self.seq)
                        self.looper()
                        yield self
                    except AttributeError:
                        pass # first line
                    self.seq = []
                    self.header = line[1:].strip("\n")
                else:
                    self.seq.append(line.strip("\n"))

            # last one
            try:
                self.seq = ''.join(self.seq)
                self.looper()
                yield self
            except AttributeError:
                sys.stderr.write("\nError: Fasta file appears to be empty. Exiting.\n\n")
                exit()


class Entry:
    def __init__(self, bins, ext, fasta, output):
        self.bins = bins
        self.ext = ext
        self.fasta = fasta
        self.output = output

        self.main()
    
    def main(self):
        self.checkout()
        if self.bins:
            if self.bins[-1] != '/':
                self.bins += '/'
            CircularBins(self.bins, self.ext, self.output)
        elif self.fasta:
            CircularFasta(self.fasta, self.output)

    def checkout(self):
        if not self.output and self.bins:
            self.output = f'vRhyme_flag_circular_bins'
        elif not self.output and self.fasta:
            self.output = f'vRhyme_flag_circular_sequences'
        if os.path.exists(self.output):
            sys.stderr.write("\nError: The output file already exists. Exiting.\n")
            sys.stderr.write(f"{self.output}\n\n")
            exit()
        os.system(f'mkdir {self.output}')
        if self.output[-1] == '/':
            self.output = self.output[:-1]


if __name__ == '__main__':

    descript = """
    Flag scaffolds that are circular (i.e., likely represent complete viral genomes).
    Input vMAG bins (-b) to filter mis-binned complete viral genomes,
    or input any sequences (-f) to identify circular scaffolds. 

    Repeats
    DTR: direct terminal repeats
    ITR: inverted terminal repeats
"""
    vRhyme = argparse.ArgumentParser(description=descript, formatter_class=argparse.RawTextHelpFormatter)
    vRhyme.add_argument('-b', metavar='', type=str, nargs=1, default=[''], help="input folder containing fasta bins files")
    vRhyme.add_argument('-e', metavar='', type=str, nargs=1, default=['fasta'], help="extension of bins files [fasta]")
    vRhyme.add_argument('-f', metavar='', type=str, nargs=1, default=[''], help='input fasta sequences file')
    vRhyme.add_argument('-o', metavar='', type=str, nargs=1, default=[''], help='name of output folder')
    vRhyme.add_argument('-l', metavar='', type=int, nargs=1, default=[20], help='minimum terminal repeat length [20]')
    vRhyme.add_argument('-m', metavar='', type=int, nargs=1, default=[5000], help='maximum length of termini windows to search for repeats [5000]')
    vRhyme.add_argument('-n', metavar='', type=int, nargs=1, default=[2], help='maximum mismatches in repeats (seed of length -l must have 0) [2]')
    #
    args = vRhyme.parse_args()
    global length, maxx, maxx_2, mis, mis_2
    bins = args.b[0]
    ext = args.e[0]
    fasta = args.f[0]
    output = args.o[0]
    length = args.l[0]
    maxx = args.m[0]
    maxx_2 = maxx*2
    mis = args.n[0]
    mis_2 = mis*2
    #
    #
    if length < 10:
        sys.stderr.write("\nError: Minimum repeat length (-l) must be an integer >= 10. Exiting.\n\n")
        exit()
    if (bins and fasta) or (not bins and not fasta):
        sys.stderr.write("\nError: Provide either bins (-b) or sequences (-f). Exiting.\n\n")
        exit()

    Entry(bins, ext, fasta, output)

#
#
#
