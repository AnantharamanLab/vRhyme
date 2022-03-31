#! /usr/bin/env python3
# vRhyme
# Author: Kristopher Kieft
# University of Wisconsin-Madison


def read_nt_stuff_circ(in_fasta, min_len, folder, base):
    '''
    Read input fasta file and generate keep/mapper dicts
    '''
    keep = {}
    mapper = {}
    k = 1
    with open(f'{folder}/{base}.circular.tsv', 'w') as out:
        out.write('scaffold\ttype\tmismatches\tlength\trepeat\n')
        for header, seq, length, Ns in fasta_parse_circ(in_fasta):
            if length-Ns >= min_len:
                circ = Circular(seq, length)
                if circ.trType:
                    out.write(f'{header}\t{circ.trType}\t{circ.mismatch}\t{circ.trLen}\t{circ.trSeq}\n')
                else:
                    keep.update({header:(length,k)})
                    mapper.update({k:header})
                    k += 1

    return keep, mapper

def read_nt_stuff_interest_circ(in_fasta, min_len, interest_list, folder, base):
    '''
    Read input fasta file and generate keep/mapper dicts;
    filter for scaffolds of interest
    '''
    keep = {}
    mapper = {}
    k = 1
    with open(f'{folder}/{base}.circular.tsv', 'w') as out:
        out.write('scaffold\ttype\tmismatches\tlength\trepeat\n')
        for header, seq, length, Ns in fasta_parse_circ(in_fasta):
            if length-Ns >= min_len and header in interest_list:
                circ = Circular(seq, length)
                if circ.trType:
                    out.write(f'{header}\t{circ.trType}\t{circ.mismatch}\t{circ.trLen}\t{circ.trSeq}\n')
                else:
                    keep.update({header:(length,k)})
                    mapper.update({k:header})
                    k += 1

    return keep, mapper


def fasta_parse_circ(infile):
    '''
    Parse fasta file
    Check if sequence is not circular
    Return header and length of sequence - Ns
    '''
    with open(infile, 'r') as fasta:
        for line in fasta:
            if line[0] == '>':
                try:
                    seq = ''.join(seq)
                    Ns = seq.count("N")
                    yield header, seq, len(seq), Ns
                except NameError:
                    pass # first line
                header = line[1:].strip("\n")
                seq = []
            else:
                seq.append(line.strip("\n"))


        # last one
        seq = ''.join(seq)
        Ns = seq.count("N")
        yield header, seq, len(seq), Ns

#####

def read_nt_stuff(in_fasta, min_len):
    '''
    Read input fasta file and generate keep/mapper dicts
    '''
    keep = {}
    mapper = {}
    k = 1
    for header,length in fasta_parse(in_fasta):
        if length >= min_len:
            keep.update({header:(length,k)})
            mapper.update({k:header})
            k += 1

    return keep, mapper


def read_nt_stuff_interest(in_fasta, min_len, interest_list):
    '''
    Read input fasta file and generate keep/mapper dicts;
    filter for scaffolds of interest
    '''
    keep = {}
    mapper = {}
    k = 1
    for header,length in fasta_parse(in_fasta):
        if length >= min_len and header in interest_list:
            keep.update({header:(length,k)})
            mapper.update({k:header})
            k += 1

    return keep, mapper


def fasta_parse(infile):
    '''
    Parse fasta file
    Return header and length of sequence - Ns
    '''
    with open(infile, 'r') as fasta:
        for line in fasta:
            if line[0] == '>':
                try:
                    yield header, length
                except NameError:
                    pass # first line
                header = line[1:].strip("\n")
                length = 0
            else:
                seq = line.strip("\n")
                Ns = seq.count("N")
                length = length + len(seq) - Ns

        # last one
        yield header, length


####

class Circular:
    '''
    Flag sequences that are circular (i.e., likely represent complete viral genomes).

    DTR: direct terminal repeats
    ITR: inverted terminal repeats
    '''
    def __init__(self, seq, l):
        self.seq = seq
        self.l = l
        self.main()
    
    def main(self):
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
            self.end_sub = self.end[self.idx+20:] # skip the seed
            self.getRemaining()
            self.trType = 'DTR'
        except ValueError:
            # error, no index
            self.trType = None
    
    def itr(self):
        try:
            self.idx = self.c_seq.index(self.seed)
            self.end_sub = self.c_seq[self.idx+20:] # skip the seed
            self.getRemaining()
            self.trType = 'ITR'
        except ValueError:
            # error, no index
            self.trType = None
    
    def getRemaining(self):
        self.start_sub = self.start[20:]
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
                if self.mismatch >= 2:
                    break
        
        if self.mismatch > 0:
            while True:
                if any(mHolder[-4:]): # no mismatches at end, check mis_2 at end
                    if mHolder[-1]:
                        self.mismatch -= 1 # popped holder is a mismatch
                    iHolder.pop()
                    mHolder.pop()
                else:
                    break
        
        self.trSeq = self.seed + ''.join(iHolder)
        self.trLen = len(self.trSeq)

    def info(self):
        if self.l > 10000:
            self.maxRep = 5000
        else:
            self.maxRep = int(self.l/2)
        self.start = self.seq[:self.maxRep] # first _maxRep_ nucleotides of sequence
        self.seed = self.start[:20] # first _length_ nucleotides of start

    def revcomp(self):
        complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'a': 't', 'c': 'g', 'g': 'c', 't': 'a', '\n': '\n'}
        self.c_seq = ''.join([complement.get(base, base) for base in self.seq[-self.maxRep:]][::-1]) # last _maxRep_ nucleotides of reverse complement sequence

#
#
#
