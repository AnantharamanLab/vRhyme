#! /usr/bin/env python3
# vRhyme
# Author: Kristopher Kieft
# University of Wisconsin-Madison


def read_nt_stuff(in_fasta, min_len):
    '''
    Read input fasta file and generate keep/mapper dicts
    '''
    keep = {}
    mapper = {}
    spaces = False
    k = 1
    with open(in_fasta, 'r') as fasta:
        seq = ''
        for line in fasta:
            if line.startswith(">"):
                Ns = seq.count("N")
                length = len(seq)
                if length-Ns >= min_len:
                    keep.update({header:(length,k)})
                    mapper.update({k:header})
                    k += 1
                seq = ''
                header = line[1:].strip("\n") # no ">", no "\n"
                if " " in header:
                    spaces = True
            else:
                seq += line.strip("\n")

        # last one
        Ns = seq.count("N")
        length = len(seq)
        if length-Ns >= min_len:
            keep.update({header:(length,k)})
            mapper.update({k:header})

    return keep, spaces, mapper

def read_nt_stuff_interest(in_fasta, min_len, interest_list):
    '''
    Read input fasta file and generate keep/mapper dicts;
    filter for scaffolds of interest
    '''
    keep = {}
    mapper = {}
    spaces = False
    k = 1
    with open(in_fasta, 'r') as fasta:
        seq = ''
        for line in fasta:
            if line.startswith(">"):
                length = len(seq)
                if seq != '' and length >= min_len and header in interest_list:
                    keep.update({header:(length,k)})
                    mapper.update({k:header})
                    k += 1
                seq = ''
                header = line[1:].strip("\n") # no ">", no "\n"
                if " " in header:
                    spaces = True
            else:
                seq += line.strip("\n")

        # last one
        length = len(seq)
        if length >= min_len and header in interest_list:
            keep.update({header:(length,k)})
            mapper.update({k:header})

    return keep, spaces, mapper


#
#
#
