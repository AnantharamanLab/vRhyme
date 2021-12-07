#! /usr/bin/env python3
# vRhyme
# Author: Kristopher Kieft
# University of Wisconsin-Madison

import numpy as np
import collections
import itertools
from numba import jit
import pickle
import sys
try:
    from scripts import save_dist_dicts as info
except (ModuleNotFoundError, ImportError):
    import save_dist_dicts as info

@jit(nopython=True)
def denomer_generate(in_list, in_avg):
    '''
    Denominator for pairwise cosine distance calculations.
    Range is for all 64 codons.
    '''
    denomer = 0
    for n in range(0,64):
        denomer += (in_list[n]-in_avg)**2
    return denomer

@jit(nopython=True)
def array_generate(datalist, avg):
    '''
    Begin calculation of numerator for cosine pairwise distance calculations.
    Range is for all 64 codons.
    '''
    for n in range(0,64):
        datalist[n] -= avg
    return datalist

def codon_counter(seq, header, codon_list):
    '''
    Count occurances of codons within genes and
    calculate numerator/denominator values for cosine pairwise distance calculations.
    Store information within transferable info.data_dict_codon
    '''
    seq = seq.upper()

    length = len(seq)
    kmer_list = [seq[i:i+3] for i in range(0,length-2,3)] # forward, -1 is so only nmers are made, step of 3
    seq = None
    kmer_list = [x for x in kmer_list if not '*' in x]
    codon_dict = dict(collections.Counter(kmer_list)) # count kmers
    kmer_list = None

    total = sum(codon_dict.values())
    codon_vals = []
    for val in codon_list:
        try:
            codon_vals.append(codon_dict[val]/total)
        except Exception:
            codon_vals.append(0)

    codon_vals = np.array(codon_vals)
    c_avg = np.mean(codon_vals)
    denomer = denomer_generate(codon_vals, c_avg)
    codon_vals = array_generate(codon_vals, c_avg)
    info.data_dict_codon.update({header:(codon_vals,denomer)})

######################## lists
def nmer_maker():
    '''
    All possible 3mers (i.e., codons) for ATCG
    '''
    codon_list = [''.join(mer) for mer in itertools.product('ATCG', repeat=3)]
    return codon_list

#######################

def codon_stuff(in_genes, pickle_folder, file_n, spaces, keep):
    '''
    Wrapper for all functions, including reading genes.
    Transferable info.data_dict_codon is saved.
    '''
    codon_list = nmer_maker()
    with open(in_genes, 'r') as fasta:
        seq = ''
        if spaces == 'no_keep':
            for line in fasta:
                if line.startswith(">"):
                    header = int(line[1:].strip("\n")) # no ">", no "\n"
                else:
                    seq = line.strip("\n").strip("*")
                    codon_counter(seq, header, codon_list)
            codon_counter(seq, header, codon_list) # include last one

        elif spaces == 'True':
            for line in fasta:
                if line.startswith(">"):
                    header = line[1:].strip("\n").replace("~(#)~"," ") # no ">", no "\n"
                    try:
                        header = header.split(" # ",1)[0]
                    except Exception:
                        pass
                    split_head = keep[header.rsplit("_",1)[0]][1]
                    length = len(seq)%3
                    if length == 1:
                        seq += '**'
                    elif length == 2:
                        seq += '*' # separate the sequences so 3-mers aren't connected
                    try:
                        if split_head != save:
                            codon_counter(seq, save, codon_list)
                            seq = ''
                    except Exception:
                        seq = ''
                        pass
                    save = split_head

                else:
                    seq += line.strip("\n").strip("*") # no ">", no "\n"
            codon_counter(seq, save, codon_list) # include last one

        else:
            for line in fasta:
                if line.startswith(">"):
                    header = line[1:].strip("\n") # no ">", no "\n"
                    try:
                        header = header.split(" # ",1)[0]
                    except Exception:
                        pass
                    split_head = keep[header.rsplit("_",1)[0]][1]
                    length = len(seq)%3
                    if length == 1:
                        seq += '**'
                    elif length == 2:
                        seq += '*' # separate the sequences so 3-mers aren't connected
                    try:
                        if split_head != save:
                            codon_counter(seq, save, codon_list)
                            seq = ''
                    except Exception:
                        seq = ''
                        pass
                    save = split_head

                else:
                    seq += line.strip("\n").strip("*") # no ">", no "\n"

            codon_counter(seq, save, codon_list) # include last one

    seq = None
    line = None

    with open(pickle_folder + file_n + "_pickle.sav", 'wb') as save_info:
         pickle.dump(info.data_dict_codon, save_info)


if __name__ == '__main__':
    pickle_folder = sys.argv[2]
    file_n = sys.argv[3]
    spaces = sys.argv[4]
    keep = None
    if spaces != 'no_keep':
        with open(pickle_folder + "keep_pickle.sav", 'rb') as read_keep:
            keep = pickle.load(read_keep)

    codon_stuff(sys.argv[1], pickle_folder, file_n, spaces, keep)

#
#
#
