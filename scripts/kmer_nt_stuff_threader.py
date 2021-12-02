#! /usr/bin/env python3
# vRhyme
# Author: Kristopher Kieft
# University of Wisconsin-Madison

import sys
import collections
import itertools
from numba import jit
import numpy as np
import save_dist_dicts as info
import pickle

@jit(nopython=True)
def denomer_generate(len_nt_list, in_list, in_avg):
    '''
    Denominator for pairwise cosine distance calculations.
    Range is for all 4mers.
    '''
    denomer = 0
    for n in range(0,len_nt_list):
        denomer += (in_list[n]-in_avg)**2
    return denomer


@jit(nopython=True)
def zom_generate(i1, i2, i3, i4, value):
    '''
    Generate ZOM pairwise distance calculations.
    '''
    zom = value/(i1*i2*i3*i4)
    return zom

@jit(nopython=True)
def expt_generate(ai,ti,ci,gi, a,t,c,g, total, value):
    '''
    Generate TUD pairwise distance calculations.
    '''
    expt_a = (ai**(2*a))
    expt_t = (ti**(2*t))
    expt_c = (ci**(2*c))
    expt_g = (gi**(2*g))
    expt = expt_a * expt_t * expt_c * expt_g * total
    tud = value*expt
    return tud

@jit(nopython=True)
def usage_generate(a,t,c,g, cpg_count):
    '''
    Calculate useful nucleotide information 
    '''
    g_plus_c = g+c
    den = a+t+g_plus_c
    GC = g_plus_c/den
    skew = (g-c)/g_plus_c
    cpg = cpg_count/den
    ai = a/den
    ti = t/den
    gi = g/den
    ci = c/den
    return GC, cpg, skew, ai, ti, gi, ci


@jit(nopython=True)
def array_generate(datalist, avg, len_nt_list):
    '''
    Begin calculation of numerator for cosine pairwise distance calculations.
    Range is for all 4mers.
    '''
    for n in range(0,len_nt_list):
        datalist[n] -= avg
    return datalist

def kmer_nt_counter(seq, header, length, nt_list, len_nt_list):
    '''
    Count occurances of 4mers within scaffolds and
    calculate numerator/denominator values for difference/cosine/ZOM/TUD pairwise distance calculations.
    Store information within transferable info.* dictionaries
    '''
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    r_seq = ''
    seq = seq.upper()

    try:
        r_seq = ''.join([complement[base] for base in seq[::-1]])
    except Exception as e:
        pass

    # counts
    g = seq.count("G")
    c = seq.count("C")
    a = seq.count("A")
    t = seq.count("T")

    cpg_count = seq.count("CG")

    GC, cpg, skew, ai, ti, gi, ci = usage_generate(a,t,c,g, cpg_count)
    composition_dict = {'A':ai, 'T':ti, 'C':ci, 'G':gi}

    info.data_dict_gc.update({header:GC})
    info.data_dict_cpg.update({header:cpg})
    info.data_dict_skew.update({header:skew})

    # tetranucleotide
    nmer_len = 4
    step = 1
    kmer_list = [seq[i:i+nmer_len] for i in range(0,length-nmer_len+1,step)] # forward
    kmer_r_list = [r_seq[i:i+nmer_len] for i in range(0,length-nmer_len+1,step)] # reverse complement
    seq = None
    r_seq = None

    master_list = kmer_list + kmer_r_list
    kmer_list = None
    kmer_r_list = None

    nt_dict = dict(collections.Counter(master_list)) # count kmers

    total = sum(nt_dict.values())
    tud_list = []
    zom_list = []
    kmer_list = []
    for i in range(0,len_nt_list):
        item = nt_list[i]
        try:
            check = nt_dict[item] # except if not in dict
            value = check/total
            kmer_list.append(value)

            g = item.count("G")
            c = item.count("C")
            a = item.count("A")
            t = item.count("T")

            i1 = composition_dict[item[0]]
            i2 = composition_dict[item[1]]
            i3 = composition_dict[item[2]]
            i4 = composition_dict[item[3]]

            zom = zom_generate(i1,i2,i3,i4, value)
            tud = expt_generate(ai,ti,ci,gi, a,t,c,g, total, value)
            tud_list.append(tud)
            zom_list.append(zom)

        except Exception:
            kmer_list.append(0)
            tud_list.append(0)
            zom_list.append(0)

    kmer_list = np.array(kmer_list)
    tud_list = np.array(tud_list)
    zom_list = np.array(zom_list)

    z_avg = np.mean(zom_list)
    denomer = denomer_generate(len_nt_list,zom_list, z_avg)
    zom_list = array_generate(zom_list, z_avg, len_nt_list)
    info.data_dict_zom.update({header:(zom_list,denomer)})

    t_avg = np.mean(tud_list)
    denomer = denomer_generate(len_nt_list,tud_list, t_avg)
    tud_list = array_generate(tud_list, t_avg, len_nt_list)
    info.data_dict_tud.update({header:(tud_list,denomer)})

    n_avg = np.mean(kmer_list)
    denomer = denomer_generate(len_nt_list,kmer_list, n_avg)
    kmer_list = array_generate(kmer_list, n_avg, len_nt_list)
    info.data_dict_nt.update({header:(kmer_list,denomer)})

    nt_dict = None

def nmer_maker():
    '''
    All possible 4mers for ATCG
    '''
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    nmer_len = 4
    nt_list = []
    temp_list = [''.join(mer) for mer in itertools.product('ATCG', repeat=nmer_len)]
    for item in temp_list:
        rc = ''.join([complement[base] for base in item[::-1]])
        if item not in nt_list and rc not in nt_list:
            nt_list.append(item)

    len_nt_list = len(nt_list)
    return nt_list,len_nt_list 


def kmer_nt_stuff(in_fasta, pickle_folder, file_n, keep):
    '''
    Wrapper for all functions, including reading scaffolds.
    '''
    nt_list, len_nt_list = nmer_maker()
    with open(in_fasta, 'r') as fasta:
        if keep == None:
            seq = ''
            for line in fasta:
                if line.startswith(">"):
                    length = len(seq)
                    if seq != '':
                        kmer_nt_counter(seq, header, length, nt_list, len_nt_list)
                    seq = ''
                    header = int(line[1:].strip("\n")) # no ">", no "\n"
                else:
                    seq += line.strip("\n")

            length = len(seq)
            kmer_nt_counter(seq, header, length, nt_list, len_nt_list) # include last one
        else:
            seq = ''
            for line in fasta:
                if line.startswith(">"):
                    length = len(seq)
                    if seq != '':
                        kmer_nt_counter(seq, header, length, nt_list, len_nt_list)
                    seq = ''
                    header = int(keep[line[1:].strip("\n")][1]) # no ">", no "\n"
                else:
                    seq += line.strip("\n")

            length = len(seq)
            kmer_nt_counter(seq, header, length, nt_list, len_nt_list) # include last one

    line = None
    seq = None
    nt_list = None

    with open(pickle_folder + file_n + "_pickle.sav", 'wb') as save_info:
         pickle.dump((info.data_dict_gc, info.data_dict_cpg, info.data_dict_skew, info.data_dict_zom, info.data_dict_tud, info.data_dict_nt), save_info)


if __name__ == '__main__':
    pickle_folder = sys.argv[2]
    file_n = sys.argv[3]
    if sys.argv[4] == 'False':
        kmer_nt_stuff(sys.argv[1], pickle_folder, file_n, None)
    else:
        with open(pickle_folder + "keep_pickle.sav", 'rb') as read_keep:
            keep = pickle.load(read_keep)
        kmer_nt_stuff(sys.argv[1], pickle_folder, file_n, keep)

#
#
#
