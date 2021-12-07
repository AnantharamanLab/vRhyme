#! /usr/bin/env python3
# vRhyme
# Author: Kristopher Kieft
# University of Wisconsin-Madison

from numba import jit
try:
    from scripts import save_dist_dicts as info
except (ModuleNotFoundError, ImportError):
    import save_dist_dicts as info


@jit(nopython=True)
def subtraction_math(a, b):
    '''
    Pairwise distance by differences
    '''
    c = abs(a-b)
    return c

@jit(nopython=True)
def comparing_gc(compare_gc, max_gc):
    '''
    Filter for GC difference maximum
    '''
    gc = True
    if compare_gc >= max_gc:
        gc = False
    return gc

@jit(nopython=True)
def comparing_kmer(compare_nt, min_kmer):
    '''
    Filter for nucleotide similarity minimum
    '''
    kmer = False
    if compare_nt >= min_kmer:
        kmer = True
    return kmer


@jit(nopython=True)
def dist_compare(first, second, denomer_f, denomer_s, length):
    '''
    Cosine similarity equation
    '''
    numer = 0
    for n in range(0,length):
        numer += (first[n])*(second[n])
    denomer = (denomer_f*denomer_s)**0.5
    value = numer/denomer
    return value


def distance_pairs(folder, base, pairs, cohen_list, max_gc, min_kmer, mapper):
    '''
    Main function for all pairwise distance calculations
    '''
    pairs_machine = []
    cohen_machine = []
    length_codon = 64
    length_nt = 136
    with open(folder + base + '.distances.tsv' ,'w') as outfile, open(folder + 'log_vRhyme_distance_names.tsv' ,'w') as namesfiles: 
        outfile.write('gc_dist\tcpg_dist\tgc_skew_dist\tZOM_dist\tTUD_dist\tnt_dist\tcodon_dist')
        namesfiles.write('query\tsubject')
        for n in range(0,len(pairs),2):
            p1 = pairs[n]
            p2 = pairs[n+1]

            compare_gc = subtraction_math(info.data_dict_gc[p1], info.data_dict_gc[p2])
            gc = comparing_gc(compare_gc, max_gc)

            if gc == True:
                n1 = info.data_dict_nt[p1]
                n1_data = n1[0]
                n1_denom = n1[1]
                n2 = info.data_dict_nt[p2]
                n2_data = n2[0]
                n2_denom = n2[1]

                compare_nt = dist_compare(n1_data, n2_data, n1_denom, n2_denom, length_nt)
                kmer = comparing_kmer(compare_nt, min_kmer)

                if kmer == True:
                    z1 = info.data_dict_zom[p1]
                    z1_data = z1[0]
                    z1_denom = z1[1]
                    z2 = info.data_dict_zom[p2]
                    z2_data = z2[0]
                    z2_denom = z2[1]

                    t1 = info.data_dict_tud[p1]
                    t1_data = t1[0]
                    t1_denom = t1[1]
                    t2 = info.data_dict_tud[p2]
                    t2_data = t2[0]
                    t2_denom = t2[1]

                    compare_zom = dist_compare(z1_data, z2_data, z1_denom, z2_denom, length_nt)
                    compare_tud = dist_compare(t1_data, t2_data, t1_denom, t2_denom, length_nt)

                    compare_cpg = subtraction_math(info.data_dict_cpg[p1], info.data_dict_cpg[p2])
                    compare_skew = subtraction_math(info.data_dict_skew[p1], info.data_dict_skew[p2])
                    try:
                        c1 = info.data_dict_codon[p1]
                        c1_data = c1[0]
                        c1_denom = c1[1]
                        c2 = info.data_dict_codon[p2]
                        c2_data = c2[0]
                        c2_denom = c2[1]
                        compare_codon = dist_compare(c1_data, c2_data, c1_denom, c2_denom, length_codon)

                    except Exception:
                        compare_codon = 0.75 # arbitrary value when no proteins

                    outfile.write('\n' + str(compare_gc) + '\t' + str(compare_cpg) + '\t' + str(compare_skew) + '\t' + str(compare_zom) + '\t' + str(compare_tud) + '\t' + str(compare_nt) + '\t' + str(compare_codon))
                    pairs_machine.extend([p1,p2])
                    namesfiles.write('\n' + str(mapper[p1]) + "\t" + str(mapper[p2]))
                    cohen_machine.append(cohen_list[int(n/2)])

    cohen_list = None
    return folder + base + '.distances.tsv', pairs_machine, cohen_machine
#
#
#
