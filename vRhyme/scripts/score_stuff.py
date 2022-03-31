#! /usr/bin/env python3
# vRhyme
# Author: Kristopher Kieft
# University of Wisconsin-Madison


import os
import subprocess
from collections import Counter


def uniques_stuff(folder):
    '''
    Compile unique bins from all iterations
    (no need to score identical bins multiple times)
    '''
    bin_files = []
    temp_files = os.listdir(folder)
    for t in temp_files:
        if t.endswith('.vRhyme-bins.tsv'):
            bin_files.append(t)
    
    if len(bin_files) == 0:
        return None, None, None
    
    uniques = {}
    total_seqs = {}
    total_bins = {}
    bins = 0
    for b in bin_files:
        iteration = int(b.split('.',1)[0])
        with open(folder + b, 'r') as file:
            total = 0
            iter_bins = 0
            for line in file: # each line is a bin
                mems = [int(m) for m in line.strip('\n').split('\t')]
                mems.sort()
                uniques.setdefault(tuple(mems),[bins]).append(iteration)
                bins += 1
                iter_bins += 1

                total += len(mems)

        total_seqs[iteration] = total
        total_bins[iteration] = iter_bins
    
    return uniques, total_seqs, total_bins


def score_stuff(uniques, iterations, prots_dict, linclust_master, red_cutoff, total_bins, total_seqs):
    '''
    Score each unique bin
    '''
    redundancy_dict = {}
    proteins_per_bin = {}
    uniques_info = {}
    for i in range(0,iterations):
        redundancy_dict[i] = 0
        proteins_per_bin[i] = 0

    for u in uniques:
        value = uniques[u]
        identifier = value[0]
        iters = value[1:]

        proteins = 0
        clusters = []
        for mem in u:
            proteins += prots_dict.get(mem,0)
            clusters.extend(linclust_master.get(mem,[]))
        
        dup_dict = Counter(clusters)
        dup_list = [val for val in dup_dict.values() if val > 1]
        duplicated_members = sum(dup_list)-len(dup_list) # total redundant proteins - representative protein in each duplicated cluster

        if duplicated_members <= red_cutoff:
            uniques_info[identifier] = (proteins,duplicated_members)
            for i in iters:
                redundancy_dict[i] += duplicated_members
                proteins_per_bin[i] += proteins
        else:
            for i in iters:
                total_bins[i] -= 1
                total_seqs[i] -= len(u)
    
    return redundancy_dict, proteins_per_bin, uniques_info, total_bins, total_seqs



def best_bin(folder, iterations, total_seqs, redundancy_dict, total_bins, proteins_per_bin, length_keep):
    '''
    Identify best bin by score
    '''
    scoring = {}
    for i in range(0,iterations):
        ts = total_seqs.get(i, 0)
        rd = redundancy_dict.get(i, 0)
        tb = total_bins.get(i, 0)
        pb = proteins_per_bin.get(i, 0)
        scoring.update({i: (ts, rd, tb, pb)})

    best = []
    for key in scoring.keys():
        score = scoring[key]
        seqs = score[0]

        if seqs > 0:
            f = seqs/length_keep
            try:
                r = (2*(score[1]/score[3]))**0.5
            except Exception:
                r = 0 # no proteins

            b = (score[2]/seqs)**2
            s = round(f - b - (3*r),4)
        else:
            s = -1

        best.append((key,s))
    best.sort(key = lambda x: x[1], reverse=True)

    subprocess.run(f'mkdir {folder}vRhyme_alternate_bins', shell=True)
    with open(f'{folder}vRhyme_alternate_bins/vRhyme_bin_scoring.tsv','w')as outscore:
        outscore.write('iteration\tsequences\tredundancy\tbins\tproteins\tscore\n')
        for item in best:
            key = item[0]
            score = item[1]
            info = scoring[key]
            outscore.write(f'{key}\t{info[0]}\t{info[1]}\t{info[2]}\t{info[3]}\t{score}\n')

    return best, scoring


def final_bin(folder, best, scoring, mapper, unqiues, uniques_info):
    '''
    Write out summary files for all bins from all iterations
    '''
    final = best[0][0]
    score = best[0][1]
    info = scoring[final]
    binned_seqs = info[0]
    binned_proteins = info[3]
    binned_redundancy = info[1]
    bins_count = info[2]

    with open(f'{folder}{final}.vRhyme-bins.tsv', 'r') as infile, open(f'{folder}vRhyme_best_bins.{final}.membership.tsv', 'w') as outfile, open(f'{folder}vRhyme_best_bins.{final}.summary.tsv', 'w') as summary:
        outfile.write('scaffold\tbin\n')
        summary.write('bin\tmembers\tproteins\tredundancy\n')
        b = 1
        for line in infile:
            line = [int(l) for l in line.strip('\n').split('\t')]
            line.sort()
            line = tuple(line)

            try:
                identifier = unqiues[line][0]
                info = uniques_info[identifier]
                proteins = info[0]
                redundancy = info[1]
                summary.write(f'{b}\t{len(line)}\t{proteins}\t{redundancy}\n')

                for l in line:
                    outfile.write(f'{mapper[l]}\t{b}\n')

                b += 1
            except Exception:
                pass # bin was removed due to redundancy

    subprocess.run(f'rm {folder}{final}.vRhyme-bins.tsv', shell=True)

    for item in best[1:]:
        item = item[0]

        if not os.path.exists(f'{folder}{item}.vRhyme-bins.tsv'):
            continue

        with open(f'{folder}{item}.vRhyme-bins.tsv', 'r') as infile, open(f'{folder}vRhyme_alternate_bins/vRhyme_alternate_bins.{item}.membership.tsv', 'w') as outfile, open(f'{folder}vRhyme_alternate_bins/vRhyme_alternate_bins.{item}.summary.tsv', 'w') as summary:
            outfile.write('scaffold\tbin\n')
            summary.write('bin\tmembers\tproteins\tredundancy\n')
            b = 1
            for line in infile:
                line = [int(l) for l in line.strip('\n').split('\t')]
                line.sort()
                line = tuple(line)
                try:
                    identifier = unqiues[line][0]
                    info = uniques_info[identifier]
                    proteins = info[0]
                    redundancy = info[1]
                    summary.write(f'{b}\t{len(line)}\t{proteins}\t{redundancy}\n')

                    for l in line:
                        outfile.write(f'{mapper[l]}\t{b}\n')

                    b += 1
                except Exception:
                    pass # bin was removed due to redundancy

        subprocess.run(f'rm {folder}{item}.vRhyme-bins.tsv', shell=True)
    
    return final, binned_seqs, bins_count, score, binned_proteins, binned_redundancy


#
#
#
