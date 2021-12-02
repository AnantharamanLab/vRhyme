#! /usr/bin/env python3
# vRhyme
# Author: Kristopher Kieft
# University of Wisconsin-Madison

import numpy as np
from numba import jit


@jit(nopython=True)
def cohenD(query, subject, cohen_threshold, cohen_norm, min_cov, cohen_range, cohen_check, max_cohen, penalty_w, penalty_n):
    '''
    Calculate Cohen's d values,
    add penalties if necessary,
    normalize and return coverage co-variance
    '''
    p = 0
    values = []
    for i in range(0,cohen_range,2): # iterate over avgs/stdevs
        query_avg = query[i]
        sub_avg = subject[i]
        sub_min = False
        query_min = False
        if sub_avg >= min_cov:
            sub_min = True
        if query_avg >= min_cov:
            query_min = True
        if sub_min or query_min:
            sub_sd = subject[i+1] 
            query_sd = query[i+1] 

            sub_sd += 0.1 # non-zero
            query_sd += 0.1
            # Cohen's d equation
            pool = ((query_sd**2+sub_sd**2)/2)**0.5
            d = abs((sub_avg-query_avg)/pool)

            if d <= cohen_threshold:
                if query_min == False:
                    p += penalty_w
                elif sub_min == False:
                    p += penalty_w
                if d > max_cohen:
                    p += penalty_w
                values.append(d)
            else:
                return False, 0
        # else: both are below min_cov, don't consider Cohen's d.
        # Should never have all comparisons be 0 due to pre-filtering

    try:
        p_penalty_w = p/penalty_w
    except Exception: # numba doesn't allow ZeroDivisionError
        p_penalty_w = 0

    if p_penalty_w <= penalty_n:
        values = np.array(values)
        avg = np.mean(values)

        avg += p
        avg /= np.log(values.size+1)

        if avg <= cohen_norm:
            return True, avg

    return False, 0


def coverage_table_stuff(cov_table, keep, max_cohen, min_cov, penalty_w, penalty_n, max_edges): 
    '''
    Parse coverage table and calculate co-variance by Cohen's d.
    pair_filter() also filters to retain only the best 'refine' hits
    '''
    cd_genome_dict = {}
    counter = {}
    full = []
    cohenout = []
    refine = max_edges*3
    with open(cov_table, 'r') as infile: 
        first_line = infile.readline().split("\t")
        cohen_check = int((len(first_line)-1)/2)
        cohen_range = cohen_check*2
        cohen_threshold = max_cohen*1.25
        cohen_norm = max_cohen*0.75
        first_line = None

        # can't have as many penalties as there are samples
        if penalty_n >= cohen_check:
            penalty_n = cohen_check-1

        for line in infile:
            line = line.strip("\n").split("\t")
            genome = line[0]

            try:
                _, map_qry = keep[genome]
                query  = [float(x) for x in line[1:]] # 1: to skip names
                query_avgs = query[::2]
                query_avgs.sort(reverse = True) # faster than max()?
                if query_avgs[0] >= min_cov:
                    query = np.array(query)
                    ##### Cohen's d
                    for map_sub,subject in cd_genome_dict.items():
                        check, avg = cohenD(query, subject, cohen_threshold, cohen_norm, min_cov, cohen_range, cohen_check, max_cohen, penalty_w, penalty_n)
                        if check:
                            counter = pair_filter(counter, map_qry, map_sub, avg, refine)
                    # add to dictionary after in order to not run a query:query check
                    cd_genome_dict[map_qry] = query

            except KeyError:
                # genome is not in keep, skip it
                pass

    counter_vals = set([item for sublist in counter.values() for item in sublist])
    for c in counter_vals:
        q,s,a = c
        full.extend([q,s])
        cohenout.append(a)

    cd_genome_dict = None
    pair_genomes = set(full)

    return pair_genomes, full, cohenout


def pair_filter(counter, map_qry, map_sub, avg, refine):
    '''
    Filter to only keep the top "refine" Cohen d values per scaffold.
    Each scaffold is considered individually. 
    More efficient than filtering at the end for larger datasets.
    '''
    q = counter.setdefault(map_qry, [])
    s = counter.setdefault(map_sub, [])
    q.append((map_qry,map_sub,avg))
    s.append((map_qry,map_sub,avg))
    q.sort(key = lambda x: x[2])
    s.sort(key = lambda x: x[2])
    counter[map_qry] = q[:refine]
    counter[map_sub] = s[:refine]

    return counter

#
#
#
