#! /usr/bin/env python3
# vRhyme
# Author: Kristopher Kieft
# University of Wisconsin-Madison

import warnings
warnings.simplefilter("ignore")
import numpy as np
from numba import jit


@jit(nopython=True)
def cohenD(query, subject, cohen_threshold, cohen_norm, min_cov, cohen_range, max_cohen, penalty_w, penalty_n, values):
    '''
    Calculate Cohen's d values,
    add penalties if necessary,
    normalize and return coverage co-variance
    '''
    p = 0
    count_penalty = 0
    values_min = 0
    for i in range(0,cohen_range,2): # iterate over avgs/stdevs
        query_avg = query[i]
        sub_avg = subject[i]
        if sub_avg >= min_cov or query_avg >= min_cov:
            values_min += 1
            sub_sd = subject[i+1] 
            query_sd = query[i+1] 

            sub_sd += 0.1 # non-zero
            query_sd += 0.1
            # Cohen's d equation
            pool = ((query_sd**2+sub_sd**2)/2)**0.5
            d = abs((sub_avg-query_avg)/pool)

            if d <= cohen_threshold:
                if d > max_cohen:
                    if np.absolute(query_avg - sub_avg) <= 1 and (sub_avg or query_avg):
                        d = cohen_norm
                    else:
                        p += penalty_w
                        count_penalty += 1

                values[i] = d
            else:
                return False, 0
        
        # else: both are below min_cov, don't consider Cohen's d.
        # Should never have all comparisons be 0 due to pre-filtering

        if count_penalty > penalty_n:
            return False, 0

    values = values[~np.isnan(values)]
    l = values.size

    if count_penalty > l/2 or l < values_min-penalty_n:
        return False, 0

    if l:
        avg = np.nanmean(values)
        avg += p
        avg /= np.log(l+1)

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
    counter_vals = []
    full = []
    cohenout = []
    refine = max_edges*5
    with open(cov_table, 'r') as infile: 
        first_line = infile.readline().split("\t")
        cohen_check = int((len(first_line)-1)/2)
        cohen_range = cohen_check*2
        values = np.empty(cohen_range)
        values[:] = np.nan
        cohen_threshold = max_cohen*1.25
        cohen_norm = max_cohen*0.75
        first_line = None

        if cohen_check == 0:
            return None, None, None, cohen_check

        # can't have 1/3 samples be penalized, can bring penalty_n = 0
        while True:
            if penalty_n/cohen_check > 0.3:
                penalty_n -= 1
            else:
                break

        for line in infile:
            line = line.strip("\n").split("\t")
            genome = line[0]

            try:
                _, map_qry = keep[genome]
                query  = [float(x) for x in line[1:]] # 1: to skip names
                query_avgs = query[::2]
                query_avgs.sort(reverse = True)
                if query_avgs[0] >= min_cov:
                    query = np.array(query)
                    ##### Cohen's d
                    for map_sub,subject in cd_genome_dict.items():
                        check, avg = cohenD(query, subject, cohen_threshold, cohen_norm, min_cov, cohen_range, max_cohen, penalty_w, penalty_n, values)
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

    counter_vals = None
    cd_genome_dict = None
    pair_genomes = set(full)

    return pair_genomes, full, cohenout, cohen_check


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
