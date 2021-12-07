#! /usr/bin/env python3
# vRhyme
# Author: Kristopher Kieft
# University of Wisconsin-Madison

import warnings
warnings.simplefilter("ignore")
import subprocess
import os
import numpy as np
from numba import jit
import pickle
import sys
import pysam


@jit(nopython=True)
def quick_stats(depth):
    '''
    Average and standard deviation
    '''
    a = 0
    s = 0
    a = np.mean(depth)
    s = np.std(depth)
    return a,s

@jit(nopython=True)
def add_depth(depth, start, end, ed, rl, read_id):
    '''
    Add depth to base if it passes the read identity threshold
    '''
    if ed/rl <= read_id:
        for i in range(start,end):
            depth[i] += 1
    return depth

@jit(nopython=True)
def add_depth_no_ed(depth, start, end):
    '''
    Add depth to base
    '''
    for i in range(start,end):
        depth[i] += 1
    return depth

def coverage_stuff(alignment, mask, keep, folder, read_id):
    '''
    Read BAM files, sort and index if necessary,
    and main wrapper for calculating coverage per scaffold,
    generate main coverage table (compatible with -c)
    '''
    if mask != None:
        mask_f = mask
        mask_r = -mask
    else:
        mask_f = None
        mask_r = None
    f = alignment.rsplit("/",1)[0]
    if int(os.stat(alignment).st_size) == 0:
        subprocess.run("echo '" + f + "' >> " + folder + "log_vRhyme_bam_noalign.txt", shell=True)
        exit()
    check_aligned = subprocess.check_output(f"samtools view -@ 1 {alignment} | head -n 1", shell=True)
    if len(check_aligned) == 0:
        subprocess.run("echo '" + f + "' >> " + folder + "log_vRhyme_bam_noalign.txt", shell=True)
        exit()

    try:
        temp = alignment.rsplit("/",1)[1]
        base = temp.rsplit(".",1)[0]
    except Exception:
        base = alignment.rsplit(".",1)[0]

    sort_check = False
    try:
        check = subprocess.check_output("samtools view -@ 1 -H " + alignment + " | grep '@HD'", shell=True)
        if "coordinate" in str(check):
            sort_check = True
    except Exception:
        # no @HD line, retain sort_check = False
        pass
    if sort_check == False:
        subprocess.run("samtools sort -@ 1 -o " + folder + "vRhyme_bam_files/" + base + ".sorted.bam " + alignment + " 2> /dev/null", shell=True)
        alignment = folder + "vRhyme_bam_files/" + base + ".sorted.bam"

    master = {}

    if not os.path.exists(alignment + '.bai'):
        subprocess.run('samtools index ' + alignment, shell=True)

    bamfile = pysam.AlignmentFile(alignment, "rb")

    if read_id != 0: 
        for x in bamfile.fetch(until_eof=True):
            genome = x.reference_name
            length = keep.get(genome,[None])[0]
            if length: # not in keep
                try:
                    if genome != prev:
                        avg, sd = quick_stats(depth)
                        master[prev] = (avg,sd)
                        depth = np.zeros(length)
                except NameError:
                    # prev not defined
                    depth = np.zeros(length)

                prev = genome

                ed = 0
                for t in x.tags:
                    if t[0] == 'NM':
                        ed = t[1]
                        break
                rl = x.query_length

                start = x.reference_start
                end = x.reference_end
                if end:
                    depth = add_depth(depth, start, end, ed, rl, read_id)

        if length: # not in keep
            # last one
            depth = depth[mask_f:mask_r]
            avg, sd = quick_stats(depth)
            master[prev] = (avg,sd) # will break if there's only a single read aligned
            depth = None
    else: # no mismatches
        for x in bamfile.fetch(until_eof=True):
            genome = x.reference_name
            length = keep.get(genome,[None])[0]
            if length != None: # not in keep
                try:
                    if genome != prev:
                        depth = depth[mask_f:mask_r]
                        avg, sd = quick_stats(depth)
                        master[prev] = (avg,sd)
                        depth = np.zeros(length)
                except NameError:
                    depth = np.zeros(length)

                prev = genome
                start = x.reference_start
                end = x.reference_end
                if end:
                    depth = add_depth_no_ed(depth, start, end)

        if length != None: # not in keep
            # last one
            depth = depth[mask_f:mask_r]
            avg, sd = quick_stats(depth)
            master[prev] = (avg,sd) # will break if there's only a single read aligned
            depth = None

    bamfile.close()

    with open(folder + "vRhyme_coverage_files/" + base + ".coverage.tsv", "w") as cov_file:
        cov_file.write("avg_" + base + "\tstdev_" + base)
        for name in keep.keys():
            try:
                stats = master[name]
                cov_file.write("\n" + str(stats[0]) + "\t" + str(stats[1]))
            except Exception:
                cov_file.write("\n0\t0")


if __name__ == '__main__':
    mask = int(sys.argv[2])
    if mask == 0:
        mask = None
    folder = sys.argv[3]
    read_id = float(sys.argv[4])
    with open(folder + "keep_pickle.sav", 'rb') as read_keep:
        keep = pickle.load(read_keep)

    coverage_stuff(sys.argv[1], mask, keep, folder, read_id)

#
#
#
