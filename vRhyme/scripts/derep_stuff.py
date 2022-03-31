#! /usr/bin/env python3
# vRhyme
# Author: Kristopher Kieft
# University of Wisconsin-Madison

import subprocess
import pandas as pd
import os
from numba import jit
import time

@jit(nopython=True)
def count_check(a,t,c,g,x,l):
    '''
    Counting if there are any non-ATCG nucleotides.
    Could just count Ns, but this makes sure there isn't anything else.
    '''
    check = False
    if a+t+c+g+x != l:
        check = True
    return check

def split_singles_stuff(in_fasta, base, derep_len, folder):
    '''
    Read in the fasta file,
    check for ambiguous bases and sequence lengths,
    split into new files with one seq per file
    '''
    letters = set(['A','T','C','G','N'])
    subprocess.run("rm -R " + folder + "vRhyme_dereplication/split_fasta_files 2> /dev/null", shell=True)
    time.sleep(0.1)
    subprocess.run("mkdir " + folder + "vRhyme_dereplication/split_fasta_files", shell=True)
    outfolder = folder + "vRhyme_dereplication/split_fasta_files/"
    seqs_dict = {}
    n = 1

    with open(in_fasta, 'r') as fasta:
        seq = ''
        for line in fasta:
            if line.startswith(">"):
                l = len(seq)
                if l >= derep_len:
                    seq = seq.upper()
                    a = seq.count('A')
                    t = seq.count('T')
                    c = seq.count('C')
                    g = seq.count('G')
                    x = seq.count('N')
                    check = count_check(a,t,c,g,x,l)
                    if check == True:
                        new_seq = ''
                        for base in seq:
                            if base in letters:
                                new_seq += base
                            else:
                                new_seq += 'N'

                        seq = new_seq
                    with open(outfolder + str(n) + ".fa", 'w') as outfasta:
                        outfasta.write(">" + header + "\n" + seq + "\n")
                        seqs_dict.update({header:seq})
                        seq = ''
                        n += 1

                seq = ''
                header = line[1:].strip("\n").replace(" ","~#~") # no ">", no "\n"
            else:
                seq += line.strip("\n")
        # count last seq
        l = len(seq)
        if l >= derep_len:
            seq = seq.upper()
            a = seq.count('A')
            t = seq.count('T')
            c = seq.count('C')
            g = seq.count('G')
            x = seq.count('N')
            check = count_check(a,t,c,g,x,l)
            if check == True:
                new_seq = ''
                for base in seq:
                    if base in letters:
                        new_seq += base
                    else:
                        new_seq += 'N'

                seq = new_seq
            with open(outfolder + str(n) + ".fa", 'w') as outfasta:
                outfasta.write(">" + header + "\n" + seq + "\n")
                seqs_dict.update({header:seq})
                seq = ''
                n += 1

    n -= 1
    return seqs_dict, n


def split_singles_stuff_interest(in_fasta, base, derep_len, folder):
    '''
    Read in the fasta file,
    make sure sequence name is in the interest list,
    check for ambiguous bases and sequence lengths,
    split into new files with one seq per file
    '''
    letters = set(['A','T','C','G','N'])
    subprocess.run("rm -R " + folder + "vRhyme_dereplication/split_fasta_files 2> /dev/null", shell=True)
    time.sleep(0.1)
    subprocess.run("mkdir " + folder + "vRhyme_dereplication/split_fasta_files", shell=True)
    outfolder = folder + "vRhyme_dereplication/split_fasta_files/"
    seqs_dict = {}
    n = 1

    with open(in_fasta, 'r') as fasta:
        seq = ''
        for line in fasta:
            if line.startswith(">"):
                l = len(seq)
                if l >= derep_len:
                    seq = seq.upper()
                    a = seq.count('A')
                    t = seq.count('T')
                    c = seq.count('C')
                    g = seq.count('G')
                    x = seq.count('N')
                    check = count_check(a,t,c,g,x,l)
                    if check == True:
                        new_seq = ''
                        for base in seq:
                            if base in letters:
                                new_seq += base
                            else:
                                new_seq += 'N'

                        seq = new_seq
                    with open(outfolder + str(n) + ".fa", 'w') as outfasta:
                        outfasta.write(">" + header + "\n" + seq + "\n")
                        seqs_dict.update({header:seq})
                        seq = ''
                        n += 1

                seq = ''
                header = line[1:].strip("\n").replace(" ","~#~") # no ">", no "\n"
            else:
                seq += line.strip("\n")
        l = len(seq)
        if l >= derep_len:
            seq = seq.upper()
            a = seq.count('A')
            t = seq.count('T')
            c = seq.count('C')
            g = seq.count('G')
            x = seq.count('N')
            check = count_check(a,t,c,g,x,l)
            if check == True:
                new_seq = ''
                for base in seq:
                    if base in letters:
                        new_seq += base
                    else:
                        new_seq += 'N'

                seq = new_seq
            with open(outfolder + str(n) + ".fa", 'w') as outfasta:
                outfasta.write(">" + header + "\n" + seq + "\n")
                seqs_dict.update({header:seq})
                seq = ''
                n += 1

    n -= 1
    return seqs_dict, n


def mash_stuff(base, threads, total, folder, sketch_k, sketch_s, iterations):
    '''
    Run Mash sketch, paste and dist
    '''
    t = int(int(threads)/int(iterations))
    if t < 1:
        t = 1
    if int(threads) < int(iterations):
        loop = threads
    else:
        loop = iterations

    time.sleep(0.1)
    subprocess.run("rm -R " + folder + "vRhyme_dereplication/mash_sketch_files 2> /dev/null", shell=True)
    time.sleep(0.1)
    subprocess.run("mkdir " + folder + "vRhyme_dereplication/mash_sketch_files", shell=True)
    i = 0
    s_list = []
    for n in range(1,total+1):
        s = subprocess.Popen("mash sketch -k " + str(sketch_k) + " -s " + str(sketch_s) + " -p " + str(t) + " -S 44 -o " + folder + "vRhyme_dereplication/mash_sketch_files/mash_sketch_" + str(n) + " " + folder + "vRhyme_dereplication/split_fasta_files/" + str(n) + ".fa 2> /dev/null", shell=True)
        s_list.append(s)
        i += 1
        if i == threads:
            for item in s_list:
                item.wait()
            s_list = []
            i = 0
    for item in s_list:
        item.wait()
    time.sleep(0.1)

    i = 0
    s_list = []
    for n in range(0,10):
        s = subprocess.Popen(f"mash paste {folder}vRhyme_dereplication/mash_sketch_temp_{n} {folder}vRhyme_dereplication/mash_sketch_files/mash_sketch_*{n}.msh 2> /dev/null", shell=True)
        s_list.append(s)
        i += 1
        if i == threads:
            for item in s_list:
                item.wait()
            s_list = []
            i = 0
    for item in s_list:
        item.wait()

    subprocess.run(f"mash paste {folder}vRhyme_dereplication/mash_sketch_combined {folder}vRhyme_dereplication/mash_sketch_temp* 2> /dev/null", shell=True)
    subprocess.run(f"rm {folder}vRhyme_dereplication/mash_sketch_temp* 2> /dev/null", shell=True) 

    time.sleep(0.1)
    #
    i = 0
    s_list = []
    for n in range(0,10):

        s = subprocess.Popen("mash dist -p " + str(t) + " -S 44 " + folder + "vRhyme_dereplication/mash_sketch_combined.msh " + folder + "vRhyme_dereplication/split_fasta_files/*" + str(n) + ".fa > " + folder + "vRhyme_dereplication/mash_dist_" + base + ".run_" + str(n) + ".tsv 2> /dev/null", shell=True) 
        s_list.append(s)
        i += 1
        if i == loop:
            for item in s_list:
                item.wait()
            s_list = []
            i = 0
    for item in s_list:
        item.wait()

    subprocess.run("cat " + folder + "vRhyme_dereplication/mash_dist_" + base + ".run_*.tsv > " + folder + "vRhyme_dereplication/mash_dist_" + base + ".tsv", shell=True)
    subprocess.run("rm -R " + folder + "vRhyme_dereplication/mash_sketch_files 2> /dev/null", shell=True)
    subprocess.Popen("rm " + folder + "vRhyme_dereplication/mash_dist_" + base + ".run_*.tsv 2> /dev/null", shell=True)
    subprocess.Popen("rm " + folder + "vRhyme_dereplication/mash_sketch_combined.msh 2> /dev/null", shell=True)
    time.sleep(0.1)

def mash_parse_stuff(base, folder):
    '''
    Parse Mash dist outputs
    '''
    with open(folder + "vRhyme_dereplication/mash_dist_" + base + ".tsv", "r") as infile, open(folder + "vRhyme_dereplication/mash_dist_" + base + ".temp.tsv", "w") as out_data:
        for line in infile:
            line = line.strip("\n").split("\t")
            if line[3] == "0" and line[0] != line[1]:
                sorted = [line[0], line[1]]
                sorted.sort()
                out_data.write(sorted[0] + "\t" + sorted[1] + "\n")
    try:
        table = pd.read_csv(folder + "vRhyme_dereplication/mash_dist_" + base + ".temp.tsv", header=None, sep="\t")
        table = table.drop_duplicates()
        table.to_csv(folder + "vRhyme_dereplication/mash_dist_" + base + ".parsed.tsv", index=False, header=False, sep="\t")
        time.sleep(0.1)
        return True
    except Exception:
        return False

def nucmer_stuff(base, threads, nuc_c, nuc_b, nuc_g, folder):
    '''
    Run Nucmer
    '''
    subprocess.run("rm -R " + folder + "vRhyme_dereplication/nucmer_delta_files 2> /dev/null", shell=True)
    time.sleep(0.1)
    subprocess.run("mkdir " + folder + "vRhyme_dereplication/nucmer_delta_files", shell=True)
    outfolder = folder + "vRhyme_dereplication/nucmer_delta_files/"
    with open(folder + "vRhyme_dereplication/mash_dist_" + base + ".parsed.tsv", "r") as infile:
        n = 1
        i = 0
        s_list = []
        for line in infile:
            line = line.strip("\n").split("\t")
            split1 = line[0]
            split2 = line[1]
            s = subprocess.Popen('nucmer -c ' + str(nuc_c) + ' -b ' + str(nuc_b) + ' -g ' + str(nuc_g) + ' -p ' + outfolder + 'nucmer_' + str(n) + ' ' + str(split1) + ' ' + str(split2) + ' 2> /dev/null', shell=True)
            i += 1
            if i == threads:
                for item in s_list:
                    item.wait()
                s_list = []
                i = 0
            n += 1
    for item in s_list:
        item.wait()

    time.sleep(0.1)
    subprocess.run("rm " + folder + "vRhyme_dereplication/mash_dist_" + base + ".parsed.tsv", shell=True)

def nucmer_parse_stuff(threads, derep_id, derep_len, folder):
    '''
    Parse Nucmer outputs
    '''
    time.sleep(0.1)
    delta_files = os.listdir(folder + "vRhyme_dereplication/nucmer_delta_files/")
    i = 0
    s_list = []
    rm_list = []
    for delta in delta_files:
        s = subprocess.Popen("show-coords -d -c -H -I " + str(derep_id*100) + " -l -L " + str(derep_len) + " -T " + folder + "vRhyme_dereplication/nucmer_delta_files/" + delta + " >> " + folder + "vRhyme_dereplication/vRhyme_nucmer_delta_parsed.tsv 2> /dev/null", shell=True)
        s_list.append(s)
        rm_list.append(delta)
        i += 1
        if i == threads:
            for item in s_list:
                item.wait()
            s_list = []
            i = 0
    for item in s_list:
        item.wait()

    time.sleep(0.1)
    subprocess.run("rm -R " + folder + "vRhyme_dereplication/nucmer_delta_files/ 2> /dev/null", shell=True)
    #

# first dereplicate everything with 100% coverage alignments, keep longest
def nucmer_filter_1_stuff(original_dict, folder):
    '''
    Filter Nucmer outputs for identical sequences
    '''
    time.sleep(0.1)
    with open(folder + "vRhyme_dereplication/vRhyme_nucmer_delta_parsed.tsv", "r") as coords:
        len_dict = {}
        len_list = []
        bins = {}
        b = 1
        for line in coords:
            line = line.strip("\n").split("\t")
            if float(line[9]) >= 99.9 or float(line[10]) >= 99.9: # calling 99.9% as identical
                qry = line[13]
                sbj = line[14]
                len_dict.update({qry:int(line[7])})
                len_dict.update({sbj:int(line[8])})
                added = False
                for key in bins:
                    value = bins[key]
                    if qry in value or sbj in value:
                        value.update([qry,sbj])
                        bins[key] = value
                        added = True
                        break
                if added == False:
                    bins[b] = set([qry,sbj])
                    b += 1

        while True:
            names = list(bins.keys())
            l = len(names)
            alter = False
            for i in range(0,l-1):
                for j in range(i+1,l):
                    x = names[i]
                    y = names[j]
                    try:
                        x_values = bins[x]
                        y_values = bins[y]

                        if x_values & y_values != set():
                            x_values = x_values.union(y_values)
                            bins[x] = x_values

                            del bins[y]
                            alter = True

                    except Exception:
                        # y was deleted
                        pass

            if alter == False: # nothing was combined
                break

        names = None

        clust_dict = {}
        remove_list = []
        for names in bins.values():
            names = list(names)
            len_list = [len_dict[i] for i in names]
            z = list(zip(names,len_list))
            z.sort(key= lambda x: x[1], reverse=True)         
            name = z[0][0]
            clust_dict.update({name:names})
        
            for i,_ in z[1:]:
                remove_list.append(i)
                del original_dict[i]

    bins = None

    remove_list = set(remove_list)
    with open(folder + "vRhyme_dereplication/vRhyme_nucmer_delta_parsed.tsv", "r") as coords, open(folder + "vRhyme_dereplication/vRhyme_nucmer_delta_filtered.tsv", "w") as filtered:
        filtered.write('S1\tE1\tS2\tE2\tLen1\tLen2\tpID\tLenQ\tLenS\tCovQ\tCovS\tFrm1\tFrm2\tQry\tSbj')
        coords_check = False
        for line in coords:
            check = line.strip("\n").split("\t")
            if check[13] not in remove_list and check[14] not in remove_list:
                filtered.write("\n" + line.strip("\n"))
                coords_check = True
    
    return clust_dict, original_dict, coords_check


@jit(nopython=True)
def filter_1(gap1, gap2, cov1, cov2, o, f15,f1,f17,f3,f9,f24,f10,f25,f4,f19,f6,f21):
    '''
    Filtering of multiple Nucmer hits per genome pairs;
    deals with gaps between hits.
    Co-directional alignments
    '''
    gap1 += f15 - f1
    gap2 += f17 - f3
    cov1 += f9 + f24
    cov2 += f10 + f25
    o += f4 + f19
    id = (f6 + f21)/2
    return gap1, gap2, cov1, cov2, o, id

@jit(nopython=True)
def filter_2(gap1, gap2, cov1, cov2, o, f15,f1,f18,f2,f9,f24,f10,f25,f4,f19,f6,f21):
    '''
    Filtering of multiple Nucmer hits per genome pairs;
    deals with gaps between hits.
    Anti-directional alignments
    '''
    gap1 += f15 - f1
    gap2 += f18 - f2
    cov1 += f9 + f24
    cov2 += f10 + f25
    o += f4 + f19
    id = (f6 + f21)/2
    return gap1, gap2, cov1, cov2, o, id


def nucmer_filter_2_stuff(iteration, nuc_split, merge_cov, nuc_g, folder):
    '''
    The main Nucmer parser for calculating identity and coverage of hits
    '''
    time.sleep(0.1)
    if iteration == 0:
        infile = folder + "vRhyme_dereplication/vRhyme_nucmer_delta_filtered.tsv"
    elif iteration > 0:
        infile = folder + "vRhyme_dereplication/vRhyme_nucmer_delta_parsed.tsv"

    with open(infile, "r") as filtered, open(folder + "vRhyme_dereplication/vRhyme_nucmer_delta_merged.tsv", "w") as outfile:
        outfile.write('S1\tE1\tS2\tE2\tLen1\tLen2\tpID\tLenQ\tLenS\tCovQ\tCovS\tFrm1\tFrm2\tQry\tSbj')
        filtered = filtered.read().replace("\n","\t").split("\t")

        if filtered[-1] == '':
            filtered.extend(['', '', '', '', '', '', '', '', '', '', '', '', '', ''])
        else:
            filtered.extend(['', '', '', '', '', '', '', '', '', '', '', '', '', '', ''])
        occur = 1
        gap1 = 0
        gap2 = 0
        cov1 = 0
        cov2 = 0
        o = 0

        if iteration == 0:
            z = 15
        else:
            z = 0

        for n in range(z, len(filtered)-15,15):
            qry = filtered[n+13]
            sbj = filtered[n+14]
            f28 = filtered[n+28]
            if (qry,sbj) == (filtered[n+28],filtered[n+29]) and f28 != '':
                occur += 1

                f15 = int(filtered[n+15])
                f1 = int(filtered[n+1])
                f9 = float(filtered[n+9])
                f24 = float(filtered[n+24])
                f10 = float(filtered[n+10])
                f25 = float(filtered[n+25])
                f4 = int(filtered[n+4])
                f19 = int(filtered[n+19])
                f6 = float(filtered[n+6])
                f21 = float(filtered[n+21])

                if filtered[n+12] == '1':
                    f17 = int(filtered[n+17])
                    f3 = int(filtered[n+3])
                    frm = 1
                    if occur == 2:
                        s1 = filtered[n]
                        s2 = filtered[n+2]
                    e1 = filtered[n+16]
                    e2 = filtered[n+18]
                    l1 = filtered[n+7]
                    l2 = filtered[n+8]
                    gap1, gap2, cov1, cov2, o, id = filter_1(gap1, gap2, cov1, cov2, o, f15,f1,f17,f3,f9,f24,f10,f25,f4,f19,f6,f21)

                elif filtered[n+12] == '-1':
                    f18 = int(filtered[n+18])
                    f2 = int(filtered[n+2])
                    frm = -1
                    if occur == 2:
                        s1 = filtered[n]
                        e2 = filtered[n+3]
                    e1 = filtered[n+16]
                    s2 = filtered[n+17]
                    l1 = filtered[n+7]
                    l2 = filtered[n+8]
                    gap1, gap2, cov1, cov2, o, id = filter_2(gap1, gap2, cov1, cov2, o, f15,f1,f18,f2,f9,f24,f10,f25,f4,f19,f6,f21)

            elif occur >= 1 and occur <= nuc_split and (cov1 >= merge_cov or cov2 >= merge_cov) and gap1 <= nuc_g and gap2 <= nuc_g: 
                outfile.write("\n" + str(s1) + "\t" + str(e1) + "\t" + str(s2) + "\t" + str(e2) + "\t" + str(o) + "\t" + str(o) + "\t" + str(id) + "\t" + str(l1) + "\t" + str(l2) + "\t" + str(cov1) + "\t" + str(cov2) + "\t1\t" + str(frm) + "\t" + qry + "\t" + sbj)
                occur = 1
                gap1 = 0
                gap2 = 0
                cov1 = 0
                cov2 = 0
                o = 0
            elif occur == 1 and f28 == '':
                outfile.write("\n" + filtered[n])
                for x in range(n+1,n+15):
                    outfile.write("\t" + filtered[x])
                occur = 1
                gap1 = 0
                gap2 = 0
                cov1 = 0
                cov2 = 0
                o = 0
            else: # occur > 1 but not written
                occur = 1
                gap1 = 0
                gap2 = 0
                cov1 = 0
                cov2 = 0
                o = 0

    subprocess.run("rm " + folder + "vRhyme_dereplication/vRhyme_nucmer_delta_filtered.tsv 2> /dev/null", shell=True)
    subprocess.run("rm " + folder + "vRhyme_dereplication/vRhyme_nucmer_delta_parsed.tsv 2> /dev/null", shell=True)
    filtered = None


@jit(nopython=True)
def sens_check_1(s1, l1, sens1, sens2, s2, l2, e1, e2, o1, o2):
    '''
    Sensitivity filtering for co-directional alignments
    '''
    check = False
    if (s1 < l1*sens1 or s2 < l2*sens1) and (l1-e1 < l1*sens1 or l2-e2 < l2*sens1) and (s2/o2 < sens2 or s1/o1 < sens2) and ((l2-e2)/o2 < sens2 or (l1-e1)/o1 < sens2):
        check = True
    return check

@jit(nopython=True)
def sens_check_2(s1, l1, sens1, sens2, s2, l2, e1, e2, o1, o2):
    '''
    Sensitivity filtering for anti-directional alignments
    '''
    check = False
    if (s1 < l1*sens1 or l2-s2 < l2*sens1) and (l1-e1 < l1*sens1 or e2 < l2*sens1) and ((l2-s2)/o2 < sens2 or s1/o1 < sens2) and (e2/o2 < sens2 or (l1-e1)/o1 < sens2):
        check = True
    return check

def nucmer_filter_3_stuff(cval, original_dict, derep_frac, sens1, sens2, derep_len, composite_dict, folder):
    '''
    Dereplication method for --method composite
    '''
    time.sleep(0.1)
    composited = []
    derep_frac = derep_frac*100
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    with open(folder + "vRhyme_dereplication/vRhyme_nucmer_delta_merged.tsv", "r") as merged:
        next(merged)
        aligned = []
        for line in merged:
            composite_seq = '' # generate composite seqs, 2 seqs at a time
            line = line.strip("\n").split("\t")

            qry = line[13]
            sbj = line[14]
            aligned.extend([qry,sbj])
            if qry not in composited and sbj not in composited and (float(line[9]) >= derep_frac or float(line[10]) >= derep_frac):
                rev_2 = False
                combo = [int(i) for i in line[0:6]+line[7:9]]
                s1 = combo[0]
                e1 = combo[1]
                s2 = combo[2]
                e2 = combo[3]
                o1 = combo[4]
                o2 = combo[5]
                l1 = combo[6] # original 7
                l2 = combo[7] # original 8
                if int(line[12]) == -1:
                    rev_2 = True

                if rev_2 == False:
                    check1 = sens_check_1(s1, l1, sens1, sens2, s2, l2, e1, e2, o1, o2)
                elif rev_2 == True:
                    check2 = sens_check_2(s1, l1, sens1, sens2, s2, l2, e1, e2, o1, o2)
                if rev_2 == False and check1 == True:

                    if s1 >= s2:
                        composite_seq += original_dict[qry][:e1]
                        composite_seq += original_dict[sbj][e2+1:]
                    else:
                        composite_seq += original_dict[sbj][:e2]
                        composite_seq += original_dict[qry][e1+1:]

                elif rev_2 == True and check2 == True: # if the seqs overlap too much in the unaligned region do not count

                    if l1-e1 >= l2-s2:
                        composite_seq += original_dict[qry][:e1]
                        composite_seq += ''.join([complement[base] for base in original_dict[sbj][:e2-1][::-1]])

                    elif l1-e1 < l2-s2:
                        composite_seq += ''.join([complement[base] for base in original_dict[sbj][s2:][::-1]])
                        composite_seq += original_dict[qry][s1+1:]

                if len(composite_seq) >= derep_len:
                    try:
                        del original_dict[qry]
                    except KeyError:
                        pass
                    try:
                        del original_dict[sbj]
                    except KeyError:
                        pass

                    original_dict.update({'vRhyme_composite_'+str(cval):composite_seq})

                    # if a composite sequence gets a hit then be sure to retain all seqs of that composite
                    header_list = []
                    composited.append(qry)
                    composited.append(sbj)

                    if qry.rsplit("_",1)[0] == 'vRhyme_composite':
                        header_list.extend(composite_dict[qry])
                    else:
                        header_list.append(qry)
                    if sbj.rsplit("_",1)[0] == 'vRhyme_composite':
                        header_list.extend(composite_dict[sbj])
                    else:
                        header_list.append(sbj)

                    header_list = list(set(header_list))
                    composite_dict.update({'vRhyme_composite_'+str(cval):header_list})
                    cval += 1

    subprocess.run("rm " + folder + "vRhyme_dereplication/vRhyme_nucmer_delta_merged.tsv 2> /dev/null", shell=True)
    return cval, composite_dict, original_dict


def nucmer_filter_4_stuff(original_dict, derep_frac, folder):
    '''
    Dereplication method for --method longest
    '''
    time.sleep(0.1)
    len_dict = {}
    longest_dict = {}
    bins = {}
    b = 1
    with open(folder + "vRhyme_dereplication/vRhyme_nucmer_delta_merged.tsv", "r") as coords:
        derep_frac = derep_frac*100
        next(coords)
        for line in coords:
            line = line.strip("\n").split("\t")
            if float(line[9]) >= derep_frac or float(line[10]) >= derep_frac:
                qry = line[13]
                sbj = line[14]
                len_dict.update({qry:int(line[7])})
                len_dict.update({sbj:int(line[8])})

                added = False
                for key in bins:
                    value = bins[key]
                    if qry in value or sbj in value:
                        value.update([qry,sbj])
                        bins[key] = value
                        added = True
                        break
                if added == False:
                    bins[b] = set([qry,sbj])
                    b += 1

        while True:
            names = list(bins.keys())
            l = len(names)
            alter = False
            for i in range(0,l-1):
                for j in range(i+1,l):
                    x = names[i]
                    y = names[j]
                    try:
                        x_values = bins[x]
                        y_values = bins[y]

                        if x_values & y_values != set():
                            x_values = x_values.union(y_values)
                            bins[x] = x_values

                            del bins[y]
                            alter = True

                    except Exception:
                        # y was deleted
                        pass

            if alter == False: # nothing was combined
                break

        names = None

        for names in bins.values():
            len_list = [len_dict[i] for i in names]
            z = list(zip(names,len_list))
            z.sort(key= lambda x: x[1], reverse=True)         
            name = z[0][0]
            longest_dict.update({name:names})

            for i,_ in z[1:]:
                del original_dict[i]

    len_dict = None
    bins = None
    return original_dict, longest_dict

#############################################################################################

def derep_stuff_composite(in_fasta, threads, interest_list, folder, base, derep_len, nuc_c, nuc_b, nuc_g, nuc_split, sens1, sens2, derep_id, derep_frac, sketch_k, sketch_s, iterations):
    '''
    Main wrapper and iterator for --method composite
    '''
    
    merge_cov = 1-derep_frac
    # start with splitting scaffolds to files
    if interest_list != None:
        original_seqs, start_total = split_singles_stuff_interest(in_fasta, base, derep_len, folder)
    else:
        original_seqs, start_total = split_singles_stuff(in_fasta, base, derep_len, folder)

    check_duplicates = len(list(set(original_seqs.keys())))
    if check_duplicates != start_total:
        with open(folder + "vRhyme_dereplication/" + base + '.vRhyme-unique.fa', 'w') as loop_fa:
            for key in original_seqs.keys():
                loop_fa.write(">" + key.replace("~#~"," ") + "\n" + original_seqs[key] + "\n")
        subprocess.run("rm -R " + folder + "vRhyme_dereplication/split_fasta_files 2> /dev/null", shell=True)
        return 'error1', check_duplicates, start_total

    iteration = 0

    composite_dict = {}
    clusters_100_dict = {}
    cval = 1
    while True:
        if iteration == 0:
            keep_start = start_total
        mash_stuff(base, threads, start_total, folder, sketch_k, sketch_s, iterations)
        aligns = mash_parse_stuff(base, folder)
        subprocess.run("rm " + folder + "vRhyme_dereplication/mash_dist_" + base + ".tsv 2> /dev/null", shell=True)
        subprocess.run("rm " + folder + "vRhyme_dereplication/mash_dist_" + base + ".temp.tsv 2> /dev/null", shell=True)

        start_total = len(original_seqs.keys())
        if aligns:
            nucmer_stuff(base, threads, nuc_c, nuc_b, nuc_g, folder)
            nucmer_parse_stuff(threads, derep_id, derep_len, folder)
            if iteration == 0:
                clusters_100_dict, original_seqs, coords_check = nucmer_filter_1_stuff(original_seqs, folder)
            
            if coords_check:
                nucmer_filter_2_stuff(iteration, nuc_split, merge_cov, nuc_g, folder)
                cval, composite_dict, original_seqs = nucmer_filter_3_stuff(cval, original_seqs, derep_frac, sens1, sens2, derep_len, composite_dict, folder)
        if aligns == False and iteration == 0:
            subprocess.run("rm -R " + folder + "vRhyme_dereplication/split_fasta_files 2> /dev/null", shell=True)
            subprocess.run("rm " + folder + "vRhyme_dereplication/vRhyme_nucmer_delta_parsed.tsv 2> /dev/null", shell=True)
            subprocess.run("rm " + folder + "vRhyme_dereplication/vRhyme_nucmer_delta_filtered.tsv 2> /dev/null", shell=True)
            return 'noalign', start_total, None

        subprocess.run("rm -R " + folder + "vRhyme_dereplication/split_fasta_files 2> /dev/null", shell=True)
        time.sleep(0.1)

        if len(original_seqs.keys()) == start_total:
            break

        subprocess.run("mkdir " + folder + "vRhyme_dereplication/split_fasta_files", shell=True)
        outfolder = folder + "vRhyme_dereplication/split_fasta_files/"

        n = 1
        for key,val in original_seqs.items():
            with open(outfolder + str(n) + ".fa", 'w') as outfasta:
                outfasta.write(">" + key + "\n" + val + "\n")
                n += 1

        iteration += 1

    total_derep = len(original_seqs.keys())

    if keep_start != total_derep:
        with open(folder + 'vRhyme_dereplication/vRhyme_derep_composite_' + base + '.fa', 'w') as loop_fa:
            for key,val in original_seqs.items():
                loop_fa.write(">" + key.replace("~#~"," ") + "\n" + val + "\n")
        in_fasta = folder + 'vRhyme_dereplication/vRhyme_derep_composite_' + base + '.fa'

        if len(clusters_100_dict.keys()) > 0:
            with open(folder + 'vRhyme_dereplication/vRhyme_derep_indentical-seqs_' + base + '.tsv', 'w') as ident_seq:
                ident_seq.write("parent\tmembers ->")
                for item,val in clusters_100_dict.items():
                    val = [v.replace("~#~"," ") for v in val]
                    ident_seq.write("\n" + item.replace("~#~"," ") + "\t" + "\t".join(val))

        if len(composite_dict.keys()) > 0:
            set_seqs = set(list(original_seqs.keys()))
            with open(folder + 'vRhyme_dereplication/vRhyme_derep_composited-seqs_' + base + '.tsv', 'w') as comp_seq, open(folder + 'vRhyme_dereplication/vRhyme_derep_composited-list_' + base + '.txt', 'w') as comp_list:
                comp_seq.write("parent\tmembers ->")
                for item,val in composite_dict.items():
                    if item in set_seqs:
                        val = [v.replace("~#~"," ") for v in val]
                        comp_seq.write("\n" + item.replace("~#~"," ") + "\t" + "\t".join(val))
                        comp_list.write("\n".join(val) + "\n")

        return 'finished', keep_start, total_derep
    else:
        return 'noalign', None, start_total

def derep_stuff_longest(in_fasta, threads, interest_list, folder, base, derep_len, nuc_c, nuc_b, nuc_g, nuc_split, derep_id, derep_frac, sketch_k, sketch_s, iterations):
    '''
    Main wrapper and iterator for --method longest
    '''
    
    merge_cov = derep_frac
    # start with splitting scaffolds to files
    if interest_list != None:
        original_seqs, start_total = split_singles_stuff_interest(in_fasta, base, derep_len, folder)
    else:
        original_seqs, start_total = split_singles_stuff(in_fasta, base, derep_len, folder)

    check_duplicates = len(list(set(original_seqs.keys())))
    if check_duplicates != start_total:
        with open(folder + "vRhyme_dereplication/" + base + '.vRhyme-unique.fa', 'w') as loop_fa:
            for key,val in original_seqs.items():
                loop_fa.write(">" + key.replace("~#~"," ") + "\n" + val + "\n")
        subprocess.run("rm -R " + folder + "vRhyme_dereplication/split_fasta_files 2> /dev/null", shell=True)
        return 'error1', check_duplicates, start_total

    mash_stuff(base, threads, start_total, folder, sketch_k, sketch_s, iterations)
    aligns = mash_parse_stuff(base, folder)
    subprocess.run("rm " + folder + "vRhyme_dereplication/mash_dist_" + base + ".tsv 2> /dev/null", shell=True)
    subprocess.run("rm " + folder + "vRhyme_dereplication/mash_dist_" + base + ".temp.tsv 2> /dev/null", shell=True)

    if aligns == True:
        nucmer_stuff(base, threads, nuc_c, nuc_b, nuc_g, folder)
        nucmer_parse_stuff(threads, derep_id, derep_len, folder)
        clusters_100_dict, original_seqs, _ = nucmer_filter_1_stuff(original_seqs, folder)
        nucmer_filter_2_stuff(0, nuc_split, merge_cov, nuc_g, folder)
        original_seqs, longest_dict = nucmer_filter_4_stuff(original_seqs, derep_frac, folder)

    final_seqs = len(original_seqs.keys())

    if start_total != final_seqs:
        with open(folder + 'vRhyme_dereplication/vRhyme_derep_longest_' + base + '.fa', 'w') as loop_fa:
            for key,val in original_seqs.items():
                loop_fa.write(">" + key.replace("~#~"," ") + "\n" + val + "\n")

        if longest_dict:
            with open(folder + 'vRhyme_dereplication/vRhyme_derep_overlap-seqs_' + base + '.tsv', 'w') as ident_seq:
                ident_seq.write("parent\tmembers ->")
                for item,val in longest_dict.items():
                    val = [v.replace("~#~"," ") for v in val]
                    ident_seq.write("\n" + item.replace("~#~"," ") + "\t" + "\t".join(val))

        if clusters_100_dict:
            with open(folder + 'vRhyme_dereplication/vRhyme_derep_indentical-seqs_' + base + '.tsv', 'w') as ident_seq:
                ident_seq.write("parent\tmembers ->")
                for item,val in clusters_100_dict.items():
                    val = [v.replace("~#~"," ") for v in val]
                    ident_seq.write("\n" + item.replace("~#~"," ") + "\t" + "\t".join(val))

        subprocess.run("rm -R " + folder + "vRhyme_dereplication/split_fasta_files 2> /dev/null", shell=True)
        subprocess.run("rm -R " + folder + "vRhyme_dereplication/vRhyme_nucmer_delta_merged.tsv 2> /dev/null", shell=True)
        return 'finished', start_total, final_seqs


    else:
        subprocess.run("rm -R " + folder + "vRhyme_dereplication/split_fasta_files 2> /dev/null", shell=True)
        subprocess.run("rm -R " + folder + "vRhyme_dereplication/vRhyme_nucmer_delta_merged.tsv 2> /dev/null", shell=True)
        return 'noalign', start_total, None


#
#
#
