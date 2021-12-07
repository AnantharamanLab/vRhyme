#! /usr/bin/env python3
# vRhyme
# Author: Kristopher Kieft
# University of Wisconsin-Madison


import argparse
import os
import sys
import uuid
import subprocess
from collections import Counter


def linclust(folder, faa, threads, spaces):
    db = f'{folder}linclust_database.db'
    result = f'{folder}linclust_results.raw'
    t = f'{folder}temp_mmseqs_dir/'

    if spaces == True:
        no_spaces = f'{folder}linclust_no-spaces.faa'
        subprocess.run("cat " + faa + " | sed 's/ # /#~#/g' | sed 's/ /!~!/g' > " + no_spaces + " 2> /dev/null", shell=True)
        subprocess.run("mmseqs createdb --dbtype 1 -v 1 " + no_spaces + " " + db + " > /dev/null 2> /dev/null", shell=True)
    else:
        subprocess.run("mmseqs createdb --dbtype 1 -v 1 " + faa + " " + db + " > /dev/null 2> /dev/null", shell=True)
    subprocess.run("mmseqs linclust " + db + " " + result + " " + t + " -v 1 --min-seq-id 0.5 -c 0.8 -e 0.01 --min-aln-len 50 --threads " + threads + " --cluster-mode 0 --seq-id-mode 0 --alignment-mode 3 --cov-mode 5 --kmer-per-seq 75 > /dev/null 2> /dev/null", shell=True) # cluster-mode 5

    out = f'{folder}linclust_results.tsv'
    subprocess.run("mmseqs createtsv " + db + " " + db + " " + result + " " + out + " -v 1 --threads " + threads + " > /dev/null 2> /dev/null", shell=True)

def parse_bins(infile, zero, header, rev, comma):
    bins = {}
    delimiter = "\t"
    if comma == True:
        delimiter = ","
    with open(infile, 'r') as binsfile:
        if header == False:
            next(binsfile)
        if rev == False:
            for line in binsfile:
                try:
                    line = line.rstrip("\n").split(delimiter)
                    bins.setdefault(int(line[1]), []).append(line[0])
                except Exception:
                    pass
        else:
            for line in binsfile:
                try:
                    line = line.rstrip("\n").split(delimiter)
                    bins.setdefault(int(line[0]), []).append(line[1])
                except Exception:
                    pass
    if zero == False:
        bins.pop(0,None)
    
    return bins


def parse_linclust(folder, spaces):

    result = f'{folder}linclust_results.tsv'

    master = {} # {genome: [rep clusters]}
    cluster = 0
    try:
        with open(result, 'r') as in_table:
            in_table = in_table.read().split("\n")
        if in_table[-1] == '':
            in_table = in_table[:-1]
        in_table.append('placeholder_1\tplaceholder_2') # lets the last line of in_table enter master
    except Exception as e:
        print(str(e))
        return master
    
    holder = set()
    for item in in_table:
        item = item.split('\t')
        if item[0] == item[1]:
            if len(holder) > 1:
                cluster += 1
                for seq in holder:
                    if spaces == True:
                        temp = seq.split('#~#',1)[0]
                        seq = temp.replace('!~!', ' ')
                    master.setdefault(seq,[]).append(cluster) 
            holder = set()
            qry = item[0].rsplit('_',1)[0]
            holder.add(qry)
        else:
            qry = item[0].rsplit('_',1)[0]
            sbj = item[1].rsplit('_',1)[0]
            holder.update([qry,sbj])

    return master

def prots_per_seq(faa):
    proteins_count = {}
    spaces = False
    with open(faa, 'r') as proteins:
        seq = ''
        for line in proteins:
            if line.startswith(">"):
                if seq != '':
                    base = header.rsplit('_',1)[0]
                    if ' ' in base:
                        spaces = True
                    try:
                        proteins_count[base] += 1
                    except Exception:
                        proteins_count[base] = 1
                    seq = ''
                header = line[1:]
                try:
                    header = header.split(" # ",1)[0]
                except Exception:
                    header = header.strip("\n")
            else:
                seq += line.strip("\n")

        base = header.rsplit('_',1)[0]
        if ' ' in base:
            spaces = True
        try:
            proteins_count[base] += 1
        except Exception:
            proteins_count[base] = 1

    return proteins_count, spaces

def scoring(outfile, seqs, original, bins, proteins_count, linclust_master, red_cutoff):
    with open(outfile,'w') as summary:
        summary.write('bin\tmembers\tproteins\tredundancy\n')
        redundancy_master = 0
        proteins_master = 0
        total_bins = 0
        for b in bins:
            proteins = 0
            clusters = []
            members = bins[b]
            len_mem = len(members)
            if len_mem > 1:
                for mem in members:
                    proteins += proteins_count.get(mem,0)
                    clusters.extend(linclust_master.get(mem,[]))

                dup_dict = Counter(clusters)
                dup_list = [val for val in dup_dict.values() if val > 1]
                duplicated_members = sum(dup_list)-len(dup_list) # total redundant proteins - representative protein in each duplicated cluster

                try:
                    red_check = duplicated_members/proteins 
                except Exception:
                    red_check = 0 # zero proteins
                
                if red_check < red_cutoff or proteins <= 5:

                    summary.write(f'{b}\t{len(members)}\t{proteins}\t{duplicated_members}\n')
                    #
                    redundancy_master += duplicated_members
                    proteins_master += proteins
                    total_bins += 1
                else:
                    seqs -= len_mem
            else:
                seqs -= 1

        if seqs > 0:
            f = seqs/original
            try:
                r = (2*(redundancy_master/proteins_master))**0.5
            except Exception:
                r = 0 # no proteins

            b = (total_bins/seqs)**2
            s = round(f - b - (3*r),4)
        else:
            s = 0

        summary.write(f'{s}\t{seqs}\t{proteins_master}\t{redundancy_master}\n')

    

if __name__ == '__main__':
    os.environ["COLUMNS"] = '100'

    descript = """
    Calculate vRhyme equivalent scores per bin.
    The outfile will contain n rows per bin;
    the final row is a summary of all bins 
    with the total score in the 'bin' column."""
    vRhyme = argparse.ArgumentParser(description=descript, formatter_class=argparse.RawTextHelpFormatter)
    vRhyme.add_argument('-i', metavar='', type=str, nargs=1, required=True, help="input tab-separated bin membership file (e.g., scaffold \\t bin)")
    vRhyme.add_argument('-p', metavar='', type=str, nargs=1, required=True, help="input proteins fasta file (input complete protein file, can input binned only)")
    vRhyme.add_argument('-o', metavar='', type=str, nargs=1, required=True, help='output tab-separated bin scores file')
    vRhyme.add_argument('-n', metavar='', type=str, nargs=1, required=True, help='number of original scaffolds before binning (binned + unbinned)')
    vRhyme.add_argument('-t', metavar='', type=str, nargs=1, default=['5'], help='number of threads to use for mmseqs2 linclust protein clustering [5]')
    vRhyme.add_argument('-r', metavar='', type=str, nargs=1, default=['1.00'], help="percent protein redundancy threshold per bin; equivalent to vRhyme --red [1.00]")
    vRhyme.add_argument('--comma', action='store_true', help="set to indicate that the -i input table is in comma-separated format [off]")
    vRhyme.add_argument('--nohead', action='store_true', help="set to indicate that the -i input table has no header [off]")
    vRhyme.add_argument('--rev', action='store_true', help="set to indicate that the -i input table columns are reversed (bin  scaffold) [off]")
    vRhyme.add_argument('--zero', action='store_true', help="set to consider bin label '0' as a vaid bin [off]")
    #
    args = vRhyme.parse_args()
    #
    faa = str(args.p[0])
    members = str(args.i[0])
    outfile = str(args.o[0])
    original = int(args.n[0])
    threads = str(args.t[0])
    red_cutoff = float(args.r[0])
    if red_cutoff > 1:
        sys.stderr.write(f"\nError: input -r as a decimal value. Example, input '--red 0.25' for 25%. Exiting.\n")
        exit(1)
    #
    if os.path.exists(str(args.o[0])):
        sys.stderr.write("\nError: output file (-o) already exists. Exiting.\n\n")
        exit()
    #
    bins = parse_bins(members, args.zero, args.nohead, args.rev, args.comma)
    if bins == {}:
        sys.stderr.write("\nNo bins were identified. Verify that the input file (-i) format is tab-separated. Exiting.\n\n")
        exit()
    seqs = 0
    for val in bins.values():
        seqs += len(val)
    if seqs > original:
        sys.stderr.write("\nError: the number of identified binned sequences (rows in -i) exceeds original binned+unbinned sequences (-n). Exiting.\n\n")
        exit()
    u = str(uuid.uuid1()).split("-")[0]
    folder = f'vRhyme_scores_temp_dir_{u}/'
    subprocess.run(f'mkdir {folder}', shell=True)
    #
    proteins_count, spaces = prots_per_seq(faa)
    linclust(folder, faa, threads, spaces)
    linclust_master = parse_linclust(folder, spaces)
    scoring(outfile, seqs, original, bins, proteins_count, linclust_master, red_cutoff)

    subprocess.run(f'rm -R {folder}', shell=True)
#
#
#
