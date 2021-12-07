#! /usr/bin/env python3
# vRhyme
# Author: Kristopher Kieft
# University of Wisconsin-Madison

import subprocess

def linclust(folder, threads):
    '''
    Mmseqs2 linclust subprocess
    '''
    subfolder = f'vRhyme_linclust_clustering/'
    faa = f'{folder}{subfolder}linclust_proteins_mapped.faa'
    db = f'{folder}{subfolder}linclust_proteins_mapped.db'
    result = f'{folder}{subfolder}vRhyme_linclust_results.raw'
    t = f'{folder}{subfolder}temp_mmseqs_dir'

    subprocess.run("mmseqs createdb --dbtype 1 -v 1 " + faa + " " + db + " > /dev/null 2> /dev/null", shell=True)
    subprocess.run("mmseqs linclust " + db + " " + result + " " + t + " -v 1 --min-seq-id 0.5 -c 0.8 -e 0.01 --min-aln-len 50 --threads " + threads + " --cluster-mode 0 --seq-id-mode 0 --alignment-mode 3 --cov-mode 5 --kmer-per-seq 75 > /dev/null 2> /dev/null", shell=True) 

    out = f'{folder}{subfolder}vRhyme_linclust_results.tsv'
    subprocess.run("mmseqs createtsv " + db + " " + db + " " + result + " " + out + " -v 1 --threads " + threads + " > /dev/null 2> /dev/null", shell=True)


def linclust_threader(folder, faa):
    '''
    Mmseqs2 linclust subprocess
    '''
    result = faa.rsplit('.',1)[0] + '_linclust.result'
    t = folder + 'temp_mmseqs_dir'
    db = faa + ".mmseqs2.db"
    subprocess.run("mmseqs createdb --dbtype 1 -v 1 " + faa + " " + db + " > /dev/null 2> /dev/null", shell=True)
    s = subprocess.Popen("mmseqs linclust " + db + " " + result + " " + t + " -v 1 --min-seq-id 0.5 -c 0.8 -e 0.01 --min-aln-len 50 --threads 1 --cluster-mode 0 --seq-id-mode 0 --alignment-mode 3 --cov-mode 5 --kmer-per-seq 55 > /dev/null 2> /dev/null", shell=True) 
    return s

def linclust_createtsv(db,result):
    '''
    Extract linclust results
    '''
    out = result.rsplit('.',1)[0] + '.tsv'
    s = subprocess.Popen("mmseqs createtsv " + db + " " + db + " " + result + " " + out + " -v 1 --threads 1  > /dev/null 2> /dev/null", shell=True)
    return s

def parse_linclust(folder):
    '''
    Parse linclust results
    '''
    subfolder = f'vRhyme_linclust_clustering/'
    results = f'{folder}{subfolder}vRhyme_linclust_results.tsv'

    master = {} # {genome: [rep clusters]}
    cluster = 0
    try:
        with open(results, 'r') as in_table:
            in_table = in_table.read().split("\n")
        if in_table[-1] == '':
            in_table = in_table[:-1]
        in_table.append('placeholder_1\tplaceholder_2') # lets the last line of in_table enter master
    except Exception:
        return master
    
    holder = set()
    for item in in_table:
        item = item.split('\t')
        if item[0] == item[1]:
            if len(holder) > 1:
                cluster += 1
                for seq in holder:
                    master.setdefault(int(seq),[]).append(cluster) 
            holder = set()
            qry = item[0].split('_')[0]
            holder.add(qry)
        else:
            qry = item[0].split('_')[0]
            sbj = item[1].split('_')[0]
            holder.update([qry,sbj])
            
    return master

#
#
#
