#! /usr/bin/env python3
# vRhyme
# Author: Kristopher Kieft
# University of Wisconsin-Madison

import subprocess
import os


def prodigal_stuff(in_fasta, folder, base, spaces, comp):
    '''
    Run Prodigal on the split runs files (parallel)
    '''
    try:
        temp = in_fasta.rsplit("/",1)[1]
        base = temp.rsplit(".",1)[0]
    except Exception:
        base = in_fasta.rsplit(".",1)[0]

    fna_list = os.listdir(folder + "vRhyme_split_runs/")
    s_list = []
    for fna in fna_list:
        if fna.rsplit(".",1)[1] == "fna":
            n = str(fna.rsplit(".",1)[0])
            s = subprocess.Popen("prodigal -m -p meta -i " + folder + "vRhyme_split_runs/" + fna + " -a " + folder + "vRhyme_split_runs/" + n + ".faa -d " + folder + "vRhyme_split_runs/" + n + ".ffn -q > /dev/null", shell=True)
            s_list.append(s)
    for item in s_list:
        item.wait()

    if comp == False:
        s1 = subprocess.Popen("cat " + folder + "vRhyme_split_runs/*.faa > " + folder + base + ".prodigal.faa", shell=True)
        s2 = subprocess.Popen("cat " + folder + "vRhyme_split_runs/*.ffn > " + folder + base + ".prodigal.ffn", shell=True)
        s1.wait()
        s2.wait()
    else:
        s1 = subprocess.Popen("cat " + folder + "vRhyme_split_runs/*.faa > " + folder + base + ".prodigal_composites.faa", shell=True)
        s2 = subprocess.Popen("cat " + folder + "vRhyme_split_runs/*.ffn > " + folder + base + ".prodigal_composites.ffn", shell=True)
        s1.wait()
        s2.wait()

    if spaces == True:
        fasta_files = os.listdir(folder + "vRhyme_split_runs/")
        s_list = []
        for fna in fasta_files:
            f = fna.rsplit(".",1)
            if f[1] == 'fna':
                fna = folder + "vRhyme_split_runs/" + fna
                b = folder + "vRhyme_split_runs/" + f[0] + ".replaced.fa"
                s = subprocess.Popen("cat " + fna + " | sed 's/~(#)~/ /g' > " + b, shell=True)
                s_list.append(s)
        for s in s_list:
            s.wait()
        subprocess.run("rm " + folder + "vRhyme_split_runs/*.fna", shell=True)

        if comp == False:
            s1 = subprocess.Popen("cat " + folder + base + ".prodigal.faa | sed 's/~(#)~/ /g' > " + folder + base + ".prodigal.replaced.faa", shell=True)
            s2 = subprocess.Popen("cat " + folder + base + ".prodigal.ffn | sed 's/~(#)~/ /g' > " + folder + base + ".prodigal.replaced.ffn", shell=True)
            s1.wait()
            s2.wait()
            subprocess.run("mv " + folder + base + ".prodigal.replaced.faa " + folder + base + ".prodigal.faa", shell=True)
            subprocess.run("mv " + folder + base + ".prodigal.replaced.ffn " + folder + base + ".prodigal.ffn", shell=True)
            return folder + base + ".prodigal.faa", folder + base + ".prodigal.ffn"
        else:
            s1 = subprocess.Popen("cat " + folder + base + ".prodigal_composites.faa | sed 's/~(#)~/ /g' > " + folder + base + ".prodigal.replaced.faa", shell=True)
            s2 = subprocess.Popen("cat " + folder + base + ".prodigal_composites.ffn | sed 's/~(#)~/ /g' > " + folder + base + ".prodigal.replaced.ffn", shell=True)
            s1.wait()
            s2.wait()
            subprocess.run("mv " + folder + base + ".prodigal.replaced.faa " + folder + base + ".prodigal_composites.faa", shell=True)
            subprocess.run("mv " + folder + base + ".prodigal.replaced.ffn " + folder + base + ".prodigal_composites.ffn", shell=True)
            return folder + base + ".prodigal_composites.faa", folder + base + ".prodigal_composites.ffn"


    return folder + base + ".prodigal.faa", folder + base + ".prodigal.ffn"


#
#
#
