#! /usr/bin/env python3
# vRhyme
# Author: Kristopher Kieft
# University of Wisconsin-Madison


import os
import subprocess


def bwa_build_stuff(in_fasta, spaces, folder):
    '''
    Build BWA index database file
    '''
    try:
        temp = in_fasta.rsplit("/",1)[1]
        base = temp.rsplit(".",1)[0]
    except Exception:
        base = in_fasta.rsplit(".",1)[0]

    if spaces == True: # Bowtie2 doesn't like spaces

        subprocess.run("cat " + in_fasta + " | sed 's/ /~(#)~/g' > " + folder + base + ".no-spaces.fasta", shell=True)

        if not os.path.exists(folder + base + ".no-spaces.bwa.index.sa"):
            subprocess.run("bwa index -p " + folder + base + ".no-spaces.bwa.index " + folder + base + ".no-spaces.fasta > /dev/null 2> /dev/null", shell=True)

        return folder + base + ".no-spaces.bwa.index"

    else:
        if not os.path.exists(folder + base + ".bwa.index.sa"):
            subprocess.run("bwa index -p " + folder + base + ".bwa.index " + in_fasta + " > /dev/null 2> /dev/null", shell=True)

        return folder + base + ".bwa.index"



def bwa_paired_stuff(threads, forward, reverse, folder, build, keep_sam, spaces):
    '''
    Run BWA on paired reads files
    '''
    try:
        temp = forward.rsplit("/",1)[1]
        base = temp.rsplit(".fastq",1)[0]
    except Exception:
        base = forward.rsplit(".fastq",1)[0]

    subprocess.run("bwa mem -t " + threads + " -o " + str(folder) + "vRhyme_sam_files/" + str(base)+".sam " + build + " " + forward + " " + reverse + " > /dev/null 2> /dev/null", shell=True)

    sam = str(folder) + "vRhyme_sam_files/" + str(base)+".sam"

    if spaces == True:
        subprocess.run("cat " + sam + " | sed 's/~(#)~/ /g' > " + folder + "vRhyme_sam_files/" + base + ".replaced.sam", shell=True)
        subprocess.run("rm " + sam, shell=True)
        subprocess.run("mv " + folder + "vRhyme_sam_files/" + base + ".replaced.sam " + sam, shell=True)

    output = folder + "vRhyme_bam_files/" + base + ".bam"
    subprocess.run("samtools view -@ " + threads + " -h -b " + sam + " > " + output + " 2> /dev/null", shell=True)
    if keep_sam == False:
        subprocess.run("rm " + sam + " 2> /dev/null", shell=True)

    return folder + "vRhyme_bam_files/" + base + ".bam"

def bwa_single_stuff(threads, single, folder, build, keep_sam, spaces):
    '''
    Run BWA on single end reads files
    '''
    try:
        temp = single.rsplit("/",1)[1]
        base = temp.rsplit(".",1)[0]
    except Exception:
        base = single.rsplit(".",1)[0]

    subprocess.run("bwa mem -t " + threads + " -o " + str(folder) + "vRhyme_sam_files/" + str(base)+".sam " + build + " " + single + " > /dev/null 2> /dev/null", shell=True) #> /dev/null 2> /dev/null

    sam = str(folder) + "vRhyme_sam_files/" + str(base)+".sam"
    if spaces == True:
        subprocess.run("cat " + sam + " | sed 's/~(#)~/ /g' > " + folder + "vRhyme_sam_files/" + base + ".replaced.sam", shell=True)
        subprocess.run("rm " + sam, shell=True)
        subprocess.run("mv " + folder + "vRhyme_sam_files/" + base + ".replaced.sam " + sam, shell=True)

    output = folder + "vRhyme_bam_files/" + base + ".bam"
    subprocess.run("samtools view -@ " + threads + " -h -b " + sam + " > " + output + " 2> /dev/null", shell=True)
    if keep_sam == False:
        subprocess.run("rm " + sam + " 2> /dev/null", shell=True)

    return folder + "vRhyme_bam_files/" + base + ".bam"


#
#
#
