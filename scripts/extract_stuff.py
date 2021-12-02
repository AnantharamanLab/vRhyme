#! /usr/bin/env python3
# vRhyme
# Author: Kristopher Kieft
# University of Wisconsin-Madison

import subprocess
import math
import copy


def kmer_split_stuff(in_fasta, folder, threads, pair_genomes, keep):
    '''
    Read fasta file and split into separate files for parallel running
    '''
    subprocess.run("rm -R " + folder + "vRhyme_split_runs/ 2> /dev/null", shell=True)
    subprocess.run("mkdir " + folder + "vRhyme_split_runs", shell=True)
    outfolder = folder + "vRhyme_split_runs/"
    total = len(pair_genomes)
    number = math.ceil(int(total)/int(threads))

    n = 1
    i = 0
    with open(in_fasta, 'r') as fasta:
        seq = ''
        for line in fasta:
            if line.startswith(">"):
                if i >= number:
                    n += 1
                    i = 0
                try: # not in keep
                    header = keep[header][1]
                    if seq != '' and header in pair_genomes:
                        i += 1
                        with open(outfolder + str(n) + ".fa", 'a') as outfasta:
                            outfasta.write(">" + str(header) + "\n" + seq + "\n")
                except Exception:
                    pass
                seq = ''
                header = line[1:].strip("\n") # no ">", no "\n"
            else:
                seq += line.strip("\n")
        try:
            header = keep[header][1]
            if header in pair_genomes:
                with open(outfolder + str(n) + ".fa", 'a') as outfasta:
                    outfasta.write(">" + str(header) + "\n" + seq + "\n")
        except Exception:
            pass


def split_stuff(in_fasta, folder, threads, pair_genomes, keep):
    '''
    Read fasta file and split into separate files for parallel running,
    replace spaces with special characters
    '''
    subprocess.run("rm -R " + folder + "vRhyme_split_runs/ 2> /dev/null", shell=True)
    subprocess.run("mkdir " + folder + "vRhyme_split_runs", shell=True)
    outfolder = folder + "vRhyme_split_runs/"
    total = len(pair_genomes)
    number = math.ceil(int(total)/int(threads))

    n = 1
    i = 0
    with open(in_fasta, 'r') as fasta:
        seq = ''
        for line in fasta:
            if line.startswith(">"):
                if i >= number:
                    n += 1
                    i = 0
                try:
                    if seq != '' and keep[header][1] in pair_genomes:
                        i += 1
                        with open(outfolder + str(n) + ".fna", 'a') as outfasta:
                            outfasta.write(">" + header.replace(" ", "~(#)~") + "\n" + seq + "\n")
                except Exception:
                    pass
                seq = ''
                header = line[1:].strip("\n") # no ">", no "\n"
            else:
                seq += line.strip("\n")
        try:
            if keep[header][1] in pair_genomes:
                with open(outfolder + str(n) + ".fna", 'a') as outfasta:
                    outfasta.write(">" + header.replace(" ", "~(#)~") + "\n" + seq + "\n")
        except Exception:
            pass


def split_stuff_composite(in_fasta, folder, threads, pair_genomes, mapper):
    '''
    Read fasta file and split into separate files for parallel running,
    handle composited sequences
    '''
    subprocess.run("rm -R " + folder + "vRhyme_split_runs 2> /dev/null", shell=True)
    subprocess.run("mkdir " + folder + "vRhyme_split_runs", shell=True)
    outfolder = folder + "vRhyme_split_runs/"
    composite_list = []
    for name in list(pair_genomes):
        map_name = mapper[name]
        if 'vRhyme_composite_' in map_name:
            composite_list.append(map_name)
    total = len(composite_list)
    number = math.ceil(int(total)/int(threads))
    composite_list = set(composite_list)

    n = 1
    i = 0
    with open(in_fasta, 'r') as fasta:
        seq = ''
        for line in fasta:
            if line.startswith(">"):
                if i >= number:
                    n += 1
                    i = 0
                if seq != '' and header in composite_list:
                    i += 1
                    with open(outfolder + str(n) + ".fna", 'a') as outfasta:
                        outfasta.write(">" + header.replace(" ", "~(#)~") + "\n" + seq + "\n")
                seq = ''
                header = line[1:].strip("\n") # no ">", no "\n"
            else:
                seq += line.strip("\n")
        if header in composite_list:
            with open(outfolder + str(n) + ".fna", 'a') as outfasta:
                outfasta.write(">" + header.replace(" ", "~(#)~") + "\n" + seq + "\n")

    composite_list = None


def split_genes_stuff(in_genes, folder, pair_genomes, keep, threads):
    '''
    Read fasta file and split into separate files for parallel running,
    functions for genes in Prodigal format
    '''
    subprocess.run("rm -R " + folder + "vRhyme_split_runs/ 2> /dev/null", shell=True)
    subprocess.run("mkdir " + folder + "vRhyme_split_runs", shell=True)
    outfolder = folder + "vRhyme_split_runs/"
    total = len(pair_genomes)
    number = math.ceil(int(total)/int(threads))

    n = 1
    i = 0
    with open(in_genes, 'r') as fasta:
        seq = ''
        for line in fasta:
            if line.startswith(">"):
                if i >= number:
                    n += 1
                    i = 0
                header = line[1:].strip("\n") # no ">", no "\n"
                try:
                    header = header.split(" # ",1)[0]
                except Exception:
                    pass
                split_head = header.rsplit("_",1)[0]
                length = len(seq)%3
                if length == 1:
                    seq += '**'
                elif length == 2:
                    seq += '*' # separate the sequences so 3-mers aren't connected
                try:
                    if split_head != save:
                        save_map = keep[save][1]
                        if save_map in pair_genomes:
                            i += 1
                            with open(outfolder + str(n) + ".ffn", 'a') as outfasta:
                                outfasta.write(">" + str(save_map) + "\n" + seq + "\n")
                        seq = ''
                except Exception:
                    seq = ''
                    pass
                save = copy.copy(split_head)

            else:
                seq += line.strip("\n").strip("*") # no ">", no "\n"

        try:
            save_map = keep[save][1]
            if save_map in pair_genomes:
                with open(outfolder + str(n) + ".ffn", 'a') as outfasta:
                    outfasta.write(">" + str(save_map) + "\n" + seq + "\n")
        except Exception:
            pass

    seq = None
    line = None


def replace_mapper(in_prots, folder, keep):
    '''
    Read fasta file and split into separate files for parallel running,
    functions for proteins in Prodigal format
    '''
    subfolder = f'vRhyme_linclust_clustering/'
    subprocess.run(f'mkdir {folder}{subfolder}', shell=True)
    with open(in_prots, 'r') as fasta, open(f'{folder}{subfolder}linclust_proteins_mapped.faa', 'w') as outfasta:
        seq = ''
        prots_dict = {}
        for line in fasta:
            if line.startswith(">"):
                if seq != '' and name != None:
                    outfasta.write(f"{header}\n{seq}\n")
                    seq = ''
                
                header = line[1:]
                try:
                    header = header.split(" # ",1)[0]
                except Exception:
                    header = header.strip("\n")
                
                headers = header.rsplit("_",1)
                try:
                    name = keep[headers[0]][1]
                    try:
                        prots_dict[name] += 1
                    except Exception:
                        prots_dict[name] = 1
                    num = headers[1]
                    header = f'>{name}_{num}'
                except Exception:
                    name = None
                    header = None
                
            else:
                seq += line.strip("\n")
        if name != None:
            outfasta.write(f"{header}\n{seq}\n")

    return prots_dict


def score_extract(folder, in_prots, uniques, keep):
    '''
    Read fasta file and count proteins per scaffold for scoring,
    also writes out proteins per unique bin for Mmseqs2
    '''
    subprocess.run(f'mkdir {folder}vRhyme_unique_bins_proteins', shell=True)

    db = {}
    prots_dict = {}
    with open(in_prots, 'r') as fasta:
        seq = ''
        for line in fasta:
            if line.startswith(">"):
                if seq != '':
                    db.setdefault(genome,[]).append((f'{genome}_{num}',seq))
                seq = ''
                header = line[1:]
                try:
                    header = header.split(" # ",1)[0]
                except Exception:
                    header = header.strip("\n")
                
                headers = header.rsplit("_",1)
                genome = keep[headers[0]][1]
                num = headers[1]

                try:
                    prots_dict[genome] += 1
                except Exception:
                    prots_dict[genome] = 1
                
            else:
                seq += line.strip("\n")
        db.setdefault(genome,[]).append((f'{genome}_{num}',seq))
        try:
            prots_dict[genome] += 1
        except Exception:
            prots_dict[genome] = 1

        for key in uniques.keys():
            identifier = uniques[key][0]

            for k in key[1:]:
                try:
                    proteins = db[k]
                except Exception:
                    proteins = []
                    # no proteins
                
                if proteins != []:
                    with open(f'{folder}vRhyme_unique_bins_proteins/{identifier}.faa', 'a') as outfile:
                        for prot in proteins:
                            outfile.write(f'>{prot[0]}\n{prot[1]}\n')            
    
    db = None
    return prots_dict


def final_fasta(in_fasta, in_genes, in_prots, folder, prefix, final):
    '''
    After identifying best iteration, write out all fasta files per bin
    '''
    subprocess.run('mkdir ' + folder + "vRhyme_best_bins_fasta", shell=True)
    with open(f'{folder}vRhyme_best_bins.{final}.membership.tsv', "r") as bins:
        bins = bins.read().replace("\n","\t").split("\t")
        if bins[-1] == '':
            bins = bins[:-1]
        bins = {bins[i]:bins[i+1] for i in range(2,len(bins),2)}
    with open(in_fasta, 'r') as fasta:
        seq = ''
        for line in fasta:
            if line.startswith(">"):
                if seq != '':
                    try:
                        bin = bins[header]
                        p = prefix.replace("#",bin)
                        with open(folder + "vRhyme_best_bins_fasta/vRhyme_bin_" + bin + ".fasta", "a") as outfile:
                            outfile.write(">" + p + header + "\n" + seq + "\n")
                    except Exception:
                        pass

                seq = ''
                header = line[1:].strip("\n") # no ">", no "\n"

            else:
                seq += line.strip("\n")

        # last one
        try:
            bin = bins[header]
            p = prefix.replace("#",bin)
            with open(folder + "vRhyme_best_bins_fasta/vRhyme_bin_" + bin + ".fasta", "a") as outfile:
                outfile.write(">" + p + header + "\n" + seq + "\n")
        except Exception:
            pass

    with open(in_genes, 'r') as fasta:
        seq = ''
        for line in fasta:
            if line.startswith(">"):
                if seq != '':
                    try:
                        bin = bins[split_header]
                        p = prefix.replace("#",bin)
                        with open(folder + "vRhyme_best_bins_fasta/vRhyme_bin_" + bin + ".ffn", "a") as outfile:
                            outfile.write(">" + p + header + "\n" + seq + "\n")
                    except Exception:
                        pass

                seq = ''
                header = line[1:].strip("\n") # no ">", no "\n"
                split_header = header.split(" # ",1)[0]
                split_header = split_header.rsplit("_",1)[0]

            else:
                seq += line.strip("\n")

        # last one
        try:
            bin = bins[split_header]
            p = prefix.replace("#",bin)
            with open(folder + "vRhyme_best_bins_fasta/vRhyme_bin_" + bin + ".ffn", "a") as outfile:
                outfile.write(">" + p + header + "\n" + seq + "\n")
        except Exception:
            pass

    with open(in_prots, 'r') as fasta:
        seq = ''
        for line in fasta:
            if line.startswith(">"):
                if seq != '':
                    try:
                        bin = bins[split_header]
                        p = prefix.replace("#",bin)
                        with open(folder + "vRhyme_best_bins_fasta/vRhyme_bin_" + bin + ".faa", "a") as outfile:
                            outfile.write(">" + p + header + "\n" + seq + "\n")
                    except Exception:
                        pass

                seq = ''
                header = line[1:].strip("\n") # no ">", no "\n"
                split_header = header.split(" # ",1)[0]
                split_header = split_header.rsplit("_",1)[0]

            else:
                seq += line.strip("\n")

        # last one
        try:
            bin = bins[split_header]
            p = prefix.replace("#",bin)
            with open(folder + "vRhyme_best_bins_fasta/vRhyme_bin_" + bin + ".faa", "a") as outfile:
                outfile.write(">" + p + header + "\n" + seq + "\n")
        except Exception:
            pass


#
#
#
