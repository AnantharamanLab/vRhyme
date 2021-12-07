#! /usr/bin/env python3
# vRhyme
# Author: Kristopher Kieft
# University of Wisconsin-Madison


import argparse
import os
import subprocess
import sys

def linker(infolder, output, extension, len_ext, n, char):
    files = os.listdir(infolder)
    files = [i for i in files if i[-len_ext:] == extension]

    if len(files) == 0:
        sys.stderr.write("\nError: No input files were identified. Verify that the input (-i) and extension (-e) are correct. Exiting.\n\n")
        exit()
    
    subprocess.run('mkdir ' + output, shell=True)

    N_string = ''.join([char] * n)
    for f in files:
        file = infolder + "/" + f
        base = f.rsplit(".",1)[0]
        with open(file, 'r') as fasta:
            seq = ''
            for line in fasta:
                if not line.startswith(">"):
                    seq += line.strip("\n")
                    seq += N_string

        with open(output + base + '.linked.' + extension, "w") as outfile:
            outfile.write(">" + base + "\n" + seq[:-n] + "\n") # -n to remove last n added characters


if __name__ == '__main__':
    descript = """
    Link scaffolds to generate a single sequence per bin (e.g., useful for CheckV input).
    Note: the order of linked sequences is arbitrary and does not reflect true connection."""
    vRhyme = argparse.ArgumentParser(description=descript, formatter_class=argparse.RawTextHelpFormatter)
    vRhyme.add_argument('-i', metavar='', type=str, nargs=1, required=True, help="folder containing input bin fasta files")
    vRhyme.add_argument('-o', metavar='', type=str, nargs=1, required=True, help='folder to deposit linked bin sequence files')
    vRhyme.add_argument('-e', metavar='', type=str, nargs=1, default = ['fasta'], help='extension of files to link [fasta]')
    vRhyme.add_argument('-n', metavar='', type=str, nargs=1, default = ['1500'], help='number of characters to add as linker [1500]')
    vRhyme.add_argument('-c', metavar='', type=str, nargs=1, default = ['N'], help='character to use as linker [N]')
    #
    args = vRhyme.parse_args()
    infolder = str(args.i[0])
    output = str(args.o[0])
    extension = str(args.e[0])
    len_ext = len(extension)
    n = int(args.n[0])
    char = str(args.c[0])
    if len(char) > 1:
        sys.stderr.write("\nError: The linker character should be a single character. Exiting.\n\n")
        exit()
    if n < 1000:
        sys.stdout.write(f"\nCaution: appending fewer than 1000 linker characters ({n}) may cause downstream analyses to artificailly bridge the link.\n\n")
    if infolder[-1] != '/':
        infolder += '/'
    if output[-1] != '/':
        output += '/'
    #
    #
    if os.path.exists(output):
        sys.stderr.write("\nError: The output folder already exists. Exiting.\n\n")
        exit()

    #
    linker(infolder, output, extension, len_ext, n, char)

#
#
#
