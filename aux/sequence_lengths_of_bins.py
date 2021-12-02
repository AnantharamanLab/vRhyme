#! /usr/bin/env python3
# vRhyme
# Author: Kristopher Kieft
# University of Wisconsin-Madison


import argparse
import os
import sys

def get_lengths(infolder, output, extension, len_ext):
    lengths = []
    files = os.listdir(infolder)
    files = [i for i in files if i[-len_ext:] == extension]

    if len(files) == 0:
        sys.stderr.write("\nError: No input files were identified. Verify that the input (-i) and extension (-e) are correct. Exiting.\n\n")
        exit()

    for f in files:
        file = infolder + "/" + f
        base = f.rsplit(".",1)[0]
        with open(file, 'r') as fasta:
            seq = 0
            for line in fasta:
                if not line.startswith(">"):
                    seq += len(line.strip("\n"))
            lengths.append((base,seq))

    lengths.sort(key=lambda x: x[1], reverse=True)

    with open(output, 'w') as outfile:
        outfile.write('bin\tlength\n')
        for item in lengths:
            key,val = item
            outfile.write(f'{key}\t{val}\n')

if __name__ == '__main__':
    descript = """
    Calculate the total nucldeotide length of each bin.
    Input unlinked bins (not from link_bin_sequences.py).
    """
    vRhyme = argparse.ArgumentParser(description=descript, formatter_class=argparse.RawTextHelpFormatter)
    vRhyme.add_argument('-i', metavar='', type=str, nargs=1, required=True, help="folder containing input bin fasta files")
    vRhyme.add_argument('-o', metavar='', type=str, nargs=1, required=True, help='name of output tab-separated file')
    vRhyme.add_argument('-e', metavar='', type=str, nargs=1, default = ['fasta'], help='extension of files to count lengths [fasta]')

    #
    args = vRhyme.parse_args()
    infolder = str(args.i[0])
    output = str(args.o[0])
    extension = str(args.e[0])
    len_ext = len(extension)
    if infolder[-1] != '/':
        infolder += '/'
    #
    #
    if os.path.exists(output):
        sys.stderr.write("\nError: The output file already exists. Exiting.\n\n")
        exit()
    #
    get_lengths(infolder, output, extension, len_ext)

#
#
#
