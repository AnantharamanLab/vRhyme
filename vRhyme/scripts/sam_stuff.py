#! /usr/bin/env python3
# vRhyme
# Author: Kristopher Kieft
# University of Wisconsin-Madison

import sys
import subprocess


def sam_stuff(sam, bam):
    subprocess.run("samtools view -@ 1 -h -b " + sam + " > " + bam + " 2> /dev/null", shell=True)


if __name__ == '__main__':
    sam = sys.argv[1]
    bam = sys.argv[2]
    sam_stuff(sam,bam)
