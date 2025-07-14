#!/usr/bin/env python
import sys

# This script calculates CpG scores, density, and CG content from a given input file.
# It expects the input file to be in FASTA format, where each sequence starts with a line beginning with '>'.
# Usage: python CpG-calling.py <input_file>

f = open(sys.argv[1], 'r')
for line in f:
    line=line.strip()
    sline=line.split()
    for n in range(0,int(sline[2]) - 500 ,250):
        print(sline[0] + "\t" +str(n) + "\t" + str(n +500) )

