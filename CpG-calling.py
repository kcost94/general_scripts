#!/usr/bin/env python
import sys

# This script calculates CpG scores, density, and CG content from a given input file.
# It expects the input file to be in FASTA format, where each sequence starts with a line beginning with '>'.
# Usage: python CpG-calling.py <input_file>

f = open(sys.argv[1], 'r')
for line in f:
    line=line.strip()
    sline=line.split()
    # Initialize variables for each fasta entry
    if sline[0][0] == ">":
        name=sline[0]
        temp="NA"
        CpG=0
        C=0
        G=0
        density=0
        score=0
        CG=0
    else:
        # Process the sequence line
        for n in range(len(sline[0])):
            if (temp + sline[0][n]).lower() == "cg":
                CpG+=1
            temp=sline[0][n]
            if (sline[0][n]).lower() == "g":
                G+=1
            if (sline[0][n]).lower() == "c":
                C+=1
        # Calculate scores and densities if C and G are greater than 0
        if C > 0 and G > 0:
            score=float(CpG*len(sline[0]))/float(G*C)
            CG=(G+C)/len(sline[0])
            density=float(CpG)/float(len(sline[0]))
        print(name +"\t" + str(score) + "\t" +str(CpG) + "\t" + str(density) + "\t" + str(CG))

f.close()

