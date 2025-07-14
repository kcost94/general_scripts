#!/usr/bin/env bash

# This script performs loops through the BED files in the 'sep-chroms' directory,
# processes each file to extract sequences, and generates CpG profiles.
for file in ../sep-chroms/chr*.bed
do
name=`basename $file .bed`

python step.py $file > $name.bed
bedtools getfasta -fi ../sep-chroms/chr*.fa -bed $name -fo $name.fa
python CpG_calling.py $name.fa > $name.txt

done
