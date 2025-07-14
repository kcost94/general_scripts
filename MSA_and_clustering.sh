#!/usr/bin/env bash

# This script performs MSA and clustering on a given sequence file.
name=$1
fasta=$2

# Example usage: ./MSA_and_clustering.sh sample_name /path/to/fasta
# this script requires MAFFT, trimAl, seqret, and PhyML to be installed.
# this script also requires R with ggtree and tidytree packages installed.
# the fasta file needs to have names shorter than 10 characters to avoid issues with the phylip format.

mkdir $name
cd $name

# run MSA using MAFFT
mafft $fasta > mafft-$name.fa

# run trimming using trimAl
trimal -in mafft-$name.fa -out trimal-$name.fa -automated1

# convert to phylip format
seqret -sequence trimal-$name.fa -outseq $name.phylip -osformat phylip

# run phylogenetic analysis using PhyML
phyml -i $name.phylip --run_id $name --quiet --no_memory_check > temp-$name.out

#run tree visualization
cp $name.phylip_phyml_tree_$name.txt in-tree.txt ; cp trimal-$name.fa in-fasta.txt
Rscript tree_grower.r ; mv outfile.png $name-tree.png

# Clean up temporary files
rm mafft-$name.fa trimal-$name.fa $name.phylip temp-$name.out in-tree.txt in-fasta.txt

cd ../

