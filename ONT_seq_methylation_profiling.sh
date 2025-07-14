#!/usr/bin/env bash

### This script is designed for ONT sequencing methylation profiling.
### It requires the input to be fast5 or pod5 files and the genome reference file.

#Input is expected to be the path to a fast5/pod5 files name is expected to be the name of the sample 
file=$1
name=$2
# The genome reference file should be provided as the third argument
genome=$3

# Example: ./ONT_seq_methylation_profiling.sh path_to_pod5s/ sample_name path_to_genome/genome.fa


### this should be run on the command line outside the script, run with GPU 
~/rds/rds-acf1004-afs-lab-rds/programs/dorado-0.6.2-linux-x64/bin/dorado basecaller --reference $genome -kit-name SQK-NBD114-24 --emit-moves hac,5mCG_5hmCG $file -r > $name-filtered.bam 

# seperate the file by barcodes
mkdir sep
~/rds/rds-acf1004-afs-lab-rds/programs/dorado-0.6.2-linux-x64/bin/dorado demux --output-dir sep --no-classify calls-ref-filtered.bam

### do this on the command line outside the script, run with CPU
# Activate the DeepMod2 environment
conda activate deepmod2

# Process each BAM file in the 'sep' directory
for bam in sep/*.bam 
 do 
 bam_name=`basename $file .bam` 
 python ~/rds/rds-acf1004-afs-lab-rds/programs/DeepMod2/deepmod2 detect --bam $bam --input $file --model bilstm_r10.4.1_4khz_v4.1 --file_type pod5 --threads 25 --ref $genome --output $bam_name --seq_type dna

 done

