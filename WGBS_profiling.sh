#!/usr/bin/env bash

# This script performs Bismark alignment for whole-genome bisulfite sequencing (WGBS) profiling.
name=$1
genome_dir=$2
fastq_dir=$3
# Example usage: ./WGBS_profiling.sh sample_name /path/to/genome_dir /path/to/fastq_dir

module load cutadapt/1.9.1-foss-2016b-Python-2.7.12
module load trimgalore/0.5.0

  trim_galore --paired --trim1 $fastq_dir/$name\_1.fastq $fastq_dir/$name\_2.fastq
  rm $fastq_dir/$name\_1.fastq $fastq_dir/$name\_2.fastq

module load Bismark/0.19.1
  bismark -p 5 --non_bs_mm -o $name $genome_dir -1 $fastq_dir/$name\_1_val_1.fq -2 $fastq_dir/$name\_2_val_2.fq

 cd $name/
  deduplicate_bismark --output_dir ./ --paired ./$name\_1_val_1_bismark_bt2_pe.bam
  bismark_methylation_extractor --parallel 5 --no_overlap --gzip -p --bedGraph --zero_based --buffer_size 25G --ignore 3 --ignore_r2 3 --cytosine_report  --genome_folder $genome_dir ./$name\_1_val_1_bismark_bt2_pe.deduplicated.sam


