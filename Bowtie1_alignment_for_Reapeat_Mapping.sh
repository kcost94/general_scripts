#!/usr/bin/env bash


# This script performs Bowtie1 alignment for unique mapping to profile repeats.
line=$1
name=$2
# Example usage: ./Bowtie1_alignment_for_Reapeat_Mapping.sh sample_name name_of_fastq_files

mkdir $line
cd $line
ln -s ../fastq/*$name* ./

# Load necessary modules
module purge
module load cutadapt/1.9.1-foss-2016b-Python-2.7.12
module load trimgalore/0.5.0

# Trim adapters and low-quality bases
trim_galore --length 10 --trim1 --paired *_R1_001.fastq *_R2_001.fastq

module load foss/2018b
module load Bowtie/1.2.3-foss-2018b

# Align reads to the reference genome using Bowtie
bowtie --quiet -m 1 -p 8 --tryhard --chunkmbs 500 /net/isi-dcnl/ifs/user_data/dschones/bioresearch/UCSC/mm9/bowtie-indexes/mm9 -1 *_val_1.fq -2 *_val_2.fq -S $line\.sam

# Convert SAM to BAM, sort, and index
module load samtools
samtools sort -@ 8 -T $line\.sorted -o $line-sorted.bam $line\.sam
samtools index $line-sorted.bam
rm $line\.sam

# Remove duplicates using Picard
module load picard
java -jar /opt/easybuild/software/picard/2.18.4-Java-1.8.0_131/picard.jar MarkDuplicates INPUT=$line-sorted.bam OUTPUT=$line-nr-sorted.bam METRICS_FILE=picard-metrics.txt VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true REMOVE_DUPLICATES=true

# Filter out unmapped reads and sort by name
samtools sort -@ 8 -T $line\.sorted -o sorted-$line-nr-sorted.bam $line-nr-sorted.bam
samtools view -F 3844 -f 2 -bh -o filtered-$line-nr-sorted.bam $line-nr-sorted.bam

samtools sort -n filtered-$line-nr-sorted.bam > filtered-name-sorted-$line.bed
bamToBed -bedpe -i filtered-name-sorted-$line.bed > $line.bed

# Generate coverage files
module load bedtools
bedtools genomecov -pc -bg -ibam filtered-$line-nr-sorted.bam -g /net/isi-dcnl/ifs/user_data/dschones/bioresearch/UCSC/mm9/mm9.chrom.sizes  > $line\.bg
module load UCSC/20180213
wigToBigWig $line\.bg /net/isi-dcnl/ifs/user_data/dschones/bioresearch/UCSC/mm9/mm9.chrom.sizes $line\.bw

# Call peaks using MACS2
module load MACS2
macs2 callpeak -t filtered-$line-nr-sorted.bam -g mm -f BAMPE -n $line

# Clean up intermediate files
rm *_val_1.fq *_val_2.fq $line\.bg $line-nr-sorted.bam filtered-$line-nr-sorted.bam picard-metrics.txt
cd ../

