#!/usr/bin/env bash

# This script performs HISAT2 alignment for RNA sequencing data.
name=$1
genome_dir=$2
genome=$3
fastq_dir=$4

# Example usage: ./HISat2_for_RNA_alignment.sh sample_name /path/to/genome_dir genome_name /path/to/fastq_dir 

#module load HISAT2
hisat2 -p 8 --dta --rf --known-splicesite-infile $genome_dir/splice-sites.txt -x $genome_dir/HISTAT2-indexes/$genome \
    -1 $fastq_dir/$name.RNAseq.R1.fastq.gz \
    -2 $fastq_dir/$name.RNAseq.R2.fastq.gz \
    -S $name.sam

# Convert SAM to BAM, sort, and index
module load samtools
samtools sort -@ 8 -T $name\.sorted -o $name-sorted.bam $name.sam
samtools index $name-sorted.bam
rm $name.sam
 
# Remove duplicates using Picard
module load picard
java -jar /opt/easybuild/software/picard/2.18.4-Java-1.8.0_131/picard.jar MarkDuplicates INPUT=$name-sorted.bam OUTPUT=$name-nr-sorted.bam METRICS_FILE=picard-metrics.txt VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true REMOVE_DUPLICATES=true

# Sort the non-redundant BAM file and create name-sorted BAM
samtools sort -@ 8 -T $name.sorted -o sorted-$name-nr-sorted.bam $name-nr-sorted.bam
samtools index $name-nr-sorted.bam
samtools sort -n -@ 8 -T $name.sorted -o name-sorted-$name-nr-sorted.bam sorted-$name-nr-sorted.bam

# Convert BAM to BED format and generate coverage files 
module load bedtools
bamToBed -i sorted-$name-nr-sorted.bam > $name-RNA.bed
genomeCoverageBed -bg -ibam sorted-$name-nr-sorted.bam -g $genome_dir/$genome.chrom.sizes > $name.bedgraph
wigToBigWig $name.bedgraph $genome_dir/$genome.chrom.sizes $name.bw

# Count reads in exonic regions
coverageBed -b $name-RNA.bed -a $genome_dir/$genome.gtf -counts > $name.counts
python get-multicov-exonic-counts.py $name.counts

