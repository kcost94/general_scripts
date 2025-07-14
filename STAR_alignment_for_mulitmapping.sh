#!/usr/bin/env bash

# This script performs STAR alignment for multi-mapping reads.

name=$1
fastq_dir=$2
genome_dir=$3
gtf_file=$4

# Example usage: ./STAR_alignment_for_mulitmapping.sh sample_name /path/to/fastq /path/to/genome_dir /path/to/gtf_file

# run STAR alignment
STAR --runThreadN 8 --outFileNamePrefix $name \
    --outSAMtype BAM SortedByCoordinate --outSAMattributes All --outSAMstrandField intronMotif \
    --genomeDir $genome_dir \
    --readFilesIn $fastq_dir/$name\.RNAseq.R1.fastq.gz $fastq_dir/$name\.RNAseq.R2.fastq.gz \
    --readFilesCommand zcat --outFilterMultimapNmax 100 --winAnchorMultimapNmax 100 \
    --sjdbGTFfile $gtf_file

mv ./$name\Aligned.sortedByCoord.out.bam $name.bam

# Convert BAM to BED format
samtools sort -@ 8 -T $name\.sorted -o $name-sorted.bam $name.bam
samtools index $name-sorted.bam

# Remove duplicates using Picard
java -jar /opt/easybuild/software/picard/2.18.4-Java-1.8.0_131/picard.jar MarkDuplicates INPUT=$name-sorted.bam OUTPUT=$name-nr-sorted.bam METRICS_FILE=picard-metrics.txt VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true REMOVE_DUPLICATES=true
samtools sort -@ 8 -T $name.sorted -o sorted-$name-nr-sorted.bam $name-nr-sorted.bam
samtools index sorted-$name-nr-sorted.bam

#rm $name.bam $name-sorted.bam $name-nr-sorted.bam

# make coverage files
genomeCoverageBed -bg -ibam sorted-$name-nr-sorted.bam -g $genome_dir/mm10.chrom.sizes  > $name.bedgraph
wigToBigWig $name.bedgraph $genome_dir/mm10.chrom.sizes $name.bw

