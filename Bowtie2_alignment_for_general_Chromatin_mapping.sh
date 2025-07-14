#!/usr/bin/env bash

line=$1
name=$2

mkdir bt2-$line
cd bt2-$line

ln -s ../fastq/*$name* ./
# Trimming adapters and low-quality bases

trim_galore --paired --length 10 --trim1 --gzip *_R1_001.fastq *_R2_001.fastq

# Load Bowtie2 for alignment
export BOWTIE2_INDEXES=/data/UCSC/mm9/bowtie2-indexes
# Align reads to the reference genome using Bowtie2
bowtie2 -p 4 --local --very-sensitive-local --phred33 -I 1 -X 500 -x $BOWTIE2_INDEXES/mm9 -1 *_val_1.fq.gz -2 *_val_2.fq.gz -S $line\.sam > align-out.txt 2>&1

# Convert SAM to BAM, sort, and index
samtools sort -@ 4 -T $line\.sorted -o $line-sorted.bam $line.sam
samtools index $line-sorted.bam
rm $line\.sam
samtools flagstat $line-sorted.bam > stat1

# Remove duplicates using Picard
java -jar /opt/easybuild/software/picard/2.18.4-Java-1.8.0_131/picard.jar MarkDuplicates INPUT=$line-sorted.bam OUTPUT=$line-nr-sorted.bam METRICS_FILE=picard-metrics.txt VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true REMOVE_DUPLICATES=true

# Filter out unmapped reads and sort by name
samtools view -hf 0x2 $line-sorted.bam | grep -v "XS:i:" > filtered-$line.sam
samtools sort -@ 8 -n -T $line-filtered.sorted -o filtered-$line-name-sorted.bam filtered-$line.sam
samtools sort -@ 4 -n -T $line.sorted -o $line-name-sorted.bam $line-sorted.bam

# Convert BAM to BEDPE format
bamToBed -bedpe -i filtered-$line-name-sorted.bam > filtered-$line\.bedpe
bamToBed -bedpe -i $line-name-sorted.bam > $line\.bedpe

# filter BEDPE files
awk '{ if ($5 >0 ) print $0 }' filtered-$line\.bedpe | awk '{ if ($2 < $5 ) start = $2 ; else start = $5 } { if ($3 > $6 ) end = $3 ; else end = $6 } { if ( end -start > 10) print $1 "\t" start "\t" end  }' > filtered-$line\.bed
awk '{ if ($5 >0 ) print $0 }' $line\.bedpe | awk '{ if ($2 < $5 ) start = $2 ; else start = $5 } { if ($3 > $6 ) end = $3 ; else end = $6 } { if ( end -start > 10) print $1 "\t" start "\t" end  }' > $line\.bed

# Generate bigWig files from BED files
bedtools genomecov -bg -i filtered-$line\.bed -g /net/isi-dcnl/ifs/user_data/dschones/bioresearch/UCSC/mm9/mm9.chrom.sizes > filtered-$line\.bg
module load UCSC/20180213
bedGraphToBigWig filtered-$line\.bg /net/isi-dcnl/ifs/user_data/dschones/bioresearch/UCSC/mm9/mm9.chrom.sizes filtered-$line\.bw

bedtools genomecov -bg -i $line\.bed -g /net/isi-dcnl/ifs/user_data/dschones/bioresearch/UCSC/mm9/mm9.chrom.sizes > $line\.bg
module load UCSC/20180213
bedGraphToBigWig $line\.bg /net/isi-dcnl/ifs/user_data/dschones/bioresearch/UCSC/mm9/mm9.chrom.sizes $line\.bw

# Call peaks using MACS2
module load MACS2
macs2 callpeak -t filtered-$line\.bed -g mm --keep-dup all -n $line
macs2 callpeak -t $line\.bed -g mm --keep-dup all -n $line

# Clean up intermediate files
rm *_val_1.fq.gz *_val_2.fq.gz $line\.bg $line-nr-sorted.bam filtered-$line-nr-sorted.bam picard-metrics.txt
rm filtered-$line\.bed filtered-$line\.bg  $line\.bg
rm filtered-$line-name-sorted.bam filtered-$line.sam $line-name-sorted.bam $line\.sam


cd ../
