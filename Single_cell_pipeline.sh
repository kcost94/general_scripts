#!/usr/bin/env bash
# This script runs the Cell Ranger ATAC pipeline for single-cell ATAC-seq data.

# Usage: ./Single_cell_pipeline.sh <sample_name> <fastq_directory> <reference_directory> <barcode_file> <peaks_file> <fragments_file>
sample_name=$1
fastq_directory=$2
reference_directory=$3
barcode_file=$4
peaks_file=$5   
fragments_file=$6   


module load cellranger/1.1.0-atac
## original run
cellranger-atac count --id=WT-rep1-mm10 --sample=WT_rep1 --fastqs=$fastq_directory --reference=$reference_directory


## reanalyze with specific parameters
# Create a parameters CSV file for reanalyzing with specific parameters
cellranger-atac reanalyze --id=WT-rep1-promoters-mm10 --params parameters.csv \
    --barcodes=$barcode_file --peaks=$peaks_file \
    --fragments=$fragments_file \
    --reference=/opt/cellranger/refdata-cellranger-atac-mm10-1.0.0


