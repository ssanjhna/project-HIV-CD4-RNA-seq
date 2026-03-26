#!/bin/bash
#SBATCH --job-name=STAR_align
#SBATCH --output=/athena/angsd/scratch/sas4096/project/STAR/star_%j.out
#SBATCH --error=/athena/angsd/scratch/sas4096/project/STAR/star_%j.err
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=40G
#SBATCH --partition=angsd_class
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=sas4096@med.cornell.edu

source /home/fs01/sas4096/miniforge3/etc/profile.d/conda.sh
conda activate angsd

THREADS=16
INPUT_DIR="/athena/angsd/scratch/sas4096/project/FASTQ/trimmed"
OUTPUT_DIR="/athena/angsd/scratch/sas4096/project/star_alignment/aligned"
STAR_INDEX="/athena/angsd/scratch/sas4096/project/star_alignment/index"

for fq in $INPUT_DIR/*_trimmed.fq.gz
do
    base=$(basename $fq _trimmed.fq.gz)
    echo "Processing $base ..."

    STAR --runMode alignReads \
         --runThreadN $THREADS \
         --genomeDir $STAR_INDEX \
         --readFilesIn $fq \
         --readFilesCommand zcat \
         --outFileNamePrefix $OUTPUT_DIR/${base}_ \
         --outSAMtype BAM SortedByCoordinate 

    samtools index $OUTPUT_DIR/${base}_Aligned.sortedByCoord.out.bam

    echo "Done with $base"
done