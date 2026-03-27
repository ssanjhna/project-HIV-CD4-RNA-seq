#!/bin/bash
#SBATCH --job-name=feature_counts
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=40G
#SBATCH --partition=angsd_class
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=sas4096@med.cornell.edu

source /home/fs01/sas4096/miniforge3/etc/profile.d/conda.sh
conda activate angsd

#!/bin/bash

BAM_DIR=/athena/angsd/scratch/sas4096/project/star_alignment/aligned
GTF=/athena/angsd/scratch/sas4096/project/star_alignment/ncbi_dataset/data/GCF_000001405.40/genomic.gtf
OUT_DIR=$BAM_DIR/featurecounts_output_all

mkdir -p "$OUT_DIR"

featureCounts \
    -a "$GTF" \
    -o "$OUT_DIR/gene_counts_all_samples.txt" \
    -s 0 \
    -T 16 \
    -M --fraction \
    -Q 0 \
    $BAM_DIR/*Aligned.sortedByCoord.out.bam