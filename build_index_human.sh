#!/bin/bash
#SBATCH --job-name=STAR_index
#SBATCH --output=/athena/angsd/scratch/sas4096/project/star_alignment/star_%j.out
#SBATCH --error=/athena/angsd/scratch/sas4096/project/star_alignment/star_%j.err
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=40G
#SBATCH --partition=angsd_class
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=sas4096@med.cornell.edu

source /home/fs01/sas4096/miniforge3/etc/profile.d/conda.sh
conda activate angsd

STAR --runMode genomeGenerate \
     --runThreadN 16 \
     --genomeDir /athena/angsd/scratch/sas4096/project/star_alignment/index \
     --genomeFastaFiles /athena/angsd/scratch/sas4096/project/star_alignment/ncbi_dataset/data/GCF_000001405.40/GCF_000001405.40_GRCh38.p14_genomic.fna \
     --sjdbGTFfile /athena/angsd/scratch/sas4096/project/star_alignment/ncbi_dataset/data/GCF_000001405.40/genomic.gtf \
     --sjdbOverhang 99