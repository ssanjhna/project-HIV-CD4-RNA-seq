#!/bin/bash
#SBATCH --job-name=qorts
#SBATCH --output=/athena/angsd/scratch/sas4096/project/star_alignment/aligned/qorts_%j.out
#SBATCH --error=/athena/angsd/scratch/sas4096/project/star_alignment/aligned/qorts_%j.err
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=40G
#SBATCH --partition=angsd_class
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=sas4096@med.cornell.edu

source /home/fs01/sas4096/miniforge3/etc/profile.d/conda.sh
conda activate qorts

GTF=/athena/angsd/scratch/sas4096/project/star_alignment/ncbi_dataset/data/GCF_000001405.40/genomic.gtf
BAMDIR=/athena/angsd/scratch/sas4096/project/star_alignment/aligned
OUTDIR=/athena/angsd/scratch/sas4096/project/star_alignment/aligned/qorts
QORTS=/athena/angsd/scratch/mef3005/share/envs/qorts/share/qorts-1.3.6-1/QoRTs.jar

for BAM in $BAMDIR/*.bam; do
    SAMPLE=$(basename $BAM .bam)
    java -Xmx8g -jar $QORTS QC \
        --singleEnded \
        --generatePlots \
        --maxReadLength 100 \
        $BAM \
        $GTF \
        $OUTDIR/$SAMPLE/
done