#!/bin/bash

tsv="/athena/angsd/scratch/sas4096/project/ENA_metadata.tsv"

folders=(
    "/athena/angsd/scratch/sas4096/project/FASTQ/trimmed"
    "/athena/angsd/scratch/sas4096/project/star_alignment/aligned"
    "/athena/angsd/scratch/sas4096/project/FASTQ/original_fastqs"
)

declare -A map
while IFS=$'\t' read -r run_accession experiment_accession experiment_title fastq_ftp; do
    [[ "$run_accession" == "run_accession" ]] && continue
    # Extract clean sample name: after second colon, remove 'Homo sapiens RNA-Seq', replace spaces with underscores
    sample_name=$(echo "$experiment_title" | awk -F': ' '{print $3}' | sed 's/ Homo sapiens RNA-Seq//g' | tr ' ' '_')
    map["$run_accession"]="$sample_name"
done < "$tsv"

echo "===== RENAMING FILES ====="
for dir in "${folders[@]}"; do
    for f in "$dir"/SRR*; do
        filename=$(basename "$f")
        # extract SRR ID at start of filename
        run=$(echo "$filename" | grep -o '^SRR[0-9]\+')
        if [ -n "$run" ] && [ -n "${map[$run]}" ]; then
            new_file="$dir/${filename/$run/${map[$run]}}"
            echo "Renaming: $f -> $new_file"
            mv "$f" "$new_file"
        else
            echo "No mapping found for $filename"
        fi
    done
done
echo "DONE"