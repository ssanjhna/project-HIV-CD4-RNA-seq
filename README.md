# project-HIV-CD4-RNA-seq
This project analyses RNA-seq data from CD4⁺ Effector Memory T cells under different stimulation conditions (Unstimulated, PMA, IL15, Bryostatin) to identify differentially expressed genes and enriched biological pathways associated with latency reversal.

ENA_metadata.tsv contains the accession IDs (run_accession	experiment_accession), per sample details including type of T cell, condition species and experiment type (experiment_title) and the fastq_ftp link to download files.

Renaming.sh contains the script that uses the ENA metadata file to rename SRR file names into descriptive names.

Multiqc_report_pretrim.html has the compiled QC analysis of the original fastq files prior to trimming. 
See Project-FinalReport.Rmd for details of trimming.
Multiqc_report_trimmed.html has the compiled QC analysis of the original fastq files post trimgalore.

Build_index_human.sh contains the script to build the STAR index for alignment. The genome used is hg38 GRCh38.p14 GCF_000001405.40 with accompanying .gtf RefSeq file.

Star.script runs the alignment using the STAR aligner. Requires fastq files and the indexed genome.

Qorts.sh runs Qorts on the aligned.bam files requiring the genommic.gtf and aligned files.
Multiqc_report_qorts.html has the compiled qorts analysis of the aligned files.

featurecountsall.sh runs featureCounts on the aligned files compiling them into 1 .txt.summary and detailed .txt files.
featureCounts_genes_all_samples.zip contains the output files of the featurecountsall.sh script.

Project_DEG_Analysis_script.R contains all the processing of feature counts data, DESeq2 object generation, normalisation, DEG analysis, GO and GSEA enrichment analysis code chunks.

DESeq2results_AllConditions.zip contains the raw tables of genes that are significantly differentially expressed (adjusted p-value < 0.05) between each treatment (PMA, IL15, Bryostatin) and unstimulated cells, including gene annotations (gene symbol, gene name, Entrez ID) and differential expression statistics (log2 fold change, p-values). See .R file for code generating these tables.

Go_txt.zip contains the results of Gene Ontology (GO) enrichment analysis (see .R file for code generating these .txt files) for genes that are differentially expressed in each condition.

Project_DEG_Analysis_script.R contains the R script for the DESeq2 object generation, normalisation and plot generation.

Project-FinalReport.Rmd and Project-FinalReport.html contain the full analysis pipeline and summarisation of results.
