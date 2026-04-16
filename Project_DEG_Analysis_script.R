library(ggplot2)
library(tidyr)
library(dplyr)

df <- read.table(
  "D:/data-sanjhna/Cornell/ANGSD/Project/featurecounts/gene_counts_all_samples.txt.summary",
  header = TRUE,
  check.names = FALSE
)

colnames(df) <- c("Status", gsub(".*/(.*?)_Aligned.*", "\\1", colnames(df)[-1]))

combined_df <- df %>%
  pivot_longer(-Status, names_to = "sample", values_to = "count")

combined_df %>%
  filter(count > 0) %>%
  ggplot(aes(x = count, y = sample, fill = Status)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  labs(
    title = "FeatureCounts Assignment Summary",
    x = "Read counts",
    y = "Sample"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

library(tidyr)
library(dplyr)
total_per_sample <- combined_df %>%
  group_by(sample) %>%
  summarise(total = sum(count))

intronic_counts <- combined_df %>%
  filter(Status == "Unassigned_NoFeatures") %>%
  rename(intronic = count) %>%
  select(sample, intronic)

fractions <- intronic_counts %>%
  left_join(total_per_sample, by = "sample") %>%
  mutate(intronic_fraction = intronic / total)

fractions <- fractions %>%
  mutate(condition = case_when(
    grepl("unstim", sample, ignore.case = TRUE) ~ "Unstimulated",
    grepl("PMA",   sample, ignore.case = TRUE) ~ "PMA",
    grepl("IL15",  sample, ignore.case = TRUE) ~ "IL15",
    grepl("bryo",  sample, ignore.case = TRUE) ~ "Bryostatin",
    TRUE ~ "Unknown"
  ))

kw_result <- kruskal.test(intronic_fraction ~ condition, data = fractions)
print(kw_result)

pairwise_result <- pairwise.wilcox.test(
  fractions$intronic_fraction,
  fractions$condition,
  p.adjust.method = "BH"
)
print(pairwise_result)

ggplot(fractions, aes(x = condition, y = intronic_fraction, fill = condition)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 2) +
  labs(
    title = "Fraction of Reads Mapping to Intronic Regions",
    x = "Condition",
    y = "Intronic Read Fraction"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

count_table <- read.table(
  "D:/data-sanjhna/Cornell/ANGSD/Project/featurecounts/gene_counts.txt",
  header = TRUE,
  skip = 1,
  check.names = FALSE
)

count_table <- count_table[, c(1, 7:ncol(count_table))]

colnames(count_table) <- c("GeneID", gsub(".*/(.*?)_Aligned.*", "\\1", colnames(count_table)[-1]))

rownames(count_table) <- count_table$GeneID
count_table <- count_table[, -1]

head(count_table)
dim(count_table)

# reading in the feature counts results and changing name
library(DESeq2)
library(vsn)
library(dplyr)

count_table <- read.table(
  "D:/data-sanjhna/Cornell/ANGSD/Project/featurecounts/gene_counts.txt",
  header = TRUE,
  skip = 1,
  check.names = FALSE
)

colnames(count_table) <- c(
  colnames(count_table)[1:6],
  gsub(".*/(EM_.*?)_Aligned.*", "\\1", colnames(count_table)[7:ncol(count_table)])
)

row.names(count_table) <- make.names(count_table$Geneid)
cts <- count_table[, 7:ncol(count_table)]

df_coldata <- fractions %>%
  select(sample, condition) %>%
  mutate(condition = factor(condition,
                            levels = c("Unstimulated", "PMA", "IL15", "Bryostatin"))) %>%
  tibble::column_to_rownames("sample")

dds <- DESeqDataSetFromMatrix(
  countData = cts[, rownames(df_coldata)],
  colData   = df_coldata,
  design    = ~ condition
)
rowData(dds) <- count_table[, 1:6]
head(counts(dds))
colSums(counts(dds))
colSums(counts(dds)) %>% barplot(las = 2, main = "Library sizes")

dim(dds)

keep_genes <- rowSums(counts(dds)) > 0
dds <- dds[keep_genes, ]
dim(dds)

dds <- estimateSizeFactors(dds)
plot(sizeFactors(dds), colSums(counts(dds)),
     ylab = "Library sizes", xlab = "Size factors", cex = 1.5,
     main = "Size factors vs library size")

par(mfrow = c(1, 2))
boxplot(log2(counts(dds) + 1), notch = FALSE, las = 2,
        main = "Non-normalized read counts",
        ylab = "log2(read counts)", cex = 0.6)
boxplot(log2(counts(dds, normalized = TRUE) + 1), notch = FALSE, las = 2,
        main = "Size-factor-normalized read counts",
        ylab = "log2(read counts)", cex = 0.6)

par(mfrow = c(1, 1))

assay(dds, "log_counts")      <- log2(counts(dds, normalized = FALSE) + 1)
assay(dds, "log_norm_counts") <- log2(counts(dds, normalized = TRUE)  + 1)

par(mfrow = c(2, 2))
dds[, grepl("Unstimulated", dds$condition)] %>%
  assay("log_norm_counts") %>% .[, 1:2] %>%
  plot(cex = 0.1, main = "Unstimulated rep1 vs rep2")

dds[, grepl("PMA", dds$condition)] %>%
  assay("log_norm_counts") %>% .[, 1:2] %>%
  plot(cex = 0.1, main = "PMA rep1 vs rep2")

dds[, grepl("IL15", dds$condition)] %>%
  assay("log_norm_counts") %>% .[, 1:2] %>%
  plot(cex = 0.1, main = "IL15 rep1 vs rep2")

dds[, grepl("Bryostatin", dds$condition)] %>%
  assay("log_norm_counts") %>% .[, 1:2] %>%
  plot(cex = 0.1, main = "Bryostatin rep1 vs rep2")
par(mfrow = c(1, 1))

vsn::meanSdPlot(assay(dds, "log_norm_counts"), ranks = FALSE, plot = FALSE)$gg +
  labs(title = "Sequencing depth normalized log2(read counts)", y = "Standard deviation")

#rlog transformation
dst_rlog <- rlog(dds, blind = TRUE)
rlog_norm_counts <- assay(dst_rlog)

par(mfrow = c(2, 2))
dst_rlog[, grepl("Unstimulated", dst_rlog$condition)] %>%
  assay() %>% .[, 1:2] %>%
  plot(cex = 0.1, main = "Unstimulated rep1 vs rep2 (rlog)")
dst_rlog[, grepl("PMA", dst_rlog$condition)] %>%
  assay() %>% .[, 1:2] %>%
  plot(cex = 0.1, main = "PMA rep1 vs rep2 (rlog)")
dst_rlog[, grepl("IL15", dst_rlog$condition)] %>%
  assay() %>% .[, 1:2] %>%
  plot(cex = 0.1, main = "IL15 rep1 vs rep2 (rlog)")
dst_rlog[, grepl("Bryostatin", dst_rlog$condition)] %>%
  assay() %>% .[, 1:2] %>%
  plot(cex = 0.1, main = "Bryostatin rep1 vs rep2 (rlog)")
par(mfrow = c(1, 1))

vsn::meanSdPlot(assay(dst_rlog), ranks = FALSE, plot = FALSE)$gg +
  labs(title = "Following rlog transformation", x = "Mean", y = "Standard deviation") +
  coord_cartesian(ylim = c(0, 6))

corr_coeff <- cor(rlog_norm_counts, method ="pearson")
as.dist(1-corr_coeff, upper=TRUE) %>%
  as.matrix %>%
  pheatmap::pheatmap(main="Pearson correlation",
                     treeheight_row=0) 

par(mfrow=c(1,1))
as.dist(1-corr_coeff) %>%
  hclust %>%
  plot(labels=colnames(.), main="rlog transformed read counts")

library(tidyverse)
library(magrittr)
library(RColorBrewer)
dds %<>% DESeq()
dds %<>% estimateSizeFactors()
dds %<>% estimateDispersions()
dds %<>% nbinomWaldTest()
dds
rowData(dds) %>% colnames

rv <- rowVars(rlog_norm_counts)
top_variable <- order(rv, decreasing=TRUE)[seq_len(500)]
pca <- prcomp(t(rlog_norm_counts[top_variable, ]))

pca_var <- pca$sdev^2
pca_var_pct <- round(pca_var / sum(pca_var) * 100, 1)

barplot(pca_var_pct[1:10],
        names.arg = paste0("PC", 1:10),
        xlab = "Principal Component",
        ylab = "Variance Explained (%)",
        main = "Scree Plot",
        col = brewer.pal(10, "Set3"),
        border = NA)

plotPCA(dst_rlog) +
  labs(color=NULL) +
  theme_bw()

rowData(dds) %>% colnames

par(mfrow = c(1, 3))
rowData(dds)$WaldPvalue_condition_PMA_vs_Unstimulated %>%
  hist(breaks = 19, main = "PMA vs Unstimulated", xlab = "p-value")
rowData(dds)$WaldPvalue_condition_IL15_vs_Unstimulated %>%
  hist(breaks = 19, main = "IL15 vs Unstimulated", xlab = "p-value")
rowData(dds)$WaldPvalue_condition_Bryostatin_vs_Unstimulated %>%
  hist(breaks = 19, main = "Bryostatin vs Unstimulated", xlab = "p-value")

library(DESeq2)
library(ggplot2)
library(dplyr)
library(gridExtra)

DGE_results_PMA  <- results(dds,
                            contrast             = c("condition", "PMA", "Unstimulated"),
                            independentFiltering = TRUE,
                            alpha                = 0.05)
DGE_results_IL15 <- results(dds,
                            contrast             = c("condition", "IL15", "Unstimulated"),
                            independentFiltering = TRUE,
                            alpha                = 0.05)
DGE_results_bryo <- results(dds,
                            contrast             = c("condition", "Bryostatin", "Unstimulated"),
                            independentFiltering = TRUE,
                            alpha                = 0.05)

plot_filter <- function(res, title) {
  filt <- metadata(res)$filterNumRej
  filt$fitted_line <- metadata(res)$lo.fit$y
  theta_cut <- metadata(res)$filterTheta
  ggplot(filt) +
    geom_point(aes(x = theta, y = numRej)) +
    geom_line(aes(x = theta, y = fitted_line), color = "red") +
    geom_vline(xintercept = theta_cut) +
    labs(title = title,
         x = "Quantile threshold (theta)",
         y = "Number of rejections") +
    theme_bw()
}

filter_plot_PMA  <- plot_filter(DGE_results_PMA,  "PMA vs Unstimulated")
filter_plot_IL15 <- plot_filter(DGE_results_IL15, "IL15 vs Unstimulated")
filter_plot_bryo <- plot_filter(DGE_results_bryo, "Bryostatin vs Unstimulated")

grid.arrange(filter_plot_PMA, filter_plot_IL15, filter_plot_bryo, ncol = 3)

par(mfrow = c(1, 3))
hist(DGE_results_PMA$padj,
     breaks = 19,
     main = "PMA vs Unstimulated (adjusted p-values)",
     xlab = "Adjusted p-value")
hist(DGE_results_IL15$padj,
     breaks = 19,
     main = "IL15 vs Unstimulated (adjusted p-values)",
     xlab = "Adjusted p-value")
hist(DGE_results_bryo$padj,
     breaks = 19,
     main = "Bryostatin vs Unstimulated (adjusted p-values)",
     xlab = "Adjusted p-value")
par(mfrow = c(1, 1))

library(dplyr)
library(tibble)
library(DESeq2)

df_PMA <- as.data.frame(DGE_results_PMA) %>%
  tibble::rownames_to_column("gene") %>%
  filter(!is.na(padj))

top_gene_PMA <- df_PMA$gene[which.min(df_PMA$padj)]

par(mfrow = c(1, 2))
plotCounts(dds,
           gene = top_gene_PMA,
           normalized = TRUE,
           xlab = "",
           main = paste0("PMA: most significant (", top_gene_PMA, ")"))
plotCounts(dds,
           gene = df_PMA$gene[which.max(df_PMA$padj)],
           normalized = TRUE,
           xlab = "",
           main = "PMA: least significant gene")

df_IL15 <- as.data.frame(DGE_results_IL15) %>%
  tibble::rownames_to_column("gene") %>%
  filter(!is.na(padj))

top_gene_IL15 <- df_IL15$gene[which.min(df_IL15$padj)]

par(mfrow = c(1, 2))
plotCounts(dds,
           gene = top_gene_IL15,
           normalized = TRUE,
           xlab = "",
           main = paste0("IL15: most significant (", top_gene_IL15, ")"))
plotCounts(dds,
           gene = df_IL15$gene[which.max(df_IL15$padj)],
           normalized = TRUE,
           xlab = "",
           main = "IL15: least significant gene")

df_bryo <- as.data.frame(DGE_results_bryo) %>%
  tibble::rownames_to_column("gene") %>%
  filter(!is.na(padj))

top_gene_bryo <- df_bryo$gene[which.min(df_bryo$padj)]

par(mfrow = c(1, 2))
plotCounts(dds,
           gene = top_gene_bryo,
           normalized = TRUE,
           xlab = "",
           main = paste0("Bryostatin: most significant (", top_gene_bryo, ")"))
plotCounts(dds,
           gene = df_bryo$gene[which.max(df_bryo$padj)],
           normalized = TRUE,
           xlab = "",
           main = "Bryostatin: least significant gene")

par(mfrow = c(1, 1))

par(mfrow = c(1, 3))

plotCounts(dds,
           gene       = top_gene_PMA,
           normalized = TRUE,
           xlab       = "",
           main       = paste0("PMA: ", top_gene_PMA))

plotCounts(dds,
           gene       = top_gene_IL15,
           normalized = TRUE,
           xlab       = "",
           main       = paste0("IL15: ", top_gene_IL15))

plotCounts(dds,
           gene       = top_gene_bryo,
           normalized = TRUE,
           xlab       = "",
           main       = paste0("Bryostatin: ", top_gene_bryo))

par(mfrow = c(1, 1))

library(EnhancedVolcano)
library(patchwork)
library(RColorBrewer)

volcano_cols <- c("grey70",
                  brewer.pal(3, "Set2")[1],
                  brewer.pal(3, "Set2")[3],
                  brewer.pal(3, "Set1")[1])

df_PMA_shrunk  <- lfcShrink(dds, coef = "condition_PMA_vs_Unstimulated",        type = "apeglm")
df_IL15_shrunk <- lfcShrink(dds, coef = "condition_IL15_vs_Unstimulated",       type = "apeglm")
df_bryo_shrunk <- lfcShrink(dds, coef = "condition_Bryostatin_vs_Unstimulated", type = "apeglm")

vp_PMA_shrunk <- EnhancedVolcano(as.data.frame(df_PMA_shrunk),
                                 lab      = rownames(df_PMA_shrunk),
                                 x        = "log2FoldChange", y = "padj",
                                 pCutoff  = 0.05, FCcutoff = 1,
                                 title    = "PMA vs Unstimulated (shrunk)",
                                 subtitle = NULL, col = volcano_cols)

vp_IL15_shrunk <- EnhancedVolcano(as.data.frame(df_IL15_shrunk),
                                  lab      = rownames(df_IL15_shrunk),
                                  x        = "log2FoldChange", y = "padj",
                                  pCutoff  = 0.05, FCcutoff = 1,
                                  title    = "IL15 vs Unstimulated (shrunk)",
                                  subtitle = NULL, col = volcano_cols)

vp_bryo_shrunk <- EnhancedVolcano(as.data.frame(df_bryo_shrunk),
                                  lab      = rownames(df_bryo_shrunk),
                                  x        = "log2FoldChange", y = "padj",
                                  pCutoff  = 0.05, FCcutoff = 1,
                                  title    = "Bryostatin vs Unstimulated (shrunk)",
                                  subtitle = NULL, col = volcano_cols)

vp_PMA_shrunk + vp_IL15_shrunk + vp_bryo_shrunk +
  plot_layout(ncol = 3, guides = "collect") &
  theme(legend.position = "bottom")

par(mfrow = c(3, 2))

plotMA(DGE_results_PMA,  alpha = 0.05, main = "PMA (no shrinkage)",          ylim = c(-4, 4))
plotMA(df_PMA_shrunk,    alpha = 0.05, main = "PMA (with shrinkage)",        ylim = c(-4, 4))

plotMA(DGE_results_IL15, alpha = 0.05, main = "IL15 (no shrinkage)",         ylim = c(-4, 4))
plotMA(df_IL15_shrunk,   alpha = 0.05, main = "IL15 (with shrinkage)",       ylim = c(-4, 4))

plotMA(DGE_results_bryo, alpha = 0.05, main = "Bryostatin (no shrinkage)",   ylim = c(-4, 4))
plotMA(df_bryo_shrunk,   alpha = 0.05, main = "Bryostatin (with shrinkage)", ylim = c(-4, 4))

par(mfrow = c(1, 1))

library(org.Hs.eg.db)

hg38 <- org.Hs.eg.db
keytypes(hg38)
columns(hg38)

head(rownames(dds))

AnnotationDbi::select(
  hg38,
  keys = rownames(dds),
  keytype = "SYMBOL",
  columns = c("ENSEMBL", "ENTREZID", "GENENAME")
)

par(mfrow = c(1, 2))
plotCounts(dds, gene = "STAT5A", intgroup = "condition",
           normalized = TRUE, xlab = "", main = "STAT5A normalised counts")
plotCounts(dds, gene = "STAT5B", intgroup = "condition",
           normalized = TRUE, xlab = "", main = "STAT5B normalised counts")

par(mfrow = c(1, 1))

stat5_id <- rownames(dds)[grep("Stat5|STAT5", rownames(dds), ignore.case = TRUE)]
stat5_id

rlog_norm_counts[c("STAT5A", "STAT5B"), ] %>%
  as.data.frame() %>%
  tibble::rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "sample", values_to = "count") %>%
  mutate(condition = dds$condition[match(sample, colnames(rlog_norm_counts))]) %>%
  ggplot(aes(x = condition, y = count, fill = condition)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 2) +
  scale_fill_brewer(palette = "Set2") +
  facet_wrap(~ gene) +
  labs(
    title = "STAT5A and STAT5B expression across conditions",
    x     = "Condition",
    y     = "rlog-normalised expression"
  ) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

stat5_data <- rlog_norm_counts[c("STAT5A", "STAT5B"), ] %>%
  as.data.frame() %>%
  tibble::rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "sample", values_to = "count") %>%
  mutate(condition = dds$condition[match(sample, colnames(rlog_norm_counts))])

for (g in c("STAT5A", "STAT5B")) {
  cat("\n====", g, "====\n")
  df <- stat5_data %>% filter(gene == g)
  
  print(kruskal.test(count ~ condition, data = df))
  
  print(pairwise.wilcox.test(df$count, df$condition, 
                             p.adjust.method = "BH"))
}

df_annots_dge <- AnnotationDbi::select(
  hg38,
  keys    = rownames(df_PMA_shrunk),
  keytype = "SYMBOL",
  columns = c("GENENAME", "ENTREZID")
)

df_annots_dge <- df_annots_dge[!duplicated(df_annots_dge$SYMBOL), ]
df_annots_dge$gene_name <- ifelse(is.na(df_annots_dge$GENENAME),
                                  df_annots_dge$SYMBOL,
                                  df_annots_dge$GENENAME)
rownames(df_annots_dge) <- df_annots_dge$SYMBOL

df_PMA_annot_final <- merge(
  df_annots_dge[, c("SYMBOL", "gene_name", "ENTREZID")],
  as.data.frame(df_PMA_shrunk),
  by.x = "SYMBOL", by.y = "row.names", all.y = TRUE)

df_IL15_annot_final <- merge(
  df_annots_dge[, c("SYMBOL", "gene_name", "ENTREZID")],
  as.data.frame(df_IL15_shrunk),
  by.x = "SYMBOL", by.y = "row.names", all.y = TRUE)

df_bryo_annot_final <- merge(
  df_annots_dge[, c("SYMBOL", "gene_name", "ENTREZID")],
  as.data.frame(df_bryo_shrunk),
  by.x = "SYMBOL", by.y = "row.names", all.y = TRUE)

df_PMA_annot_final  %>% arrange(padj) %>% head(10)
df_IL15_annot_final %>% arrange(padj) %>% head(10)
df_bryo_annot_final %>% arrange(padj) %>% head(10)

write.table(subset(df_PMA_annot_final,  padj < 0.05),
            file = "DESeq2results_PMA_vs_Unstimulated.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(subset(df_IL15_annot_final, padj < 0.05),
            file = "DESeq2results_IL15_vs_Unstimulated.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(subset(df_bryo_annot_final, padj < 0.05),
            file = "DESeq2results_Bryostatin_vs_Unstimulated.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

library(DESeq2)
library(magrittr)

DGE_genes_PMA  <- as.data.frame(DGE_results_PMA)  %>% filter(!is.na(padj) & padj < 0.05) %>% arrange(padj)
DGE_genes_IL15 <- as.data.frame(DGE_results_IL15) %>% filter(!is.na(padj) & padj < 0.05) %>% arrange(padj)
DGE_genes_bryo <- as.data.frame(DGE_results_bryo) %>% filter(!is.na(padj) & padj < 0.05) %>% arrange(padj)

nrow(DGE_genes_PMA)
nrow(DGE_genes_IL15)
nrow(DGE_genes_bryo)

library(pheatmap)

deg_PMA  <- rownames(DGE_genes_PMA)
deg_IL15 <- rownames(DGE_genes_IL15)
deg_bryo <- rownames(DGE_genes_bryo)

genes_dge <- unique(c(deg_PMA, deg_IL15, deg_bryo))

rlog_dge <- rlog_norm_counts[genes_dge, ]

anno_col <- data.frame(
  condition = dds$condition,
  row.names = colnames(rlog_norm_counts)
)

pheatmap(rlog_dge,
         scale          = "none",
         show_rownames  = FALSE,
         annotation_col = anno_col,
         main           = "Combined DEG heatmap (rlog counts, unscaled)",
         color          = colorRampPalette(RColorBrewer::brewer.pal(7, "Greens"))(100))

pheatmap(rlog_dge,
         scale          = "row",
         show_rownames  = FALSE,
         annotation_col = anno_col,
         main           = "Combined DEG heatmap (row-scaled z-score)")

lfc_matrix <- data.frame(
  PMA  = as.data.frame(DGE_results_PMA)[genes_dge,  "log2FoldChange"],
  IL15 = as.data.frame(DGE_results_IL15)[genes_dge, "log2FoldChange"],
  Bryo = as.data.frame(DGE_results_bryo)[genes_dge, "log2FoldChange"],
  row.names = genes_dge
)

top50_genes <- lfc_matrix %>%
  mutate(mean_abs_lfc = rowMeans(abs(across(everything())))) %>%
  arrange(desc(mean_abs_lfc)) %>%
  head(50) %>%
  rownames()

mat_scaled <- t(scale(t(rlog_dge[top50_genes, ])))

pheatmap(mat_scaled,
         annotation_col = anno_col,
         cluster_rows   = TRUE,
         cluster_cols   = TRUE,
         show_rownames  = TRUE,
         show_colnames  = TRUE,
         fontsize_row   = 7,
         border_color   = NA,
         breaks         = seq(-2, 2, length.out = 101),
         main           = "Top 50 DE Genes (padj < 0.05, z-scored rlog counts)")

ptefb_genes <- c("CDK9", "CCNT1", "CCNT2")

DGE_genes_PMA  %>% filter(rownames(.) %in% ptefb_genes)
DGE_genes_IL15 %>% filter(rownames(.) %in% ptefb_genes)
DGE_genes_bryo %>% filter(rownames(.) %in% ptefb_genes)

rlog_norm_counts[ptefb_genes[ptefb_genes %in% rownames(dds)], ] %>%
  as.data.frame() %>%
  tibble::rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "sample", values_to = "count") %>%
  mutate(condition = dds$condition[match(sample, colnames(rlog_norm_counts))]) %>%
  ggplot(aes(x = condition, y = count, fill = condition)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 2) +
  scale_fill_brewer(palette = "Set2") +
  facet_wrap(~ gene, scales = "free_y") +
  labs(title = "P-TEFb complex gene expression across conditions",
       x = "Condition", y = "rlog-normalised expression") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(DOSE)

organism <- org.Hs.eg.db
ego_PMA <- enrichGO(
  gene          = rownames(DGE_genes_PMA),
  universe      = rownames(dds),
  ont           = "ALL",
  keyType       = "SYMBOL",
  minGSSize     = 3,
  maxGSSize     = 800,
  pvalueCutoff  = 0.05,
  OrgDb         = organism,
  pAdjustMethod = "BH"
)

ego_IL15 <- enrichGO(
  gene          = rownames(DGE_genes_IL15),
  universe      = rownames(dds),
  ont           = "ALL",
  keyType       = "SYMBOL",
  minGSSize     = 3,
  maxGSSize     = 800,
  pvalueCutoff  = 0.05,
  OrgDb         = organism,
  pAdjustMethod = "BH"
)

ego_bryo <- enrichGO(
  gene          = rownames(DGE_genes_bryo),
  universe      = rownames(dds),
  ont           = "ALL",
  keyType       = "SYMBOL",
  minGSSize     = 3,
  maxGSSize     = 800,
  pvalueCutoff  = 0.05,
  OrgDb         = organism,
  pAdjustMethod = "BH"
)

ego_PMA[1, ] %>% str
ego_IL15[1, ] %>% str
ego_bryo[1, ] %>% str

rnp_biogenesis_genes <- unlist(strsplit(ego_PMA[1, "geneID"], "/"))
head(rnp_biogenesis_genes)

mitotic_cycle_genes <- unlist(strsplit(ego_IL15[1, "geneID"], "/"))
head(mitotic_cycle_genes)

chromosome_org_genes <- unlist(strsplit(ego_bryo[1, "geneID"], "/"))
head(mitotic_cycle_genes)

write.table(ego_PMA@result[ , c("ID", "pvalue")],
            file="enrichGO_PMA.txt", sep="\t",
            quote=FALSE, row.names=FALSE)
write.table(ego_IL15@result[ , c("ID", "pvalue")],
            file="enrichGO_IL15.txt", sep="\t",
            quote=FALSE, row.names=FALSE)
write.table(ego_bryo@result[ , c("ID", "pvalue")],
            file="enrichGO_Bryo.txt", sep="\t",
            quote=FALSE, row.names=FALSE)

gene_list_PMA  <- sort(setNames(DGE_results_PMA$log2FoldChange,  rownames(DGE_results_PMA)),  decreasing = TRUE)

gene_list_IL15 <- sort(setNames(DGE_results_IL15$log2FoldChange, rownames(DGE_results_IL15)), decreasing = TRUE)

gene_list_bryo <- sort(setNames(DGE_results_bryo$log2FoldChange, rownames(DGE_results_bryo)), decreasing = TRUE)

gene_list_PMA  <- gene_list_PMA[!is.na(gene_list_PMA)]
gene_list_IL15 <- gene_list_IL15[!is.na(gene_list_IL15)]
gene_list_bryo <- gene_list_bryo[!is.na(gene_list_bryo)]

gse_PMA <- gseGO(
  geneList      = gene_list_PMA,
  ont           = "ALL",
  keyType       = "SYMBOL",      
  minGSSize     = 3,
  maxGSSize     = 800,
  pvalueCutoff  = 0.05,
  verbose       = TRUE,
  OrgDb         = organism,
  pAdjustMethod = "BH"
)

gse_IL15 <- gseGO(
  geneList      = gene_list_IL15,
  ont           = "ALL",
  keyType       = "SYMBOL",
  minGSSize     = 3,
  maxGSSize     = 800,
  pvalueCutoff  = 0.05,
  verbose       = TRUE,
  OrgDb         = organism,
  pAdjustMethod = "BH"
)

gse_bryo <- gseGO(
  geneList      = gene_list_bryo,
  ont           = "ALL",
  keyType       = "SYMBOL",
  minGSSize     = 3,
  maxGSSize     = 800,
  pvalueCutoff  = 0.05,
  verbose       = TRUE,
  OrgDb         = organism,
  pAdjustMethod = "BH"
)

dotplot(gse_PMA,  showCategory = 10, split = ".sign", title = "GSEA PMA vs Unstimulated") +
  facet_grid(. ~ .sign)

dotplot(gse_IL15, showCategory = 10, split = ".sign", title = "GSEA IL15 vs Unstimulated") +
  facet_grid(. ~ .sign)

dotplot(gse_bryo, showCategory = 10, split = ".sign", title = "GSEA Bryostatin vs Unstimulated") +
  facet_grid(. ~ .sign)

gse_PMA$Description[1]
gse_IL15$Description[1]
gse_bryo$Description[1]

gseaplot(gse_PMA,  by = "all", title = gse_PMA$Description[1],  geneSetID = 1)
gseaplot(gse_IL15, by = "all", title = gse_IL15$Description[1], geneSetID = 1)
gseaplot(gse_bryo, by = "all", title = gse_bryo$Description[1], geneSetID = 1)