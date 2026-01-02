#!/usr/bin/env Rscript
library(data.table)
library(DESeq2)
library(dplyr)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(ggplot2)
library(dplyr)
library(tibble)
library(ggrepel)
library(RColorBrewer)
library(clusterProfiler)
library(enrichplot)
library(DOSE)


setwd("<OUTPUT_DIR_OF_LNCPIPE.SBATCH>")
#----------------------------------------------------------------
# create a count matrix - merged_featureCounts.txt

# find all featureCounts files
count_files <- list.files("subread", pattern = "*.featureCounts.tsv$", full.names = TRUE)
cat("Found", length(count_files), "count files\n")

# read first file
first_file <- fread(count_files[1], skip = 1)
count_matrix <- first_file[, 1:6]  # Keep annotation columns

# add counts from all samples
for (file in count_files) {
  sample_name <- gsub(".featureCounts.tsv", "", basename(file))
  dt <- fread(file, skip = 1)
  count_matrix <- cbind(count_matrix, dt[[7]])
  colnames(count_matrix)[ncol(count_matrix)] <- sample_name
}

cat("Count matrix dimensions:", dim(count_matrix), "\n")
fwrite(count_matrix, "merged_featureCounts.txt", sep = "\t")
cat("✓ Created merged_featureCounts.txt\n")

#----------------------------------------------------------------
# DESeq2 analysis

cat("=== Step 1: Load Data ===\n")

# load count matrix and metadata
count_matrix <- fread("merged_featureCounts.txt")
# create sample_metadata.txt
sample_info <- fread("sample_metadata.txt")
cat("Samples:", nrow(sample_info), "\n")
cat("Genes:", nrow(count_matrix), "\n")

# prepare counts
counts_only <- as.matrix(count_matrix[, 7:ncol(count_matrix)])
rownames(counts_only) <- count_matrix$Geneid
colnames(counts_only) <- sample_info$sample

# filter genes (keep if detected in at least 3 samples)
keep <- rowSums(counts_only > 0) >= 3
counts_filtered <- counts_only[keep, ]
cat("Genes after filtering:", nrow(counts_filtered), "\n")

# create DESeq2 object
# define design (here it's 2 factors being accounted for: cell line and treatment groups (C for Control samples and T for treated samples))
# the current design(dds) function is made to answer: What genes change due to Treatment (T vs C), after adjusting for differences between cell lines?

dds <- DESeqDataSetFromMatrix(
  countData = counts_filtered,
  colData = sample_info,
  design = ~ CellLine + Treatment
)

cat("\n=== Step 2: Run DESeq2 ===\n")
dds <- DESeq(dds)
res <- results(dds, contrast = c("Treatment", "T", "C"))

# convert to data frame
res_df <- as.data.frame(res) %>%
  rownames_to_column(var = "ENSEMBL")

# clean ENSEMBL IDs
res_df$ENSEMBL_clean <- sub("\\..*", "", res_df$ENSEMBL)

# add gene symbols and ENTREZ IDs
res_df$SYMBOL <- mapIds(org.Hs.eg.db,
                        keys = res_df$ENSEMBL_clean,
                        column = "SYMBOL",
                        keytype = "ENSEMBL",
                        multiVals = "first")
res_df$ENTREZID <- mapIds(org.Hs.eg.db,
                          keys = res_df$ENSEMBL_clean,
                          column = "ENTREZID",
                          keytype = "ENSEMBL",
                          multiVals = "first")

cat("\n=== Step 3: Annotate with GENCODE Gene Types ===\n")

# extract gene types from GENCODE GTF
gtf_file <- "<FULL_PATH_TO_GTF_FILE>"
gtf <- fread(cmd = paste0("zcat ", gtf_file, " | grep -v '^#' | awk '$3==\"gene\"'"),
             sep = "\t", header = FALSE)
colnames(gtf) <- c("chr", "source", "type", "start", "end", "score", "strand", "frame", "attributes")

# parse gene_id and gene_type
parse_attr <- function(attr_str) {
  gene_id <- sub(".*gene_id \"([^\"]+)\".*", "\\1", attr_str)
  gene_type <- sub(".*gene_type \"([^\"]+)\".*", "\\1", attr_str)
  data.frame(gene_id = gene_id, gene_type = gene_type, stringsAsFactors = FALSE)
}
parsed <- do.call(rbind, lapply(gtf$attributes, parse_attr))
parsed$gene_id_clean <- sub("\\..*", "", parsed$gene_id)

# merge with results
res_annotated <- merge(res_df, parsed[, c("gene_id_clean", "gene_type")],
                       by.x = "ENSEMBL_clean",
                       by.y = "gene_id_clean",
                       all.x = TRUE)
cat("Gene type distribution:\n")
print(head(sort(table(res_annotated$gene_type, useNA = "always"), decreasing = TRUE), 10))

# save final results
fwrite(res_annotated, "DESeq2_annotated.csv")

# extract significant lncRNAs
sig_lncrna <- res_annotated %>%
  filter(gene_type == "lncRNA") %>%
  filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 1.0) %>%
  arrange(padj) %>%
  mutate(SYMBOL = ifelse(is.na(SYMBOL), ENSEMBL_clean, SYMBOL))  # Fix NA symbols

fwrite(sig_lncrna, "DESeq2_significant_lncRNAs.csv")

cat("\n=== RESULTS ===\n")
cat("Total genes analyzed:", nrow(res_annotated), "\n")
cat("Total lncRNAs:", sum(res_annotated$gene_type == "lncRNA", na.rm = TRUE), "\n")
cat("Significant lncRNAs:", nrow(sig_lncrna), "\n")
cat("\n✓ Analysis complete!\n")
cat("Files created:\n")
cat("  - DESeq2_FINAL_annotated.csv\n")
cat("  - DESeq2_significant_lncRNAs_FINAL.csv\n")

#----------------------------------------------------------------
# filtering for top significant RNA and lncRNA genes

# get TOP 5 significant lncRNAs
top5_lncrna <- sig_lncrna %>%
  arrange(padj) %>%
  head(5)

cat("\n=== TOP 5 Most Significant lncRNAs ===\n")
print(top5_lncrna[, c("ENSEMBL_clean", "SYMBOL", "baseMean", "log2FoldChange", "padj")])

top5_ids <- top5_lncrna$ENSEMBL_clean
top5_labels <- ifelse(!is.na(top5_lncrna$SYMBOL),
                     top5_lncrna$SYMBOL,
                     top5_lncrna$ENSEMBL_clean)

# get all significant genes for GO/KEGG
sig_all <- res_annotated %>%
  filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 1.0)

sig_protein <- sig_all %>%
  filter(gene_type == "protein_coding") %>%
  filter(!is.na(ENTREZID))

cat("\nSignificant genes for enrichment:\n")
cat("  Protein-coding:", nrow(sig_protein), "\n")
cat("  lncRNAs:", sum(sig_all$gene_type == "lncRNA"), "\n")

# -----------------------------------------------------
# VISUALIZATION 1: PCA Plot (TOP 5 labeled)

cat("\nStep 3: Creating PCA plot...\n")

vsd <- vst(dds, blind = FALSE)
pcaData <- plotPCA(vsd, intgroup = c("CellLine", "Treatment"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

pdf("PCA_plot.pdf", width = 10, height = 8)
ggplot(pcaData, aes(x = PC1, y = PC2, color = Treatment, shape = CellLine)) +
  geom_point(size = 5, alpha = 0.8) +
  scale_color_manual(values = c("C" = "#00BFC4", "T" = "#F8766D"),
                    labels = c("Control", "Treated")) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA Analysis: Sample Clustering by Treatment") +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
        legend.position = "right",
        legend.title = element_text(face = "bold"),
        panel.grid.major = element_line(color = "gray90"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1))
dev.off()

# -----------------------------------------------------
# VISUALIZATION 2: Volcano Plot (TOP 5 labeled)

cat("Step 4: Creating Volcano plot with TOP 5 lncRNAs...\n")

volcano_data <- res_annotated %>%
  filter(!is.na(padj), !is.na(log2FoldChange)) %>%
  mutate(
    minuslog10padj = -log10(padj),
    gene_category = case_when(
      ENSEMBL_clean %in% top5_ids ~ "TOP 5 lncRNA",
      gene_type == "lncRNA" & padj < 0.05 & abs(log2FoldChange) > 1.0 ~ "Other sig lncRNA",
      padj < 0.05 & abs(log2FoldChange) > 1.0 ~ "Significant",
      TRUE ~ "Not significant"
    ),
    label = ifelse(ENSEMBL_clean %in% top5_ids,
                  ifelse(!is.na(SYMBOL), SYMBOL, ENSEMBL_clean),
                  "")
  )

# Order for plotting (TOP 5 on top)
volcano_data$gene_category <- factor(volcano_data$gene_category,
                                     levels = c("Not significant", "Significant",
                                               "Other sig lncRNA", "TOP 5 lncRNA"))

pdf("Volcano_TOP5_lncRNA.pdf", width = 12, height = 10)
ggplot(volcano_data, aes(x = log2FoldChange, y = minuslog10padj)) +
    geom_point(aes(color = gene_category, size = gene_category, alpha = gene_category)) +
    scale_color_manual(values = c("Not significant" = "gray70",
                                  "Significant" = "#F9C107",
                                  "Other sig lncRNA" = "#FF6A00",
                                  "TOP 5 lncRNA" = "red")) +
    scale_size_manual(values = c("Not significant" = 0.5,
                                 "Significant" = 1,
                                 "Other sig lncRNA" = 2,
                                 "TOP 5 lncRNA" = 4)) +
    scale_alpha_manual(values = c("Not significant" = 0.3,
                                  "Significant" = 0.5,
                                  "Other sig lncRNA" = 0.7,
                                  "TOP 5 lncRNA" = 1)) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue", linewidth = 0.5) +
    geom_vline(xintercept = c(-1.0, 1.0), linetype = "dashed", color = "blue", linewidth = 0.5) +
    geom_text_repel(data = volcano_data %>% filter(label != ""),
                    aes(label = label),
                    size = 5,
                    fontface = "bold",
                    box.padding = 1,
                    point.padding = 0.5,
                    max.overlaps = 100,
                    min.segment.length = 0,
                    segment.color = "red",
                    segment.size = 0.5) +
    labs(title = "Volcano Plot: Differential Expression Analysis\n(TOP 5 lncRNAs Highlighted in Red)",
         subtitle = paste0("Total significant genes: ", nrow(sig_all),
                           " | Significant lncRNAs: ", sum(sig_all$gene_type == "lncRNA")),
         x = "Log2 Fold Change (Treatment vs Control)",
         y = "-Log10 Adjusted P-value",
         color = "Gene Category") +
    guides(size = "none", alpha = "none") +
    theme_bw(base_size = 14) +
    theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5, size = 12),
          legend.position = "right",
          legend.title = element_text(face = "bold"),
          panel.grid.major = element_line(color = "gray90"))
dev.off()

# -----------------------------------------------------
# VISUALIZATION 3: Log2 Fold Change Barplot (TOP 5 only)

cat("Step 5: Creating Log2FC barplot for TOP 5...\n")

top5_plot <- top5_lncrna %>%
  mutate(gene_label = ifelse(!is.na(SYMBOL), SYMBOL, ENSEMBL_clean),
         direction = ifelse(log2FoldChange > 0, "Up", "Down"))

pdf("TOP5_Log2FC_barplot.pdf", width = 10, height = 6)
ggplot(top5_plot, aes(x = reorder(gene_label, log2FoldChange),
                      y = log2FoldChange,
                      fill = direction)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(aes(label = sprintf("%.2f", log2FoldChange),
                y = log2FoldChange + sign(log2FoldChange) * 0.2),
            size = 5, fontface = "bold") +
  scale_fill_manual(values = c("Up" = "#F8766D", "Down" = "#00BFC4")) +
  coord_flip() +
  labs(title = "TOP 5 Most Significant lncRNAs",
       subtitle = "Log2 Fold Change (Treated vs Control)",
       x = "lncRNA",
       y = "Log2 Fold Change",
       fill = "Regulation") +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "right",
        axis.text.y = element_text(face = "bold", size = 12),
        panel.grid.major.y = element_blank())
dev.off()
