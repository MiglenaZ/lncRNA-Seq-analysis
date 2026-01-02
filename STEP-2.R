library(data.table)
library(DESeq2)
library(data.table)
library(dplyr)
library(org.Hs.eg.db)
library(AnnotationDbi)

#----------------------------------------------------------------
# create a count matrix - merged_featureCounts.txt

setwd("<OUTPUT_DIR_OF_LNCPIPE.SBATCH>")
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

setwd("<OUTPUT_DIR_OF_LNCPIPE.SBATCH>")
cat("=== Step 1: Load Data ===\n")

# load count matrix and metadata
count_matrix <- fread("merged_featureCounts.txt")
# create sample_metadata.txt
sample_info <- fread("sample_metadata.txt")
cat("Samples:", nrow(sample_info), "\n")
cat("Genes:", nrow(count_matrix), "\n")

# Prepare counts
counts_only <- as.matrix(count_matrix[, 7:ncol(count_matrix)])
rownames(counts_only) <- count_matrix$Geneid

# Filter genes (keep if detected in at least 3 samples)
keep <- rowSums(counts_only > 0) >= 3
counts_filtered <- counts_only[keep, ]
cat("Genes after filtering:", nrow(counts_filtered), "\n")

# Create DESeq2 object
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

# Convert to data frame
res_df <- as.data.frame(res) %>%
  rownames_to_column(var = "ENSEMBL")

# Clean ENSEMBL IDs
res_df$ENSEMBL_clean <- sub("\\..*", "", res_df$ENSEMBL)

# Add gene symbols and ENTREZ IDs
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
# Extract gene types from GENCODE GTF
gtf_file <- "<FULL_PATH_TO_GTF_FILE>"
gtf <- fread(cmd = paste0("zcat ", gtf_file, " | grep -v '^#' | awk '$3==\"gene\"'"),
             sep = "\t", header = FALSE)
colnames(gtf) <- c("chr", "source", "type", "start", "end", "score", "strand", "frame", "attributes")

# Parse gene_id and gene_type
parse_attr <- function(attr_str) {
  gene_id <- sub(".*gene_id \"([^\"]+)\".*", "\\1", attr_str)
  gene_type <- sub(".*gene_type \"([^\"]+)\".*", "\\1", attr_str)
  data.frame(gene_id = gene_id, gene_type = gene_type, stringsAsFactors = FALSE)
}
parsed <- do.call(rbind, lapply(gtf$attributes, parse_attr))
parsed$gene_id_clean <- sub("\\..*", "", parsed$gene_id)

# Merge with results
res_annotated <- merge(res_df, parsed[, c("gene_id_clean", "gene_type")],
                       by.x = "ENSEMBL_clean",
                       by.y = "gene_id_clean",
                       all.x = TRUE)
cat("Gene type distribution:\n")
print(head(sort(table(res_annotated$gene_type, useNA = "always"), decreasing = TRUE), 10))

# Save final results
fwrite(res_annotated, "DESeq2_annotated.csv")

# Extract significant lncRNAs
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
