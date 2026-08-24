#!/usr/bin/env Rscript
# DESeq2 analysis: normalize, estimate dispersions,
# fit model, and extract results for each contrast.

library(DESeq2)
library(dplyr)
library(readr)
library(purrr)

# Capture log
logfile <- snakemake@log[[1]]
dir.create(dirname(logfile), showWarnings = FALSE, recursive = TRUE)
logfile <- file(logfile, open = "w")
sink(logfile, type = "output")
sink(logfile, type = "message", append = TRUE)

cat("=== DESeq2 Analysis ===\n")
cat("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    "\n\n")

# Read count matrix
cat("Reading count matrix...\n")
counts_file <- snakemake@input[["counts"]]
counts_df <- read_tsv(counts_file)

cat("Input file:", counts_file, "\n")
cat("Dimensions:", nrow(counts_df), "genes x",
    ncol(counts_df) - 2, "samples\n")

# Extract gene IDs and count matrix
gene_ids <- counts_df[[1]]
counts_matrix <- as.matrix(counts_df[, -c(1, 2)])
rownames(counts_matrix) <- gene_ids

cat("\nCount matrix summary:\n")
cat("- Min:", min(counts_matrix), "\n")
cat("- Max:", max(counts_matrix), "\n")
cat("- Mean:", round(mean(counts_matrix), 2), "\n")

# Build sample metadata
samples_info <- snakemake@params[["samples_info"]]
# TODO: Add input/config checks.
counts_matrix <- counts_matrix[, names(samples_info)]

cat("\nSample mapping:\n")
# TODO: Add input/config checks.
samples_info <- strsplit(unlist(samples_info), '_', fixed=TRUE)
coldata <- data.frame(sample = colnames(counts_matrix))
coldata$background <- sapply(samples_info, "[[", 1)[coldata$sample]
coldata$condition <- sapply(samples_info, "[[", 2)[coldata$sample]
# TODO: Cleanup.
#print(coldata)
#foo
#
#for (i in seq_len(nrow(coldata))) {
#  sample_name <- coldata$sample[i]
#  if (sample_name %in% names(samples_info)) {
#    condition_label <- samples_info[[sample_name]]
#    parts <- strsplit(condition_label, "_")[[1]]
#    if (length(parts) == 2) {
#      coldata$background[i] <- parts[1]
#      coldata$condition[i] <- parts[2]
#    }
#  }
#}

print(coldata)

# Validate metadata
if (anyNA(coldata$condition) || anyNA(coldata$background)) {
  stop("ERROR: Could not map all samples to conditions/backgrounds.
       Check config.yaml 'samples' section.")
}

# Create DESeqDataSet
cat("\nCreating DESeqDataSet...\n")
dds <- DESeqDataSetFromMatrix(
  countData = round(counts_matrix),
  colData = coldata,
  design = ~ condition
)

cat("DESeqDataSet created with", nrow(dds), "genes and",
    ncol(dds), "samples\n")

# Pre-filtering: remove genes with very low counts
cat("\nPre-filtering low-count genes...\n")
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]
cat("Genes retained after filtering:", nrow(dds), "\n")

# Run DESeq2
cat("\nRunning DESeq2 (normalize, estimate dispersions, fit)...\n")
dds <- DESeq(dds)

cat("✓ Normalization complete\n")
cat("✓ Dispersion estimation complete\n")
cat("✓ Model fitting complete\n")

# Save DDS object
cat("\nSaving DESeqDataSet object...\n")
output_dds <- snakemake@output[["dds_rds"]]
saveRDS(dds, output_dds)
cat("Saved:", output_dds, "\n")

# Extract results for each contrast
cat("\n=== Extracting Results for Each Contrast ===\n")
contrasts <- snakemake@params[["contrasts"]]

for (contrast_name in names(contrasts)) {
  cat("\n---", contrast_name, "---\n")
  contrast <- contrasts[[contrast_name]]
  cat("Background:", contrast$background, "\n")
  cat("Condition1 (numerator):", contrast$condition1, "\n")
  cat("Condition2 (denominator):", contrast$condition2, "\n")

  # Filter samples by background
  samples_for_contrast <- which(
    coldata$background == contrast$background
  )

  if (length(samples_for_contrast) == 0) {
    warning("No samples found for background: ",
            contrast$background)
    next
  }

  # Check replicate count
  n_samples <- length(samples_for_contrast)
  if (n_samples < 2) {
    warning("WARNING: Only ", n_samples,
            " sample(s) for contrast ",
            contrast$name,
            ". DESeq2 requires n >= 2.")
  }

  # Extract results
  tryCatch({
    res <- results(
      dds,
      contrast = c("condition", contrast$condition1, contrast$condition2),
      independentFiltering = TRUE,
      alpha = snakemake@params[["padj_threshold"]]
    )

    # Convert to data frame
    res_df <- as.data.frame(res)
    res_df$gene_id <- rownames(res)
    res_df <- res_df[, c("gene_id", "log2FoldChange",
                         "pvalue", "padj")]

    # Sort by padj
    res_df <- res_df[order(res_df$padj), ]

    # Count significant genes
    n_sig <- sum(res_df$padj < snakemake@params[[
      "padj_threshold"]], na.rm = TRUE)
    cat("Significant genes (padj <",
        snakemake@params[["padj_threshold"]], "):",
        n_sig, "\n")

    # Save results
    output_file <- snakemake@output[["results"]][[
      match(contrast_name, names(contrasts))
    ]]
    write_csv(res_df, output_file)
    cat("Saved:", output_file, "\n")

  }, error = function(e) {
    warning("Error extracting contrast ", contrast$name,
            ": ", e$message)
  })
}

cat("\n=== DESeq2 Analysis Complete ===\n")
sink(type="message")
sink()
close(logfile)
