# Install BiocManager and ANCOMBC if needed
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("ANCOMBC")

# Load libraries
library(biomformat)
library(ANCOMBC)
library(phyloseq)

# Set logging to Snakemake log file
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

# Log versions
cat("R version:", R.version.string, "\n")
cat("ANCOMBC version:", as.character(packageVersion("ANCOMBC")), "\n")

# Load the input table
print("Loading table...")
table <- biomformat::read_biom(snakemake@input[["table"]])
table <- as.matrix(biomformat::biom_data(table))

# Load the metadata
print("Loading metadata...")
metadata <- read.table(snakemake@input[["metadata"]], sep = "\t", header = TRUE, row.names = 1)

# Get parameters from Snakemake
covariate <- snakemake@params[[1]][["factor_name"]]
target <- snakemake@params[[1]][["target_level"]]
reference <- snakemake@params[[1]][["reference_level"]]
confounders <- snakemake@params[[1]][["confounders"]]

# Harmonize table and metadata samples
print("Harmonizing table and metadata samples...")
samples <- colnames(table)
metadata <- subset(metadata, rownames(metadata) %in% samples)
metadata[[covariate]] <- as.factor(metadata[[covariate]])
metadata[[covariate]] <- relevel(metadata[[covariate]], reference)
sample_order <- rownames(metadata)
table <- table[, sample_order]
rownames(table) <- paste0("F_", rownames(table))  # Prevent feature renaming by R

# Create phyloseq object
print("Converting to phyloseq...")
taxa <- otu_table(table, taxa_are_rows = TRUE)
meta <- sample_data(metadata)
physeq <- phyloseq(taxa, meta)

# Run ANCOMBC2
print("Running ANCOMBC2...")
fix_formula <- if (length(confounders) != 0) {
  paste(c(covariate, confounders), collapse = " + ")
} else {
  covariate
}

print(fix_formula)

ancombc.results <- ancombc2(
  data = physeq,
  assay_name = "counts",
  tax_level = NULL,
  fix_formula = fix_formula,
  group = covariate,
  p_adj_method = "BH",
  struc_zero = TRUE,
  neg_lb = TRUE,
  alpha = 0.05,
  global = FALSE,
  lib_cut = 0,
  s0_perc = 0.05
)


print("Saved RDS!")
saveRDS(ancombc.results, snakemake@output[[2]])

# Extract results
print("Extracting results...")
res <- ancombc.results$res

# Dynamically construct the column names
coef_col <- paste0("lfc_", covariate, target)
pval_col <- paste0("p_", covariate, target)
qval_col <- paste0("q_", covariate, target)

# Extract columns from the result dataframe
coefs <- res[[coef_col]]
pvals <- res[[pval_col]]
qvals <- res[[qval_col]]

# Combine results
print("Combining results...")
results_all <- data.frame(
  taxon = gsub("^F_", "", res$taxon),
  coef = coefs,
  pval = pvals,
  qval = qvals
)

# Save results
write.table(results_all, file = snakemake@output[[1]], sep = "\t", quote = FALSE, row.names = FALSE)
print("Saved differentials!")
