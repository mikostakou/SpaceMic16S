# 16S phyloseq object construction and rarefaction
#
# This script prepares the 16S phyloseq objects used for downstream analyses.
#
# Main steps:
# - Load ASV and taxonomy tables from the DADA2 pipeline
# - Retain bacterial ASVs and remove chloroplast and mitochondrial sequences
# - Combine ASV table, taxonomy table, and metadata into a phyloseq object
# - Summarize sequencing depth per sample
# - Rarefy samples to 8,000 reads
# - Save unrarefied and rarefied phyloseq objects


library(phyloseq)

sample_depth <- 8000
seed <- 1

dada2_path <- "results/dada2_16S"
out_path <- "results/phyloseq_16S"

dir.create(out_path, showWarnings = FALSE, recursive = TRUE)


# Load DADA2 outputs

seqtab <- readRDS(file.path(dada2_path, "seqtab16S.nochim.rds"))
taxa <- readRDS(file.path(dada2_path, "taxtable16S.rds"))


# Keep bacterial ASVs and remove organelles

tax <- as.data.frame(taxa, stringsAsFactors = FALSE)

tax_string <- apply(
  tax, 1,
  function(x) paste(x[!is.na(x)], collapse = ";")
)

chloroplast <- grepl("chloroplast", tax_string, ignore.case = TRUE)
mitochondria <- grepl("mitochond", tax_string, ignore.case = TRUE)

keep <- !is.na(tax$Kingdom) &
  tax$Kingdom == "Bacteria" &
  !chloroplast &
  !mitochondria

seqtab <- seqtab[, keep, drop = FALSE]
taxa <- taxa[colnames(seqtab), , drop = FALSE]


# Load metadata

meta <- read.table(
  "data/meta_16S_Bexis.txt",
  header = TRUE,
  row.names = 1,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

meta[] <- lapply(
  meta,
  function(x) if (is.character(x)) trimws(x) else x
)



# Keep samples shared between ASV table and metadata

common_samples <- intersect(rownames(meta), rownames(seqtab))

seqtab <- seqtab[common_samples, , drop = FALSE]
meta <- meta[common_samples, , drop = FALSE]


# Remove ASVs absent from the matched samples

present_asvs <- colSums(seqtab) > 0

seqtab <- seqtab[, present_asvs, drop = FALSE]
taxa <- taxa[colnames(seqtab), , drop = FALSE]


# Create unrarefied phyloseq object

physeq16S <- phyloseq(
  otu_table(seqtab, taxa_are_rows = FALSE),
  sample_data(meta),
  tax_table(as.matrix(taxa))
)

saveRDS(
  physeq16S,
  file.path(out_path, "physeq16S_bacteria_unrarefied.rds")
)


# Summarize sequencing depth

read_counts <- sort(sample_sums(physeq16S))

write.table(
  read_counts,
  file.path(out_path, "sample_read_counts_sorted.txt"),
  sep = "\t",
  quote = FALSE,
  col.names = FALSE
)


# Rarefy to 8,000 reads

physeq_8000 <- prune_samples(
  sample_sums(physeq16S) >= sample_depth,
  physeq16S
)

set.seed(seed)

rarefied16S <- rarefy_even_depth(
  physeq_8000,
  sample.size = sample_depth,
  rngseed = seed,
  replace = FALSE
)

rarefied16S <- prune_taxa(
  taxa_sums(rarefied16S) > 0,
  rarefied16S
)


# Save rarefied phyloseq object

saveRDS(
  rarefied16S,
  file.path(out_path, "physeq16S_rarefied_8000.rds")
)
