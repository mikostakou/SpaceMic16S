# 16S phyloseq object construction and rarefaction
#
# This script prepares the 16S phyloseq object used for downstream analyses.
#
# Main steps:
# - Load ASV table and taxonomy table from the DADA2 pipeline
# - Combine ASV table, taxonomy table, and metadata into a phyloseq object
# - Export OTU abundance, taxonomy, and metadata tables
# - Summarize sequencing depth per sample
# - Rarefy the phyloseq object to the selected sequencing depth

library(phyloseq)
library(openxlsx)

# Rarefaction depth selected after inspecting sequencing depth distributions.
sample_depth <- 8000
seed <- 1

# Load DADA2 outputs

seqtab_nochim <- readRDS("seqtab16S_paired.nochim.rds")
taxa <- readRDS("taxtable16S_paired.rds")

# Load metadata

meta <- read.table(
  "meta16sf.txt",
  header = TRUE,
  row.names = 1,
  sep = "\t",
  stringsAsFactors = FALSE
)

# Keep only samples shared between ASV table and metadata

common_samples <- intersect(rownames(seqtab_nochim), rownames(meta))

seqtab_nochim <- seqtab_nochim[common_samples, , drop = FALSE]
meta <- meta[common_samples, , drop = FALSE]

# Prepare taxonomy table

taxa <- tax_table(taxa)
taxa_names(taxa) <- colnames(seqtab_nochim)

# Create phyloseq object

physeq16S <- phyloseq(
  otu_table(seqtab_nochim, taxa_are_rows = FALSE),
  sample_data(meta),
  tax_table(taxa)
)

print(physeq16S)

# Export OTU table, taxonomy table, and metadata

otu_export <- as.data.frame(t(otu_table(physeq16S)))
tax_export <- as.data.frame(tax_table(physeq16S))

write.table(
  otu_export,
  file = "asv_abundance.txt",
  sep = "\t",
  quote = FALSE,
  row.names = TRUE,
  col.names = NA
)

write.table(
  tax_export,
  file = "taxonomy.txt",
  sep = "\t",
  quote = FALSE,
  row.names = TRUE,
  col.names = NA
)

wb <- createWorkbook()

addWorksheet(wb, "Metadata")
writeData(wb, "Metadata", meta, rowNames = TRUE)

addWorksheet(wb, "OTU_Abundance")
writeData(wb, "OTU_Abundance", otu_export, rowNames = TRUE)

addWorksheet(wb, "Taxonomy")
writeData(wb, "Taxonomy", tax_export, rowNames = TRUE)

saveWorkbook(
  wb,
  "soil_bacteria_dataset.xlsx",
  overwrite = TRUE
)

# Save unrarefied phyloseq object

saveRDS(physeq16S, "physeq16S.rds")

# Sequencing depth summary before rarefaction

read_counts <- sort(sample_sums(physeq16S))

write.table(
  read_counts,
  file = "sample_read_counts_sorted.txt",
  sep = "\t",
  quote = FALSE,
  col.names = FALSE
)

print(summary(sample_sums(physeq16S)))

# Rarefy phyloseq object

set.seed(seed)

rarefied16S <- rarefy_even_depth(
  physeq16S,
  sample.size = sample_depth,
  rngseed = seed
)

# Check samples removed during rarefaction

removed_samples <- setdiff(
  sample_names(physeq16S),
  sample_names(rarefied16S)
)

if (length(removed_samples) > 0) {
  print(removed_samples)
}

print(summary(sample_sums(rarefied16S)))

# Save rarefied phyloseq object for downstream sampling-design analyses

saveRDS(rarefied16S, "rarefied16S.rds")
