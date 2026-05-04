# 16S paired-end DADA2 pipeline
#
# This script processes paired-end 16S amplicon sequencing data using DADA2
#
# Main steps:
# - Filter and trim paired-end reads
# - Learn error rates
# - Dereplicate reads
# - Infer ASVs
# - Merge paired reads
# - Remove chimeras
# - Assign taxonomy using SILVA
# - Save ASV table, taxonomy table, sample names, and read tracking table


library(dada2)

path <- "data/fastq_16S"

silva_train_set <- "silva_nr99_v138.1_train_set.fa.gz"

trunc_len <- c(170, 150)
max_ee <- c(3, 3)
n_threads <- 8

# List paired-end FASTQ files

fnFs <- sort(list.files(path, pattern = "_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2_001.fastq.gz", full.names = TRUE))

# Extract sample names from forward-read filenames

sample_names <- sapply(
  strsplit(basename(fnFs), "_R1_001.fastq.gz"),
  `[`,
  1
)

# Create filtered-read output paths

filt_path <- file.path(path, "filtered16S")
dir.create(filt_path, showWarnings = FALSE, recursive = TRUE)

filtFs <- file.path(filt_path, paste0(sample_names, "_R1.filtered.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample_names, "_R2.filtered.fastq.gz"))

# Filter and trim reads

filter_out <- filterAndTrim(
  fnFs,
  filtFs,
  fnRs,
  filtRs,
  truncLen = trunc_len,
  maxEE = max_ee,
  multithread = n_threads
)

# Learn error rates

errF <- learnErrors(filtFs, multithread = n_threads)
errR <- learnErrors(filtRs, multithread = n_threads)

# Dereplicate reads

derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)

names(derepFs) <- sample_names
names(derepRs) <- sample_names

# Infer ASVs

dadaFs <- dada(derepFs, err = errF, multithread = n_threads)
dadaRs <- dada(derepRs, err = errR, multithread = n_threads)

# Merge paired reads

mergers <- mergePairs(
  dadaFs,
  derepFs,
  dadaRs,
  derepRs,
  verbose = FALSE
)

# Create ASV table and remove chimeras

seqtab <- makeSequenceTable(mergers)

seqtab_nochim <- removeBimeraDenovo(
  seqtab,
  method = "consensus",
  multithread = n_threads
)

# Track read counts through the pipeline

getN <- function(x) {
  sum(getUniques(x))
}

track <- cbind(
  filter_out,
  denoisedF = sapply(dadaFs, getN),
  denoisedR = sapply(dadaRs, getN),
  merged = sapply(mergers, getN),
  nonchim = rowSums(seqtab_nochim)
)

colnames(track)[1:2] <- c("input", "filtered")
rownames(track) <- sample_names

write.table(
  track,
  "tracking16S_paired.table.txt",
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

# Assign taxonomy

taxa <- assignTaxonomy(
  seqtab_nochim,
  silva_train_set,
  multithread = n_threads,
  tryRC = TRUE
)

taxa <- as.data.frame(taxa)
taxa[] <- lapply(taxa, as.character)

# Fill missing taxonomy with the nearest available higher rank

fill_na <- function(x, parent) {
  x[is.na(x)] <- parent[is.na(x)]
  x
}

taxa$Phylum <- fill_na(taxa$Phylum, taxa$Kingdom)
taxa$Class  <- fill_na(taxa$Class, taxa$Phylum)
taxa$Order  <- fill_na(taxa$Order, taxa$Class)
taxa$Family <- fill_na(taxa$Family, taxa$Order)
taxa$Genus  <- fill_na(taxa$Genus, taxa$Family)

taxa[] <- lapply(taxa, as.factor)
taxa <- as.matrix(taxa)

# Save outputs
saveRDS(seqtab_nochim, "seqtab16S_paired.nochim.rds")
saveRDS(taxa, "taxtable16S_paired.rds")

write.table(
  rownames(seqtab_nochim),
  "sample.names16S_paired.txt",
  sep = "\t",
  quote = FALSE,
  col.names = FALSE,
  row.names = FALSE
)
