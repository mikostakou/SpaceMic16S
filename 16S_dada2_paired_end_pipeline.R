# 16S paired-end DADA2 pipeline
#
# This script processes paired-end 16S amplicon sequencing data using DADA2.
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
out_path <- "results/dada2_16S"

dir.create(out_path, showWarnings = FALSE, recursive = TRUE)

filt_path <- file.path(out_path, "filtered16S")
dir.create(filt_path, showWarnings = FALSE, recursive = TRUE)

silva_train_set <- "data/reference/silva_nr99_v138.1_train_set.fa.gz"

trunc_len <- c(170, 150)
trim_left <- c(19, 20)
max_ee <- c(3, 3)
n_threads <- 8


# List paired-end FASTQ files

fnFs <- sort(list.files(
  path,
  pattern = "_R1_001\\.fastq\\.gz$",
  full.names = TRUE
))

fnRs <- sub("_R1_001\\.fastq\\.gz$", "_R2_001.fastq.gz", fnFs)

sample_names <- sub(
  "_R1_001\\.fastq\\.gz$",
  "",
  basename(fnFs)
)


# Create filtered-read output paths

filtFs <- file.path(filt_path, paste0(sample_names, "_R1.filtered.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample_names, "_R2.filtered.fastq.gz"))

names(filtFs) <- sample_names
names(filtRs) <- sample_names


# Filter and trim reads

filter_out <- filterAndTrim(
  fnFs, filtFs, fnRs, filtRs,
  truncLen = trunc_len,
  trimLeft = trim_left,
  maxN = 0,
  maxEE = max_ee,
  truncQ = 2,
  rm.phix = TRUE,
  compress = TRUE,
  multithread = n_threads,
  verbose = TRUE
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
  dadaFs, derepFs,
  dadaRs, derepRs,
  verbose = FALSE
)


# Create ASV table and remove chimeras

seqtab <- makeSequenceTable(mergers)

seqtab_nochim <- removeBimeraDenovo(
  seqtab,
  method = "consensus",
  multithread = n_threads,
  verbose = TRUE
)


# Track read counts through the pipeline

getN <- function(x) sum(getUniques(x))

track <- cbind(
  input = filter_out[, 1],
  filtered = filter_out[, 2],
  denoisedF = sapply(dadaFs, getN),
  denoisedR = sapply(dadaRs, getN),
  merged = sapply(mergers, getN),
  nonchim = rowSums(seqtab_nochim)
)

rownames(track) <- sample_names

write.table(
  track,
  file.path(out_path, "tracking16S.txt"),
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

saveRDS(seqtab_nochim, file.path(out_path, "seqtab16S.nochim.rds"))
saveRDS(taxa, file.path(out_path, "taxtable16S.rds"))

write.table(
  rownames(seqtab_nochim),
  file.path(out_path, "sample.names16S.txt"),
  sep = "\t",
  quote = FALSE,
  col.names = FALSE,
  row.names = FALSE
)
