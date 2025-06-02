#### 0101 analysis for freebayes ####
library(VariantAnnotation)
library(GenomicRanges)
library(rtracklayer)
library(GenomicFeatures)
library(dplyr)

rm(list = ls())

# === SETTINGS ===
vcf_dir <- "data/vcf/transcriptome/freebayes/"
vcf_files <- list.files(vcf_dir, pattern = "\\.vcf$", full.names = TRUE)
gtf_file <- "data/annotation/grch38/genes.gtf"
output_dir <- "results/vcf/plot/"

# === Load GTF Annotation Once ===
gtf <- import(gtf_file)

# Define annotation features
exons <- gtf[gtf$type == "exon"]
cds <- gtf[gtf$type == "CDS"]
utr <- gtf[gtf$type %in% c("five_prime_utr", "three_prime_utr")]
transcripts <- gtf[gtf$type == "transcript"]
tss <- resize(transcripts, width=1, fix="start")
promoters <- promoters(tss, upstream=2000, downstream=200)

txdb <- makeTxDbFromGFF(gtf_file, format="gtf")
introns <- unlist(intronsByTranscript(txdb, use.names=FALSE))

annotated_regions <- reduce(c(promoters, exons, cds, utr, introns))

# === Annotation Function ===
annotate_vcf <- function(vcf_path) {
  sample_name <- tools::file_path_sans_ext(basename(vcf_path))
  vcf <- readVcf(vcf_path, "GRCh38")
  snp_gr <- rowRanges(vcf)
  n_total <- length(snp_gr)
  
  get_overlap_count <- function(region_gr) {
    hits <- findOverlaps(snp_gr, region_gr)
    length(unique(queryHits(hits)))
  }
  
  counts <- list(
    promoter = get_overlap_count(promoters),
    exon = get_overlap_count(exons),
    cds = get_overlap_count(cds),
    utr = get_overlap_count(utr),
    intron = get_overlap_count(introns)
  )
  
  intergenic_hits <- GenomicRanges::setdiff(snp_gr, annotated_regions)
  counts$intergenic <- length(intergenic_hits)
  counts$total <- n_total
  
  proportions <- sapply(counts, function(x) x / n_total)
  
  df <- data.frame(
    Sample = sample_name,
    Region = names(proportions),
    Count = unlist(counts),
    Proportion = round(unlist(proportions), 4)
  )
  
  out_file <- file.path(output_dir, paste0(sample_name, "_snp_location_summary.csv"))
  write.csv(df, out_file, row.names = FALSE)
  
  return(df)
}

# === Run Batch Annotation ===
summary_list <- lapply(vcf_files, annotate_vcf)
merged_summary <- do.call(rbind, summary_list)

# Optionally write merged summary
write.csv(merged_summary, file.path(output_dir, "summary_snp_freebayes.csv"), row.names = FALSE)

#### 0102 analysis for bcftools ####
library(VariantAnnotation)
library(GenomicRanges)
library(rtracklayer)
library(GenomicFeatures)
library(dplyr)

rm(list = ls())

# === SETTINGS ===
vcf_dir <- "data/vcf/transcriptome/bcftools/"
vcf_files <- list.files(vcf_dir, pattern = "\\.vcf$", full.names = TRUE)
gtf_file <- "data/annotation/grch38/genes.gtf"
output_dir <- "results/vcf/plot/"

# === Load GTF Annotation Once ===
gtf <- import(gtf_file)

# Define annotation features
exons <- gtf[gtf$type == "exon"]
cds <- gtf[gtf$type == "CDS"]
utr <- gtf[gtf$type %in% c("five_prime_utr", "three_prime_utr")]
transcripts <- gtf[gtf$type == "transcript"]
tss <- resize(transcripts, width=1, fix="start")
promoters <- promoters(tss, upstream=2000, downstream=200)

txdb <- makeTxDbFromGFF(gtf_file, format="gtf")
introns <- unlist(intronsByTranscript(txdb, use.names=FALSE))

annotated_regions <- reduce(c(promoters, exons, cds, utr, introns))

# === Annotation Function ===
annotate_vcf <- function(vcf_path) {
  sample_name <- tools::file_path_sans_ext(basename(vcf_path))
  vcf <- readVcf(vcf_path, "GRCh38")
  snp_gr <- rowRanges(vcf)
  n_total <- length(snp_gr)
  
  get_overlap_count <- function(region_gr) {
    hits <- findOverlaps(snp_gr, region_gr)
    length(unique(queryHits(hits)))
  }
  
  counts <- list(
    promoter = get_overlap_count(promoters),
    exon = get_overlap_count(exons),
    cds = get_overlap_count(cds),
    utr = get_overlap_count(utr),
    intron = get_overlap_count(introns)
  )
  
  intergenic_hits <- GenomicRanges::setdiff(snp_gr, annotated_regions)
  counts$intergenic <- length(intergenic_hits)
  counts$total <- n_total
  
  proportions <- sapply(counts, function(x) x / n_total)
  
  df <- data.frame(
    Sample = sample_name,
    Region = names(proportions),
    Count = unlist(counts),
    Proportion = round(unlist(proportions), 4)
  )
  
  out_file <- file.path(output_dir, paste0(sample_name, "_snp_location_summary.csv"))
  write.csv(df, out_file, row.names = FALSE)
  
  return(df)
}

# === Run Batch Annotation ===
summary_list <- lapply(vcf_files, annotate_vcf)
merged_summary <- do.call(rbind, summary_list)

# Optionally write merged summary
write.csv(merged_summary, file.path(output_dir, "summary_snp_bcftools.csv"), row.names = FALSE)


#### 0103 analysis for cellsnp ####
library(VariantAnnotation)
library(GenomicRanges)
library(rtracklayer)
library(GenomicFeatures)
library(dplyr)

rm(list = ls())

# === SETTINGS ===
vcf_dir <- "data/vcf/transcriptome/cellsnp/"
vcf_files <- list.files(vcf_dir, pattern = "\\.vcf$", full.names = TRUE)
gtf_file <- "data/annotation/grch38/genes.gtf"
output_dir <- "results/vcf/plot/"

# === Load GTF Annotation Once ===
gtf <- import(gtf_file)

# Define annotation features
exons <- gtf[gtf$type == "exon"]
cds <- gtf[gtf$type == "CDS"]
utr <- gtf[gtf$type %in% c("five_prime_utr", "three_prime_utr")]
transcripts <- gtf[gtf$type == "transcript"]
tss <- resize(transcripts, width=1, fix="start")
promoters <- promoters(tss, upstream=2000, downstream=200)

txdb <- makeTxDbFromGFF(gtf_file, format="gtf")
introns <- unlist(intronsByTranscript(txdb, use.names=FALSE))

annotated_regions <- reduce(c(promoters, exons, cds, utr, introns))

# === Annotation Function ===
annotate_vcf <- function(vcf_path) {
  sample_name <- tools::file_path_sans_ext(basename(vcf_path))
  vcf <- readVcf(vcf_path, "GRCh38")
  snp_gr <- rowRanges(vcf)
  n_total <- length(snp_gr)
  
  get_overlap_count <- function(region_gr) {
    hits <- findOverlaps(snp_gr, region_gr)
    length(unique(queryHits(hits)))
  }
  
  counts <- list(
    promoter = get_overlap_count(promoters),
    exon = get_overlap_count(exons),
    cds = get_overlap_count(cds),
    utr = get_overlap_count(utr),
    intron = get_overlap_count(introns)
  )
  
  intergenic_hits <- GenomicRanges::setdiff(snp_gr, annotated_regions)
  counts$intergenic <- length(intergenic_hits)
  counts$total <- n_total
  
  proportions <- sapply(counts, function(x) x / n_total)
  
  df <- data.frame(
    Sample = sample_name,
    Region = names(proportions),
    Count = unlist(counts),
    Proportion = round(unlist(proportions), 4)
  )
  
  out_file <- file.path(output_dir, paste0(sample_name, "_snp_location_summary.csv"))
  write.csv(df, out_file, row.names = FALSE)
  
  return(df)
}

# === Run Batch Annotation ===
summary_list <- lapply(vcf_files, annotate_vcf)
merged_summary <- do.call(rbind, summary_list)

# Optionally write merged summary
write.csv(merged_summary, file.path(output_dir, "summary_snp_cellsnp.csv"), row.names = FALSE)

#### 0104 analysis for SNP array ####
library(VariantAnnotation)
library(GenomicRanges)
library(rtracklayer)
library(GenomicFeatures)
library(dplyr)

rm(list = ls())

# === SETTINGS ===
vcf_dir <- "data/vcf/array/"
vcf_files <- list.files(vcf_dir, pattern = "\\.vcf$", full.names = TRUE)
gtf_file <- "data/annotation/grch38/genes.gtf"
output_dir <- "results/vcf/plot/"

# === Load GTF Annotation Once ===
gtf <- import(gtf_file)

# Define annotation features
exons <- gtf[gtf$type == "exon"]
cds <- gtf[gtf$type == "CDS"]
utr <- gtf[gtf$type %in% c("five_prime_utr", "three_prime_utr")]
transcripts <- gtf[gtf$type == "transcript"]
tss <- resize(transcripts, width=1, fix="start")
promoters <- promoters(tss, upstream=2000, downstream=200)

txdb <- makeTxDbFromGFF(gtf_file, format="gtf")
introns <- unlist(intronsByTranscript(txdb, use.names=FALSE))

annotated_regions <- reduce(c(promoters, exons, cds, utr, introns))

vcf_path <- vcf_files[1]
# === Annotation Function ===
annotate_vcf <- function(vcf_path) {
  sample_name <- tools::file_path_sans_ext(basename(vcf_path))
  vcf <- readVcf(vcf_path, "GRCh38")
  snp_gr <- rowRanges(vcf)
  n_total <- length(snp_gr)
  
  get_overlap_count <- function(region_gr) {
    hits <- findOverlaps(snp_gr, region_gr)
    length(unique(queryHits(hits)))
  }
  
  counts <- list(
    promoter = get_overlap_count(promoters),
    exon = get_overlap_count(exons),
    cds = get_overlap_count(cds),
    utr = get_overlap_count(utr),
    intron = get_overlap_count(introns)
  )
  
  intergenic_hits <- GenomicRanges::setdiff(snp_gr, annotated_regions)
  counts$intergenic <- length(intergenic_hits)
  counts$total <- n_total
  
  proportions <- sapply(counts, function(x) x / n_total)
  
  df <- data.frame(
    Sample = sample_name,
    Region = names(proportions),
    Count = unlist(counts),
    Proportion = round(unlist(proportions), 4)
  )
  
  out_file <- file.path(output_dir, paste0(sample_name, "_snp_location_summary.csv"))
  write.csv(df, out_file, row.names = FALSE)
  
  return(df)
}

# === Run Batch Annotation ===
summary_list <- lapply(vcf_files, annotate_vcf)
merged_summary <- do.call(rbind, summary_list)

# Optionally write merged summary
write.csv(merged_summary, file.path(output_dir, "summary_snp_snparray.csv"), row.names = FALSE)



#### 0201 concordance for freebayes ####
library(VariantAnnotation)
library(GenomicRanges)
library(dplyr)
rm(list = ls())

# Define sample IDs to process
sample_ids <- c("IDCM1", "IDCM2", "IDCM3","IDCM4", "IDCM5", "IDCM6")  # Replace with your actual sample names

# Base paths
vcf_array_base <- "data/vcf/array/snparray_"
vcf_rnaseq_base <- "data/vcf/transcriptome/freebayes/snp_rna_freebayes_"

# Output path
output_dir <- "results/vcf/concordance/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Function to compute genotype concordance for a single sample
process_sample <- function(sample_id) {
  # File paths
  vcf_array_path <- paste0(vcf_array_base, sample_id, ".vcf")
  vcf_rnaseq_path <- paste0(vcf_rnaseq_base, sample_id, ".vcf")
  
  # Load VCFs
  vcf_array <- readVcf(vcf_array_path, genome = "GRCh38")
  vcf_rnaseq <- readVcf(vcf_rnaseq_path, genome = "GRCh38")
  
  # Extract SNPs
  gr_array <- rowRanges(vcf_array)[isSNV(vcf_array)]
  gr_rnaseq <- rowRanges(vcf_rnaseq)[isSNV(vcf_rnaseq)]
  
  # Strip seqinfo to avoid overlap errors
  seqinfo(gr_array) <- Seqinfo(seqlevels(gr_array))
  seqinfo(gr_rnaseq) <- Seqinfo(seqlevels(gr_rnaseq))
  
  # Find overlapping SNPs
  hits <- findOverlaps(gr_array, gr_rnaseq)
  
  # Genotype comparison
  geno_array <- geno(vcf_array)$GT[queryHits(hits)]
  geno_rnaseq <- geno(vcf_rnaseq)$GT[subjectHits(hits)]
  
  comparison_df <- data.frame(
    chr = as.character(seqnames(gr_array)[queryHits(hits)]),
    pos = start(gr_array)[queryHits(hits)],
    GT_array = geno_array,
    GT_rnaseq = geno_rnaseq
  )
  
  comparison_df$match <- comparison_df$GT_array == comparison_df$GT_rnaseq
  concordance <- mean(comparison_df$match, na.rm = TRUE)
  
  # Save per-sample results
  write.csv(comparison_df,
            file = file.path(output_dir, paste0("snp_concordance_freebayes_", sample_id, ".csv")),
            row.names = FALSE)
  
  write.csv(data.frame(sample = sample_id, concordance = concordance),
            file = file.path(output_dir, paste0("snp_concordance_freebayes_", sample_id, "_summary.csv")),
            row.names = FALSE)
  
  # Optionally return concordance
  return(data.frame(sample = sample_id, concordance = concordance, n_overlap = nrow(comparison_df)))
}

# Run the loop
summary_list <- lapply(sample_ids, process_sample)

# Combine and save full summary
final_summary <- do.call(rbind, summary_list)
write.csv(final_summary,
          file = file.path(output_dir, "snp_concordance_freebayes_summary.csv"),
          row.names = FALSE)

print(final_summary)



#### 0202 concordance for cellsnp ####
library(VariantAnnotation)
library(GenomicRanges)
library(dplyr)
rm(list = ls())

# Define sample IDs to process
sample_ids <- c("IDCM1", "IDCM2", "IDCM3","IDCM4", "IDCM5", "IDCM6")  # Replace with your actual sample names

# Base paths
vcf_array_base <- "data/vcf/array/snparray_"
vcf_rnaseq_base <- "data/vcf/transcriptome/cellsnp/snp_rna_cellSNP_"

# Output path
output_dir <- "results/vcf/concordance/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Function to compute genotype concordance for a single sample
process_sample <- function(sample_id) {
  # File paths
  vcf_array_path <- paste0(vcf_array_base, sample_id, ".vcf")
  vcf_rnaseq_path <- paste0(vcf_rnaseq_base, sample_id, "_GT", ".vcf")
  
  # Load VCFs
  vcf_array <- readVcf(vcf_array_path, genome = "GRCh38")
  vcf_rnaseq <- readVcf(vcf_rnaseq_path, genome = "GRCh38")
  
  # Extract SNPs
  gr_array <- rowRanges(vcf_array)[isSNV(vcf_array)]
  gr_rnaseq <- rowRanges(vcf_rnaseq)[isSNV(vcf_rnaseq)]
  
  # Strip seqinfo to avoid overlap errors
  # seqinfo(gr_array) <- Seqinfo(seqlevels(gr_array))
  # seqinfo(gr_rnaseq) <- Seqinfo(seqlevels(gr_rnaseq))
  
  # Find overlapping SNPs
  hits <- findOverlaps(gr_array, gr_rnaseq)
  
  # Genotype comparison
  geno_array <- geno(vcf_array)$GT[queryHits(hits)]
  geno_rnaseq <- geno(vcf_rnaseq)$GT[subjectHits(hits)]
  
  comparison_df <- data.frame(
    chr = as.character(seqnames(gr_array)[queryHits(hits)]),
    pos = start(gr_array)[queryHits(hits)],
    GT_array = geno_array,
    GT_rnaseq = geno_rnaseq
  )
  
  comparison_df$match <- comparison_df$GT_array == comparison_df$GT_rnaseq
  concordance <- mean(comparison_df$match, na.rm = TRUE)
  
  # Save per-sample results
  write.csv(comparison_df,
            file = file.path(output_dir, paste0("snp_concordance_cellSNP_", sample_id, "_GT", ".csv")),
            row.names = FALSE)
  
  write.csv(data.frame(sample = sample_id, concordance = concordance),
            file = file.path(output_dir, paste0("snp_concordance_cellSNP_", sample_id, "_GT", "_summary.csv")),
            row.names = FALSE)
  
  # Optionally return concordance
  return(data.frame(sample = sample_id, concordance = concordance, n_overlap = nrow(comparison_df)))
}

# Run the loop
summary_list <- lapply(sample_ids, process_sample)

# Combine and save full summary
final_summary <- do.call(rbind, summary_list)
write.csv(final_summary,
          file = file.path(output_dir, "snp_concordance_cellSNP_GT_summary.csv"),
          row.names = FALSE)

print(final_summary)

#### 0203 concordance for bcftools####
library(VariantAnnotation)
library(GenomicRanges)
library(dplyr)
rm(list = ls())

# Define sample IDs to process
sample_ids <- c("IDCM1", "IDCM2", "IDCM3","IDCM4", "IDCM5", "IDCM6")  # Replace with your actual sample names

# Base paths
vcf_array_base <- "data/vcf/array/snparray_"
vcf_rnaseq_base <- "data/vcf/transcriptome/bcftools/snp_rna_bcftools_"

# Output path
output_dir <- "results/vcf/concordance/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Function to compute genotype concordance for a single sample
process_sample <- function(sample_id) {
  # File paths
  vcf_array_path <- paste0(vcf_array_base, sample_id, ".vcf")
  vcf_rnaseq_path <- paste0(vcf_rnaseq_base, sample_id, "_GT", ".vcf")
  
  # Load VCFs
  vcf_array <- readVcf(vcf_array_path, genome = "GRCh38")
  vcf_rnaseq <- readVcf(vcf_rnaseq_path, genome = "GRCh38")
  
  # Extract SNPs
  gr_array <- rowRanges(vcf_array)[isSNV(vcf_array)]
  gr_rnaseq <- rowRanges(vcf_rnaseq)[isSNV(vcf_rnaseq)]
  
  # Strip seqinfo to avoid overlap errors
  seqinfo(gr_array) <- Seqinfo(seqlevels(gr_array))
  seqinfo(gr_rnaseq) <- Seqinfo(seqlevels(gr_rnaseq))
  
  # Find overlapping SNPs
  hits <- findOverlaps(gr_array, gr_rnaseq)
  
  # Genotype comparison
  geno_array <- geno(vcf_array)$GT[queryHits(hits)]
  geno_rnaseq <- geno(vcf_rnaseq)$GT[subjectHits(hits)]
  
  comparison_df <- data.frame(
    chr = as.character(seqnames(gr_array)[queryHits(hits)]),
    pos = start(gr_array)[queryHits(hits)],
    GT_array = geno_array,
    GT_rnaseq = geno_rnaseq
  )
  
  comparison_df$match <- comparison_df$GT_array == comparison_df$GT_rnaseq
  concordance <- mean(comparison_df$match, na.rm = TRUE)
  
  # Save per-sample results
  write.csv(comparison_df,
            file = file.path(output_dir, paste0("snp_concordance_bcftools_", sample_id, "_GT", ".csv")),
            row.names = FALSE)
  
  write.csv(data.frame(sample = sample_id, concordance = concordance),
            file = file.path(output_dir, paste0("snp_concordance_bcftools_", sample_id, "_GT", "_summary.csv")),
            row.names = FALSE)
  
  # Optionally return concordance
  return(data.frame(sample = sample_id, concordance = concordance, n_overlap = nrow(comparison_df)))
}

# Run the loop
summary_list <- lapply(sample_ids, process_sample)

# Combine and save full summary
final_summary <- do.call(rbind, summary_list)
write.csv(final_summary,
          file = file.path(output_dir, "snp_concordance_bcftools_GT_summary.csv"),
          row.names = FALSE)

print(final_summary)

#### 03 concordance plot Total SNP vs Total SNP ####
library(knitr)
library(kableExtra)
library(ggplot2)

final_summary <- read.csv("results/vcf/concordance/snp_concordance.csv")

ggplot(final_summary, aes(x = sample, y = concordance, fill = method)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = round(concordance, 2)),
            position = position_dodge(width = 0.8),
            vjust = -0.4, size = 3.0) +
  ylim(0, 1) +
  scale_fill_manual(
    values = c(
      "cellSNP" = "#1b9e77",     # teal green
      "FreeBayes" = "#d95f02",   # orange
      "BCFtools" = "#7570b3"     # purple
    )
  ) +
  labs(title = "Genotype Concordance: RNA-Seq vs SNP Array",
       x = "Sample",
       y = "Concordance Rate",
       fill = "Variant Caller") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(size = 10, angle = 0, hjust = 0.5, vjust = 5),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 15, face = "bold"),
    axis.title.y = element_text(size = 15, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
  )
ggsave("results/vcf/concordance/snp_concordance_plot.pdf", width = 6, height = 4, dpi = 300)

#### 04 concordance plot Genic SNP vs Genic SNP ####
library(knitr)
library(kableExtra)
library(ggplot2)

final_summary <- read.csv("results/vcf/concordance/Genic/snp_concordance.csv")

ggplot(final_summary, aes(x = sample, y = concordance, fill = method)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = round(concordance, 2)),
            position = position_dodge(width = 0.8),
            vjust = -0.4, size = 3.0) +
  ylim(0, 1) +
  scale_fill_manual(
    values = c(
      "cellSNP" = "#1b9e77",     # teal green
      "FreeBayes" = "#d95f02",   # orange
      "BCFtools" = "#7570b3"     # purple
    )
  ) +
  labs(title = "Genotype Concordance: RNA-Seq vs SNP Array",
       x = "Sample",
       y = "Concordance Rate",
       fill = "Variant Caller") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(size = 10, angle = 0, hjust = 0.5, vjust = 5),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 15, face = "bold"),
    axis.title.y = element_text(size = 15, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
  )
ggsave("results/vcf/concordance/Genic/snp_concordance_plot.pdf", width = 6, height = 4, dpi = 300)
