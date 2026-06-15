#!/usr/bin/env Rscript
# ============================================================
# 01_prepare_data.R — Load & preprocess all shared data objects
#
# Outputs (saved to proj1/prepared_data.Rdata):
#   genome        DNAStringSet — reference genome sequences
#   txdb          TxDb object
#   sampleinfo    data.frame — sample metadata
#   files         named character vector — per-sample BED paths
#   files.color   named character vector — per-sample plot colours
#   rna           long-format transcript TPM data.frame
#   gene.rna      wide-format gene TPM matrix (samples × genes)
# ============================================================

suppressPackageStartupMessages({
    library(dplyr)
    library(Biostrings)
    library(GenomicFeatures)
    library(reshape2)
})

source("00_config.R")

OUT_DIR <- "proj"
dir.create(OUT_DIR, showWarnings = FALSE)
setwd(OUT_DIR)

# ─── Reference genome & annotation ──────────────────────────
message("Loading genome FASTA...")
genome <- readDNAStringSet(GENOME_FA)

message("Building TxDb from GTF...")
txdb <- makeTxDbFromGFF(GTF_PATH)

# ─── Sample metadata ─────────────────────────────────────────
sampleinfo <- read.delim(SAMPLE_INFO)

# ─── m1A site BED files ──────────────────────────────────────
all_beds  <- list.files(SITE_BED_DIR, pattern = "\\.filter\\.site\\.bed$", full.names = TRUE)
bed_labels <- gsub("\\.filter\\.site\\.bed$", "", basename(all_beds))
names(all_beds) <- bed_labels

files       <- all_beds[IDS]
names(files) <- sampleinfo$ID_name[match(IDS, sampleinfo$ID)]

files.color <- sampleinfo$color[match(IDS, sampleinfo$ID)]
names(files.color) <- sampleinfo$ID_name[match(IDS, sampleinfo$ID)]

message("Sample files:")
print(files)

# ─── Transcript-level TPM ────────────────────────────────────
message("Loading transcript TPM...")
rna.raw <- read.csv(TRANS_TPM, row.names = 1) |> as.matrix()

rna <- rna.raw |>
    reshape2::melt() |>
    mutate(ID = gsub("_|Input.+", "", Var2)) |>
    group_by(Var1, ID) |>
    summarise(value = mean(value, na.rm = TRUE), .groups = "drop") |>
    mutate(gene = gsub("\\..+", "", Var1)) |>
    left_join(sampleinfo, by = "ID")

# ─── Gene-level TPM ──────────────────────────────────────────
message("Loading gene TPM...")
gene.rna.raw <- read.csv(GENE_TPM, row.names = 1) |> as.matrix()
gene.rna.input <- gene.rna.raw[, grepl("Input", colnames(gene.rna.raw))]

gene.rna <- gene.rna.input |>
    reshape2::melt() |>
    mutate(ID = gsub("_|Input.+", "", Var2)) |>
    dplyr::filter(ID %in% IDS) |>
    group_by(Var1, ID) |>
    summarise(value = mean(value, na.rm = TRUE), .groups = "drop") |>
    mutate(gene = gsub("\\..+", "", Var1)) |>
    left_join(sampleinfo, by = "ID")

gene.rna <- reshape2::acast(gene.rna, Var1 ~ ID_name, value.var = "value")

# ─── Save shared objects ─────────────────────────────────────
message("Saving prepared data to proj1/prepared_data.Rdata")
save(genome, txdb, sampleinfo, files, files.color, rna, gene.rna,
     file = "prepared_data.Rdata")

message("Done: 01_prepare_data.R")