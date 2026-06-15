#!/usr/bin/env Rscript
# ============================================================
# 08_pca.R — PCA of IP samples on peak coverage profiles
#
# Requires: proj1/prepared_data.Rdata
#           *.allcoverage.bed files in BAMS_DIR (00_config.R)
#
# Outputs (in proj1/):
#   IP.PCA.pdf   — PCA scatter plot coloured by sample
# ============================================================

suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
})

source("00_config.R")
setwd("proj1")
load("prepared_data.Rdata")

# ─── Load coverage BED files ─────────────────────────────────
message("Loading coverage BED files from: ", BAMS_DIR)
cov_files <- list.files(BAMS_DIR, pattern = "\\.allcoverage\\.bed$", full.names = TRUE)
if (length(cov_files) == 0) stop("No *.allcoverage.bed files found in: ", BAMS_DIR)

# Reference for column names (use first file)
ref_bed <- read.delim(cov_files[1], header = FALSE)

dd <- do.call(rbind, lapply(cov_files, function(f) {
    read.delim(f, header = FALSE)$V5
}))
rownames(dd) <- gsub("\\.allcoverage\\.bed$", "", basename(cov_files))
colnames(dd) <- ref_bed$V4
dd <- t(dd)

# ─── CPM normalisation ───────────────────────────────────────
dd_norm <- 1e6 * t(t(dd) / colSums(dd))

# Keep only IP samples for selected IDS
dd_norm <- dd_norm[, gsub("IP[0-9]+$", "", colnames(dd_norm)) %in% IDS]

# ─── PCA ─────────────────────────────────────────────────────
message("Running PCA on ", ncol(dd_norm), " IP samples...")
pca_res <- prcomp(dd_norm, scale. = TRUE)
summary(pca_res)

pct <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 1)

p <- pca_res$rotation[, 1:2] |>
    data.frame() |>
    mutate(sample_id = rownames(pca_res$rotation),
           ID = gsub("IP[12]", "", sample_id)) |>
    left_join(sampleinfo, by = "ID") |>
    ggplot(aes(PC1, PC2, color = ID_name)) +
    geom_point(size = 3) +
    theme_bw() +
    scale_color_manual(values = files.color) +
    labs(
        x     = paste0("PC1 (", pct[1], "%)"),
        y     = paste0("PC2 (", pct[2], "%)"),
        color = "Sample"
    )

ggsave(p, filename = "IP.PCA.pdf", height = 6, width = 7)
message("Done: 08_pca.R")