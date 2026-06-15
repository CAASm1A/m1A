#!/usr/bin/env Rscript
# ============================================================
# 04_motif_analysis.R — m1A sequence motif discovery
#
# Requires: proj1/prepared_data.Rdata
#           PLOT_FNS (get_fa, get_motif) defined in 00_config.R
#
# For each sample:
#   1. Extracts ±MOTIF_EXT bp sequences around m1A summit
#   2. Builds a matched background set
#   3. Runs motif discovery (HOMER/similar via get_motif)
#
# Outputs (in proj1/m1A_motif/):
#   <sample>.fa        — foreground sequences
#   <sample>.bg.fa     — background sequences
#   <sample>/          — motif results directory
# ============================================================

suppressPackageStartupMessages({
    library(Biostrings)
})

source("00_config.R")
source(PLOT_FNS)
setwd("proj")
load("prepared_data.Rdata")

motif_dir <- "m1A_motif"
dir.create(motif_dir, showWarnings = FALSE)

message("Running motif analysis for ", length(files), " samples...")
for (sample_name in names(files)) {
    message("  Processing: ", sample_name)
    fa_prefix  <- file.path(motif_dir, sample_name)
    fa_file    <- paste0(fa_prefix, ".fa")
    bg_file    <- paste0(fa_prefix, ".bg.fa")
    out_dir    <- fa_prefix

    seq_result <- get_fa(files[sample_name], genome, ext = MOTIF_EXT)
    write(seq_result$fa, fa_file)
    write(seq_result$bg, bg_file)
    get_motif(fa_file, bg_file, out_dir)
}

message("Done: 04_motif_analysis.R")