#!/usr/bin/env Rscript
# ============================================================
# 07_go_enrichment.R — GO enrichment & expression heatmap
#
# Requires: proj1/prepared_data.Rdata
#           proj1/all_m1A.info.Rdata
#
# Outputs (in proj1/):
#   <sample>.GO.csv       — GO enrichment results per sample
#   m1A.EXP.heatmap.pdf   — expression heatmap of m1A-modified genes
# ============================================================

suppressPackageStartupMessages({
    library(dplyr)
    library(clusterProfiler)
    library(org.At.tair.db)
    library(pheatmap)
})

source("00_config.R")
setwd("proj")
load("prepared_data.Rdata")
load("all_m1A.info.Rdata")

# ─── 1. Per-sample GO enrichment ─────────────────────────────
message("Running GO enrichment...")
lapply(names(all_m1A.info), function(sample_name) {
    message("  Sample: ", sample_name)
    gene_ids <- unique(all_m1A.info[[sample_name]]$anno@anno$geneId)

    ego <- enrichGO(
        gene          = gene_ids,
        OrgDb         = org.At.tair.db,
        keyType       = "TAIR",
        ont           = "ALL",
        pAdjustMethod = "none"
    )
    write.csv(ego@result, paste0(sample_name, ".GO.csv"), row.names = FALSE)
    ego
})

# ─── 2. Expression heatmap for m1A-modified genes ────────────
message("Plotting expression heatmap for m1A genes...")
m1A_genes <- sapply(names(all_m1A.info), function(sample_name) {
    all_m1A.info[[sample_name]]$anno@anno$geneId
}) |> unlist() |> unique()

pdata <- gene.rna[intersect(m1A_genes, rownames(gene.rna)), ]
pdata <- pdata[rowSums(pdata > 0) > 3, ]

pheatmap::pheatmap(
    pdata,
    scale         = "row",
    show_rownames = FALSE,
    breaks        = seq(-2, 2, length.out = 100),
    filename      = "m1A.EXP.heatmap.pdf"
)

message("Done: 07_go_enrichment.R")