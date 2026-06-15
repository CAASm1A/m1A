#!/usr/bin/env Rscript
# ============================================================
# 02_annotate_peaks.R — Annotate m1A sites & compute AUG distance
#
# Requires: proj1/prepared_data.Rdata (from 01_prepare_data.R)
#
# Outputs (per sample, in proj1/):
#   <sample>.m1A.annotation.csv   — annotated peak table
# Saved object:
#   proj1/all_m1A.info.Rdata      — list of AnnoResult objects
# ============================================================

suppressPackageStartupMessages({
    library(dplyr)
    library(ChIPseeker)
    library(GenomicFeatures)
    library(GenomicRanges)
})

source("00_config.R")
setwd("proj")
load("prepared_data.Rdata")

# ─── Build start-codon (AUG) GRanges ────────────────────────
message("Building CDS / AUG reference...")
cds <- GenomicFeatures::cdsBy(txdb, by = "tx", use.name = TRUE)

allcds <- parallel::mclapply(names(cds), function(tx_id) {
    nar <- range(cds[[tx_id]])
    data.frame(
        seqnames    = as.character(seqnames(nar)),
        tx          = tx_id,
        start_codon = ifelse(strand(nar) == "+", start(nar), end(nar)),
        stop_codon  = ifelse(strand(nar) == "+", end(nar),   start(nar)),
        strand      = as.character(strand(nar))
    )
}, mc.cores = N_CORES) |> do.call(rbind, args = _)

atg <- GRanges(
    seqnames = allcds$seqnames,
    ranges   = IRanges(start = allcds$start_codon, end = allcds$start_codon + 1),
    strand   = allcds$strand,
    txname   = allcds$tx
)

# ─── Annotate each sample ────────────────────────────────────
message("Annotating m1A peaks...")
all_m1A.info <- lapply(names(files), function(sample_name) {
    message("  Processing: ", sample_name)
    anno <- annotatePeak(
        files[sample_name],
        TxDb    = txdb,
        tssRegion = c(0, 0),
        genomicAnnotationPriority = c("5UTR", "3UTR", "Exon", "Intron",
                                      "Downstream", "Intergenic", "Promoter")
    )

    # Distance to nearest AUG codon (matched by transcript ID)
    atg_idx <- match(anno@anno$transcriptId, atg$txname)
    anno@anno$distanceToAUG <- distance(anno@anno, atg[atg_idx])

    write.csv(
        as.data.frame(anno@anno),
        paste0(sample_name, ".m1A.annotation.csv"),
        row.names = FALSE
    )

    list(anno = anno, sample = sample_name)
})
names(all_m1A.info) <- names(files)

# ─── Save ────────────────────────────────────────────────────
save(all_m1A.info, atg, allcds, file = "all_m1A.info.Rdata")
message("Done: 02_annotate_peaks.R")