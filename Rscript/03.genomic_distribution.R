#!/usr/bin/env Rscript
# ============================================================
# 03_genomic_distribution.R — Genomic distribution of m1A sites
#
# Requires: proj1/prepared_data.Rdata
#           proj1/all_m1A.info.Rdata
#           PLOT_FNS (plot_fns.R) defined in 00_config.R
#
# Outputs (in proj1/):
#   distribution.pdf              — mRNA region distribution
#   peaks.distribution.pdf        — plotPeakProf2 gene body profile
#   peaks.distribution2.pdf       — CDS-centred profile with flanks
#   AUG.distribution.pdf          — per-sample distance to AUG
#   stop_codon.distribution.pdf   — per-sample distance to stop codon
#   TSS.distribution.pdf          — per-sample distance to TSS
#   AUG.all.distribution.pdf      — all samples merged, AUG
#   stop_codon.all.distribution.pdf
#   TSS.all.distribution.pdf
#   AUG.exons.distribution.pdf    — AUG distance split by exon rank
# ============================================================

suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(ChIPseeker)
    library(GenomicFeatures)
    library(GenomicRanges)
})

source("00_config.R")
source(PLOT_FNS)
setwd("proj")
load("prepared_data.Rdata")
load("all_m1A.info.Rdata")   # provides: all_m1A.info, atg, allcds

# ─── 1. mRNA region distribution (dis_mRNA) ─────────────────
message("Plotting mRNA region distribution...")
pdf("distribution.pdf", 8, 5)
p <- dis_mRNA(files = files, txdb = txdb, label = names(files), ext = 0, downstream = 0)
print(p + scale_color_manual(values = files.color))
dev.off()

# ─── 2. Gene-body profile (ChIPseeker plotPeakProf2) ─────────
message("Plotting gene body profiles...")
pdf("peaks.distribution.pdf")
for (sample_name in names(files)) {
    plotPeakProf2(
        files[sample_name],
        TxDb       = txdb,
        upstream   = rel(0.2),
        downstream = rel(0.2),
        conf       = 0.95,
        by         = "gene",
        type       = "body",
        nbin       = 20
    )
}
dev.off()

# ─── 3. CDS-centred profile ──────────────────────────────────
message("Plotting CDS-centred profile...")
allcds_gr <- GRanges(
    seqnames = allcds$seqnames,
    ranges   = IRanges(start = allcds$start, end = allcds$end),
    strand   = allcds$strand
)

upstream_bp   <- 2000
downstream_bp <- 2000
bins          <- 800
breaks <- c(
    1,
    bins * 0.1 * upstream_bp / 1000,
    bins + bins * 0.1 * upstream_bp / 1000,
    bins + 2 * bins * 0.1 * upstream_bp / 1000
)

getTagMatrix_cds <- function(peak_file) {
    getTagMatrix(
        peak_file,
        windows    = makeBioRegionFromGranges(allcds_gr, type = "body", by = "gene"),
        nbin       = bins,
        upstream   = upstream_bp,
        downstream = downstream_bp
    )
}
tagMatrixList <- lapply(as.list(files), getTagMatrix_cds)

p_cds <- plotPeakProf(tagMatrixList) +
    scale_x_continuous(
        breaks = breaks,
        labels = c("5UTR", "CDS_start", "CDS_end", "3UTR")
    )
pdf("peaks.distribution2.pdf", 8, 6)
print(p_cds)
dev.off()

# ─── 4. Landmark distance plots ──────────────────────────────

# Stop codon GRanges
stop_codon <- GRanges(
    seqnames = allcds$seqnames,
    ranges   = IRanges(start = allcds$stop_codon, end = allcds$stop_codon + 1),
    strand   = allcds$strand,
    txname   = allcds$tx
)

# TSS GRanges
tx <- transcripts(txdb)
tss <- GRanges(
    seqnames = seqnames(tx),
    ranges   = IRanges(
        start = ifelse(strand(tx) == "+", start(tx),     end(tx) - 1),
        end   = ifelse(strand(tx) == "+", start(tx) + 1, end(tx))
    ),
    strand  = strand(tx),
    txname  = tx$tx_name
)

# Filter AUG sites away from chromosome edges
atg$chr_len <- width(genome[seqnames(atg)])
atg_filtered <- atg[(atg$chr_len - end(atg)) > 300 & start(atg) > 300]

message("Plotting distance to AUG, stop codon, TSS (per sample + merged)...")

# Per-sample plots
pdf("AUG.distribution.pdf", 10, 6)
print(plot_aug(files, atg_filtered, extend = AUG_EXTEND,
               xlab = "Distance from AUG codon (nts)", colors = files.color) +
      scale_x_continuous(n.breaks = 10))
dev.off()

pdf("stop_codon.distribution.pdf", 10, 6)
print(plot_aug(files, stop_codon, extend = STOP_EXTEND,
               xlab = "Distance from stop codon (nts)", colors = files.color) +
      scale_x_continuous(n.breaks = 10))
dev.off()

pdf("TSS.distribution.pdf", 7, 5)
print(plot_aug(files, tss, extend = TSS_EXTEND,
               xlab = "Distance from TSS (nts)", colors = files.color) +
      coord_cartesian(xlim = c(-200, 1000)) +
      scale_x_continuous(n.breaks = 10))
dev.off()

# All-samples-merged plots
all_peaks_gr <- unlist(GRangesList(lapply(files, readPeakFile)), use.names = FALSE)

pdf("AUG.all.distribution.pdf", 10, 6)
print(plot_aug(all_peaks_gr, atg_filtered, extend = AUG_EXTEND,
               xlab = "Distance from AUG codon (nts)") +
      scale_x_continuous(n.breaks = 10))
dev.off()

pdf("stop_codon.all.distribution.pdf", 10, 6)
print(plot_aug(all_peaks_gr, stop_codon, extend = STOP_EXTEND,
               xlab = "Distance from stop codon (nts)") +
      scale_x_continuous(n.breaks = 10))
dev.off()

pdf("TSS.all.distribution.pdf", 7, 5)
print(plot_aug(all_peaks_gr, tss, extend = TSS_EXTEND,
               xlab = "Distance from TSS (nts)") +
      coord_cartesian(xlim = c(-200, 1000)) +
      scale_x_continuous(n.breaks = 10))
dev.off()

# ─── 5. AUG distance split by exon rank ──────────────────────
message("Plotting AUG distance by exon rank...")
cds_obj <- GenomicFeatures::cdsBy(txdb, by = "tx", use.name = TRUE)
allAUG.pos <- parallel::mclapply(names(cds_obj), function(tx_id) {
    data.frame(txname = tx_id, exon_rank = cds_obj[[tx_id]]$exon_rank[1])
}, mc.cores = N_CORES) |> do.call(rbind, args = _)

plot_aug_by_exon <- function(peak_file, atg_gr, title = "",
                             extend = AUG_EXTEND, xlab = "") {
    atg_ext  <- resize(atg_gr, width = width(atg_gr) + extend, fix = "center")
    over     <- findOverlaps(readPeakFile(peak_file), atg_ext, type = "within")
    query    <- as.data.frame(readPeakFile(peak_file)[from(over)])
    colnames(query) <- paste0("q.", colnames(query))
    target   <- as.data.frame(atg_ext[to(over)])

    pdata <- cbind(query, target) |>
        distinct(q.start, q.end, .keep_all = TRUE) |>
        mutate(pos = ifelse(strand == "+",
                            q.start - (start + extend / 2),
                            (end - extend / 2) - q.start)) |>
        left_join(allAUG.pos, by = c("txname")) |>
        dplyr::filter(exon_rank < 3) |>
        group_by(exon_rank) |>
        mutate(total = n()) |>
        dplyr::filter(total > 5)

    ggplot(pdata, aes(pos, color = factor(exon_rank))) +
        geom_density(aes(y = 100 * after_stat(count) / sum(after_stat(count)))) +
        theme_classic() +
        labs(x = xlab, y = "m1A peaks density") +
        BoldTheme() +
        ggtitle(title)
}

pdf("AUG.exons.distribution.pdf", 7, 5)
for (sample_name in names(files)) {
    print(plot_aug_by_exon(files[sample_name], atg_filtered,
                           title = sample_name,
                           xlab = "Distance from AUG codon (nts)"))
}
dev.off()

message("Done: 03_genomic_distribution.R")