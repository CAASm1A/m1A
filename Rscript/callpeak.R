#!/usr/bin/env Rscript
# ============================================================
# callpeak.R — MACS2 peak calling & annotation for m1A-seq
# Organism: Arabidopsis thaliana (TAIR10)
#
# Usage:
#   Rscript callpeak.R <input_prefix> <IP_prefix>
#
# Arguments:
#   input_prefix   Prefix matching Input BAM(s) in current directory
#   IP_prefix      Prefix matching IP BAM(s) in current directory
#
# Output (written to ./<IP_prefix>/):
#   peaks.annotation.csv       Per-peak genomic annotation table
#   peaks.distribution.pdf     Peak profile over gene body (plotPeakProf2)
#   peaks.distribution2.pdf    Peak profile over CDS with flanks (plotPeakProf)
# ============================================================

suppressPackageStartupMessages({
    library(dplyr)
    library(ChIPseeker)
    library(GenomicFeatures)
    library(ggplot2)
})

# ─── User-configurable paths ────────────────────────────────
GTF_PATH   <- "/path/to/TAIR10/genes/genes.gtf"   # GTF annotation file
GENOME_SIZE <- 119667750                            # Arabidopsis effective genome size (bp)

# ─── Parse command-line arguments ───────────────────────────
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
    stop("Usage: Rscript callpeak.R <input_prefix> <IP_prefix>")
}
input_prefix <- args[1]
ip_prefix    <- args[2]

# ─── Build TxDb from GTF ────────────────────────────────────
message("Building TxDb from: ", GTF_PATH)
txdb <- makeTxDbFromGFF(GTF_PATH)

# ─── Locate BAM files ────────────────────────────────────────
all_bams    <- list.files(".", pattern = "\\.bam$")
input_bams  <- grep(input_prefix, all_bams, value = TRUE)
ip_bams     <- grep(ip_prefix,    all_bams, value = TRUE)

if (length(input_bams) == 0) stop("No Input BAMs found matching: ", input_prefix)
if (length(ip_bams)    == 0) stop("No IP BAMs found matching: ",    ip_prefix)

Input.file <- paste(input_bams, collapse = " ")
IP.file    <- paste(ip_bams,    collapse = " ")

message("Input BAMs : ", Input.file)
message("IP BAMs    : ", IP.file)

# ─── Output directory ────────────────────────────────────────
BamDir <- getwd()
OutDir <- file.path(BamDir, ip_prefix)
dir.create(OutDir, showWarnings = FALSE)

# ─── MACS2 peak calling ──────────────────────────────────────
macs2_cmd <- paste0(
    "macs2 callpeak",
    " -t ", IP.file,
    " -c ", Input.file,
    " -n ", ip_prefix,
    " -g ", GENOME_SIZE,
    " -f BAM",
    " --nomodel --extsize 100",
    " --outdir ", OutDir
)
message("Running MACS2:\n  ", macs2_cmd)
system(macs2_cmd)

# ─── Process summit BED → ±100 bp windows ────────────────────
setwd(OutDir)
summit_file <- paste0(ip_prefix, "_summits.bed")
if (!file.exists(summit_file)) stop("Summit file not found: ", summit_file)

peaks <- read.delim(summit_file, header = FALSE)
colnames(peaks)[1:2] <- c("chr", "summit")

peaks %>%
    dplyr::select(chr, summit) %>%
    group_by(chr, summit) %>%
    mutate(summit = mean(
        as.numeric(unlist(strsplit(as.character(summit), "_")))
    )) %>%
    ungroup() %>%
    mutate(start = summit - 100, end = summit + 100) %>%
    dplyr::select(chr, start, end) %>%
    write.table("summit_200.bed",
                sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# ─── Peak annotation ─────────────────────────────────────────
message("Annotating peaks...")
peakAnno <- annotatePeak(
    summit_file,
    TxDb    = txdb,
    tssRegion = c(0, 0),
    genomicAnnotationPriority = c("5UTR", "3UTR", "Exon", "Intron",
                                  "Downstream", "Intergenic", "Promoter")
)
write.csv(as.data.frame(peakAnno@anno), "peaks.annotation.csv", row.names = FALSE)

# ─── Plot 1: gene body profile (plotPeakProf2) ───────────────
pdf("peaks.distribution.pdf")
plotPeakProf2(
    summit_file,
    TxDb       = txdb,
    upstream   = rel(0.2),
    downstream = rel(0.2),
    conf       = 0.95,
    by         = "gene",
    type       = "body",
    nbin       = 20
)
dev.off()

# ─── Plot 2: CDS-centred profile with flanks ─────────────────

# Theme helpers
BoldTheme <- function() {
    theme(
        axis.text  = element_text(color = 1, size = 12, face = "bold"),
        plot.title = element_text(hjust = 0.5, color = 1, size = 15, face = "bold"),
        axis.title = element_text(color = 1, size = 14, face = "bold")
    )
}

# Build merged CDS GRanges (gene-level extent)
allcds <- cdsBy(txdb, by = "gene")
allcds <- unlist(allcds)
allcds_df <- data.frame(allcds, gene = names(allcds)) %>%
    group_by(gene) %>%
    mutate(start = min(start), end = max(end)) %>%
    dplyr::select(seqnames, start, end, gene, strand) %>%
    distinct()
allcds_gr <- GRanges(
    seqnames = allcds_df$seqnames,
    ranges   = IRanges(start = allcds_df$start, end = allcds_df$end),
    strand   = allcds_df$strand
)

upstream   <- 2000
downstream <- 2000
bins       <- 800

# x-axis break positions corresponding to: 5UTR | CDS start | CDS end | 3UTR
breaks <- c(
    1,
    bins * 0.1 * upstream / 1000,
    bins + bins * 0.1 * upstream / 1000,
    bins + 2 * bins * 0.1 * upstream / 1000
)

getTagMatrix2 <- function(peak_file, upstream, downstream, bins) {
    getTagMatrix(
        peak_file,
        windows    = makeBioRegionFromGranges(allcds_gr, type = "body", by = "gene"),
        nbin       = bins,
        upstream   = upstream,
        downstream = downstream
    )
}

tagMatrixList <- lapply(list(summit_file), getTagMatrix2,
                        upstream = upstream, downstream = downstream, bins = bins)

p <- plotPeakProf(tagMatrixList) +
    scale_x_continuous(
        breaks = breaks,
        labels = c("5UTR", "CDS_start", "CDS_end", "3UTR")
    )

pdf("peaks.distribution2.pdf", width = 8, height = 6)
print(p)
dev.off()

message("Done. Results written to: ", OutDir)