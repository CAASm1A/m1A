# ============================================================
# 00_config.R — Shared configuration for m1A-seq analysis
#
# Source this file at the top of each analysis script:
#   source("00_config.R")
# ============================================================

# ─── Paths ──────────────────────────────────────────────────
GENOME_FA    <- "/path/to/TAIR10/tair10.noChr.fa"
GTF_PATH     <- "/path/to/TAIR10/genes/genes.gtf"
SITE_BED_DIR <- "/path/to/analysis/merge_dup"          # Directory containing *.filter.site.bed
BAMS_DIR     <- "/path/to/bams"                         # Directory containing *.allcoverage.bed

# External function library (dis_mRNA, plot_aug, get_fa, get_motif, etc.)
PLOT_FNS     <- "/path/to/plot_fns.R"

# Pre-computed files
MFE_RES      <- "/path/to/all5UTR.MFE.res"             # RNALfold output
SAMPLE_INFO  <- "/path/to/sampleinfo.txt"               # Tab-separated sample metadata
TRANS_TPM    <- "/path/to/trans.TPM.csv"                # Transcript-level TPM matrix
GENE_TPM     <- "/path/to/genes.TPM.csv"                # Gene-level TPM matrix
INTENSITY    <- "/path/to/bams/merge_peaks/allintensity.csv"

# ─── Sample IDs ─────────────────────────────────────────────
# IDs must match the ID column in SAMPLE_INFO
IDS <- c("A", "H", "E", "F", "G", "C", "D", "P", "I", "J")

# ─── Analysis parameters ────────────────────────────────────
N_CORES      <- 30      # Cores for parallel::mclapply
AUG_EXTEND   <- 600     # Flanking region around AUG codon (bp)
STOP_EXTEND  <- 2000    # Flanking region around stop codon (bp)
TSS_EXTEND   <- 2000    # Flanking region around TSS (bp)
MOTIF_EXT    <- 10      # Extension around m1A site for motif extraction (bp)
BG_SAMPLE_N  <- 2000    # Number of non-m1A AUG sites to sample as background

# ─── ggplot2 theme helpers ───────────────────────────────────
library(ggplot2)

BoldTheme <- function() {
    theme(
        axis.text  = element_text(color = 1, size = 12, face = "bold"),
        plot.title = element_text(hjust = 0.5, color = 1, size = 15, face = "bold"),
        axis.title = element_text(color = 1, size = 14, face = "bold")
    )
}

NoLegend <- function() {
    theme(legend.position = "none", validate = TRUE)
}

#!/usr/bin/env Rscript
# ============================================================
# plot_fns.R — Helper functions for m1A-seq downstream analysis
#
# Source this file after 00_config.R:
#   source("plot_fns.R")
#
# Functions:
#   anno_sites()     — annotate a character vector of "chr_start_end" sites
#   build_gr()       — build region GRanges (5UTR/CDS/3UTR/flanks) from TxDb
#   get_dis()        — compute positional coordinates of peaks within regions
#   dis_mRNA()       — plot m1A density across mRNA architecture
#   get_fa()         — extract foreground & background FASTA sequences
#   get_motif()      — run MEME motif discovery
#   get_randSeq()    — sample random background sequences from genome
#   plot_aug()       — density plot of m1A sites relative to a landmark
#   gc_content()     — per-position GC% across a DNAStringSet
#   get_EV_steps()   — sliding-window RNAfold MFE along sequences
# ============================================================

suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(GenomicRanges)
    library(GenomicFeatures)
    library(ChIPseeker)
    library(Biostrings)
    library(BSgenome)
})

# ─── User-configurable tool paths ───────────────────────────
MEME_BIN   <- "/path/to/meme/bin/meme"       # MEME executable
RNAFOLD    <- "/path/to/bin/RNAfold"          # RNAfold executable
N_CORES    <- 30                              # Cores for parallel::mclapply

# ============================================================
# anno_sites — annotate "chr_start_end" formatted site IDs
# ============================================================
anno_sites <- function(sites, txdb) {
    suppressPackageStartupMessages({
        require(clusterProfiler)
        require(ChIPseeker)
        require(GenomicRanges)
    })
    coords <- strsplit(sites, "_") |> do.call(rbind, args = _)
    peaks  <- GRanges(
        seqnames = coords[, 1],
        ranges   = IRanges(as.numeric(coords[, 2]), as.numeric(coords[, 3]))
    )
    annotatePeak(
        peaks,
        TxDb    = txdb,
        tssRegion = c(0, 0),
        genomicAnnotationPriority = c("5UTR", "3UTR", "Exon", "Intron",
                                      "Downstream", "Intergenic", "Promoter")
    )
}

# ============================================================
# build_gr — assemble per-transcript region GRanges
#
# Returns a GRanges with columns: rank, txname, type, index, length
# type ∈ {"5UTR", "CDS", "3UTR", "pre_5UTR", "pre_3UTR"}
# ============================================================
build_gr <- function(txdb, downstream = 100, ext = 100, cds = 200, utr = 100) {
    alltrans <- transcriptsBy(txdb) |> unlist()
    cdsall   <- cdsBy(txdb, by = "tx", use.name = TRUE) |> unlist()
    all5utr  <- fiveUTRsByTranscript(txdb, use.names = TRUE) |> unlist()
    all3utr  <- threeUTRsByTranscript(txdb, use.names = TRUE) |> unlist()

    make_region <- function(gr, type_label) {
        data.frame(gr, txname = names(gr)) |>
            dplyr::select(seqnames, start, end, strand,
                          rank = exon_rank, txname) |>
            mutate(type = type_label) |>
            arrange(txname)
    }

    r1 <- make_region(all5utr, "5UTR")
    r2 <- make_region(all3utr, "3UTR")
    r3 <- make_region(cdsall,  "CDS")

    r4 <- data.frame(alltrans, txname = alltrans$tx_name) |>
        mutate(
            start1 = ifelse(strand == "+", start - downstream, end + 1),
            end1   = ifelse(strand == "+", start - 1,          end + downstream)
        ) |>
        dplyr::select(seqnames, start = start1, end = end1, strand,
                      rank = tx_id, txname) |>
        mutate(type = "pre_5UTR")

    r5 <- data.frame(alltrans, txname = alltrans$tx_name) |>
        mutate(
            start1 = ifelse(strand == "+", end + 1,          start - downstream),
            end1   = ifelse(strand == "+", end + downstream,  start - 1)
        ) |>
        dplyr::select(seqnames, start = start1, end = end1, strand,
                      rank = tx_id, txname) |>
        mutate(type = "pre_3UTR")

    allr <- rbind(r1, r2, r3, r4, r5)

    add_index <- function(df, decreasing = FALSE) {
        df |>
            arrange(txname, if (decreasing) desc(start) else start) |>
            group_by(txname, type) |>
            mutate(index = cumsum(end - start + 1))
    }

    allregion <- rbind(
        allr |> dplyr::filter(strand == "+") |> add_index(FALSE),
        allr |> dplyr::filter(strand == "-") |> add_index(TRUE)
    ) |>
        arrange(seqnames, start) |>
        group_by(txname, type) |>
        mutate(length = sum(end - start + 1))

    GRanges(
        seqnames = allregion$seqnames,
        ranges   = IRanges(start = allregion$start, end = allregion$end),
        strand   = allregion$strand,
        rank     = allregion$rank,
        txname   = allregion$txname,
        type     = allregion$type,
        index    = allregion$index,
        length   = allregion$length
    )
}

# ============================================================
# get_dis — map peaks to positions within transcript regions
# ============================================================
get_dis <- function(peakfile, allregion.gr,
                    cus_order = c("5UTR", "CDS", "3UTR", "pre_5UTR", "pre_3UTR"),
                    ext = 100, cds = 200, utr = 100) {
    peaks <- readPeakFile(peakfile)
    over  <- findOverlaps(peaks, allregion.gr, type = "within")
    query <- as.data.frame(peaks[from(over)])
    colnames(query) <- paste0("q.", colnames(query))
    target <- as.data.frame(allregion.gr[to(over)])

    cbind(query, target) |>
        mutate(type = factor(type, levels = cus_order)) |>
        arrange(q.seqnames, q.start, type) |>
        group_by(q.seqnames, q.start) |>
        dplyr::slice(1) |>
        ungroup() |>
        mutate(
            pos1 = ifelse(
                strand == "+",
                (index - (end - q.end)) / length,
                (index - (q.end - start)) / length
            ),
            pos = case_when(
                type == "pre_5UTR" ~ ext * pos1,
                type == "5UTR"     ~ ext + utr * pos1,
                type == "CDS"      ~ ext + utr + cds * pos1,
                type == "3UTR"     ~ ext + utr + cds + utr * pos1,
                type == "pre_3UTR" ~ ext + utr + cds + utr + ext * pos1
            )
        )
}

# ============================================================
# dis_mRNA — density plot of m1A sites across mRNA architecture
# ============================================================
dis_mRNA <- function(files, txdb, label = "",
                     downstream = 100, ext = 100, utr = 100, cds = 200,
                     size = 1, bw = 0) {
    cus_order   <- c("5UTR", "CDS", "3UTR", "pre_5UTR", "pre_3UTR")
    allregion.gr <- build_gr(txdb, downstream, ext, cds, utr)

    alldata <- lapply(seq_along(files), function(n) {
        dis.data <- get_dis(files[n], allregion.gr, cus_order, ext, cds, utr)
        data.frame(dis.data, sample = label[n])
    }) |> do.call(rbind, args = _)

    p <- if (bw > 0) {
        ggplot(alldata, aes(pos, color = sample)) + geom_density(linewidth = size, bw = bw)
    } else {
        ggplot(alldata, aes(pos, color = sample)) + geom_density(linewidth = size)
    }

    ymin2 <- max(density(alldata$pos)$y) / 25
    e <- ext; u <- utr; cd <- cds

    p <- p +
        # Upstream flank
        annotate("rect",  xmin = 0,          xmax = e,             ymin = -ymin2 * 0.51, ymax = -ymin2 * 0.49, colour = "black") +
        annotate("text",  x = e * 0.5,        y = -ymin2,           label = paste0("-", downstream, "bp"), vjust = 1) +
        # 5'UTR
        annotate("rect",  xmin = e,            xmax = e + u,         ymin = -ymin2 * 0.6,  ymax = -ymin2 * 0.4,  colour = "black") +
        annotate("text",  x = e + u * 0.5,     y = -ymin2,           label = "5'UTR", vjust = 1) +
        # CDS
        annotate("rect",  xmin = e + u,        xmax = e + u + cd,    ymin = -ymin2,         ymax = 0,              colour = "grey") +
        annotate("text",  x = e + u + cd * 0.5, y = -ymin2,          label = "CDS", vjust = 1) +
        # 3'UTR
        annotate("rect",  xmin = e + u + cd,   xmax = e + u + cd + u, ymin = -ymin2 * 0.6, ymax = -ymin2 * 0.4, colour = "black") +
        annotate("text",  x = e + u + cd + u * 0.5, y = -ymin2,      label = "3'UTR", vjust = 1) +
        # Downstream flank
        annotate("rect",  xmin = e + u + cd + u, xmax = 2*e + 2*u + cd, ymin = -ymin2 * 0.51, ymax = -ymin2 * 0.49, colour = "black") +
        annotate("text",  x = e + u + cd + u + e * 0.5, y = -ymin2,  label = paste0("+", downstream, "bp"), vjust = 1) +
        theme_classic() +
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
              axis.line.x = element_blank(), legend.position = "bottom") +
        labs(x = "", y = "m1A density") +
        coord_cartesian(xlim = c(0, 2*e + 2*u + cd))

    return(p)
}

# ============================================================
# get_fa — extract foreground & chromosome-matched background FASTA
# ============================================================
get_fa <- function(site.bed.file, genome, ext = 10, fix = "center") {
    genome_gr <- GRanges(
        seqnames = names(genome),
        ranges   = IRanges(start = 1, end = width(genome))
    )

    sites_gr <- if (inherits(site.bed.file, "GRanges")) {
        site.bed.file
    } else {
        readPeakFile(site.bed.file)
    }
    sites_gr <- resize(sites_gr, width = 2 * ext, fix = fix)

    # Keep only sites fully within chromosome bounds
    over    <- findOverlaps(sites_gr, genome_gr)
    keep_idx <- which(
        seqnames(sites_gr) %in% seqnames(genome_gr) &
        start(sites_gr) > start(genome_gr[subjectHits(over)]) &
        end(sites_gr)   < end(genome_gr[subjectHits(over)])
    )
    sites_gr <- sites_gr[keep_idx]

    sites_seq  <- BSgenome::getSeq(genome, sites_gr)
    sites_fa   <- paste0(">", seq_along(sites_seq), "\n", as.character(sites_seq))

    # Background: chromosome-matched random positions
    bg_counts  <- table(seqnames(sites_gr))
    bg_counts  <- bg_counts[bg_counts > 0]
    genome_df  <- data.frame(names = names(genome), width = width(genome))

    bg_df <- lapply(names(bg_counts), function(chr) {
        n   <- bg_counts[[chr]]
        chr_len <- genome_df$width[genome_df$names == chr]
        pos <- sample(seq_len(chr_len), n)
        data.frame(seqname = chr, start = pos, end = pos + 1)
    }) |> do.call(rbind, args = _)

    bg_gr  <- resize(
        GRanges(seqnames = bg_df$seqname, ranges = IRanges(bg_df$start, bg_df$end)),
        width = 2 * ext, fix = fix
    )
    bg_seq <- BSgenome::getSeq(genome, bg_gr)
    bg_fa  <- paste0(">", seq_along(bg_seq), "\n", as.character(bg_seq))

    list(fa = sites_fa, bg = bg_fa, index = keep_idx, fas = sites_seq, bgs = bg_seq)
}

# ============================================================
# get_motif — run MEME on foreground vs background FASTA
# ============================================================
get_motif <- function(input.fa, bg.fa, label = "") {
    cmd <- paste(
        MEME_BIN, input.fa,
        "-oc", label,
        "-mod zoops",
        "-nmotifs 25",
        "-minw 8 -maxw 15",
        "-minsites 5",
        "-p 5",
        "-dna -revcomp -nostatus"
    )
    system(cmd)
}

# ============================================================
# get_randSeq — random background sequences matched to input sites
# ============================================================
get_randSeq <- function(sites, genome) {
    genome_df <- data.frame(names = names(genome), width = width(genome))
    bg_counts <- table(seqnames(sites))
    bg_counts <- bg_counts[bg_counts > 0]

    bg_df <- lapply(names(bg_counts), function(chr) {
        n       <- bg_counts[[chr]]
        chr_len <- genome_df$width[genome_df$names == chr]
        pos     <- sample(100:(chr_len - 100), n)
        data.frame(seqname = chr, start = pos, end = pos + 100)
    }) |> do.call(rbind, args = _)

    bg_gr <- GRanges(
        seqnames = bg_df$seqname,
        ranges   = IRanges(bg_df$start, bg_df$end)
    )
    BSgenome::getSeq(genome, bg_gr)
}

# ============================================================
# plot_aug — density plot of m1A sites relative to a genomic landmark
#
# peakfiles: named character vector of BED paths, OR a GRanges object
# atg:       GRanges of landmark positions (AUG, stop codon, TSS, ...)
# ============================================================
plot_aug <- function(peakfiles, atg, extend = 600, fix = "center",
                     xlab = "", colors = character(0)) {
    atg_ext <- resize(atg, width = width(atg) + extend, fix = fix)

    compute_pos <- function(peaks_gr, sample_name = "") {
        over   <- findOverlaps(peaks_gr, atg_ext, type = "within")
        query  <- as.data.frame(peaks_gr[from(over)])
        colnames(query) <- paste0("q.", colnames(query))
        target <- as.data.frame(atg_ext[to(over)])
        cbind(query, target) |>
            distinct(q.start, q.end, .keep_all = TRUE) |>
            mutate(
                pos    = ifelse(strand == "+",
                                q.start - (start + extend / 2),
                                (end - extend / 2) - q.start),
                sample = sample_name
            )
    }

    allpdata <- if (is.character(peakfiles)) {
        lapply(names(peakfiles), function(nm) {
            compute_pos(readPeakFile(peakfiles[[nm]]), nm)
        }) |> do.call(rbind, args = _)
    } else {
        compute_pos(peakfiles)
    }

    p <- ggplot(allpdata, aes(pos, color = sample)) +
        geom_density() +
        theme_classic() +
        labs(x = xlab, y = "m1A sites density") +
        BoldTheme() +
        theme(legend.position = "bottom")

    if (length(colors) > 1) p <- p + scale_color_manual(values = colors)
    return(p)
}

# ============================================================
# gc_content — per-position GC% across a DNAStringSet
# ============================================================
gc_content <- function(dna_set) {
    seq_matrix  <- as.matrix(dna_set)
    gc_counts   <- rowSums(seq_matrix == "G") + rowSums(seq_matrix == "C")
    100 * gc_counts / ncol(seq_matrix)
}

# ============================================================
# get_EV_steps — sliding-window RNAfold MFE along sequences
#
# Returns a matrix: rows = sequences, cols = windows
# ============================================================
get_EV_steps <- function(sequence, window_size = 30, step_size = 10) {
    parallel::mclapply(seq_along(sequence), function(i) {
        seq_i   <- as.character(sequence[i])
        seq_len <- nchar(seq_i)
        starts  <- seq(1, seq_len - window_size + 1, by = step_size)
        windows <- sapply(starts, function(s) substr(seq_i, s, s + window_size - 1))

        # Write windows to a temp file, run RNAfold, clean up
        tmp_in  <- tempfile(fileext = ".txt")
        tmp_out <- tempfile(fileext = ".res")
        writeLines(windows, tmp_in)
        system(paste(RNAFOLD, "-i", tmp_in, "--noPS >", tmp_out))
        on.exit({ unlink(tmp_in); unlink(tmp_out) }, add = TRUE)

        lines <- readLines(tmp_out)
        grep("\\(", lines, value = TRUE) |>
            gsub(pattern = ".+\\( ", replacement = "") |>
            gsub(pattern = "\\)", replacement = "") |>
            as.numeric()
    }, mc.cores = N_CORES) |> do.call(rbind, args = _)
}