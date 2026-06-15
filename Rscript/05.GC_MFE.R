#!/usr/bin/env Rscript
# ============================================================
# 05_sequence_features.R — GC content & MFE of 5' UTR
#
# Requires: proj1/prepared_data.Rdata
#           proj1/all_m1A.info.Rdata
#           MFE_RES (RNALfold output) defined in 00_config.R
#           PLOT_FNS (get_fa, gc_content, get_EV_steps) defined in 00_config.R
#
# Outputs (in proj1/):
#   GC.m1A.sites.pdf          — GC content at m1A sites vs background
#   GC.5UTR.boxplot.pdf       — 5'UTR GC: m1A vs non-m1A genes
#   GC.5UTR.boxplot.stat.csv  — statistical test results
#   EV.boxplot.pdf            — 5'UTR MFE: m1A vs non-m1A genes
#   EV.5UTR.boxplot.stat.csv
#   GC.perSample.pdf          — per-sample GC profile around AUG
#   EV.perSample.pdf          — per-sample MFE profile around AUG
# ============================================================

suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(ggpubr)
    library(Biostrings)
    library(GenomicFeatures)
    library(GenomicRanges)
    library(reshape2)
})

source("00_config.R")
source(PLOT_FNS)
setwd("proj")
load("prepared_data.Rdata")
load("all_m1A.info.Rdata")   # provides: all_m1A.info, atg, allcds

# ─── Build unified transcript annotation table ───────────────
alltrans.info <- lapply(names(all_m1A.info), function(x) {
    data.frame(all_m1A.info[[x]]$anno@anno, sample = x)
}) |>
    do.call(rbind, args = _) |>
    mutate(
        txname = transcriptId,
        type   = gsub(" \\(.+", "", annotation)
    ) |>
    dplyr::select(txname, type, sample) |>
    distinct(sample, txname, .keep_all = TRUE)

# ─── 1. GC content at m1A sites (per annotation region) ─────
message("Plotting GC content at m1A sites...")
pdf("GC.m1A.sites.pdf", 8, 5)
lapply(names(all_m1A.info), function(sample_name) {
    target    <- all_m1A.info[[sample_name]]$anno@anno
    seqs      <- get_fa(target, genome, ext = 50)
    site_seq  <- seqs$fas
    bg_seq    <- seqs$bgs
    site_type <- gsub(" \\(.+", "", target$annotation[seqs$index])

    pdata <- data.frame(
        GC     = as.numeric(letterFrequency(site_seq, letters = "GC", as.prob = TRUE)),
        GC_neg = as.numeric(letterFrequency(bg_seq,   letters = "GC", as.prob = TRUE)),
        type   = site_type
    ) |>
        dplyr::filter(type %in% c("3' UTR", "5' UTR", "Exon")) |>
        reshape2::melt()

    pdata <- rbind(pdata, mutate(pdata, type = "All")) |>
        mutate(type = factor(type, levels = c("All", "Exon", "5' UTR", "3' UTR")))

    p <- ggplot(pdata, aes(variable, value)) +
        geom_boxplot() +
        theme_classic() +
        stat_compare_means() +
        facet_wrap(~type, nrow = 1) +
        labs(x = "", y = "GC content") +
        BoldTheme() +
        ggtitle(sample_name)
    print(p)
})
dev.off()

# ─── 2. 5'UTR GC content: m1A vs non-m1A genes ───────────────
message("Computing 5'UTR GC content...")
all5utr     <- fiveUTRsByTranscript(txdb, use.names = TRUE) |> unlist()
all5utr_seq <- getSeq(genome, all5utr)

all5utr_df <- data.frame(id = names(all5utr_seq), seq = as.character(all5utr_seq)) |>
    group_by(id) |>
    summarise(seq = paste(seq, collapse = ""), .groups = "drop")

all5utr_fa <- DNAStringSet(setNames(all5utr_df$seq, all5utr_df$id))

allGC <- data.frame(
    GC     = as.numeric(letterFrequency(all5utr_fa, letters = "GC", as.prob = TRUE)),
    txname = all5utr_df$id
)

gc_alldata <- lapply(names(all_m1A.info), function(sample_name) {
    allGC |>
        left_join(alltrans.info |> dplyr::filter(sample == sample_name),
                  by = "txname") |>
        mutate(
            m1A    = ifelse(is.na(type), "non-m1A", "m1A"),
            type2  = ifelse(is.na(type), "non-m1A", type),
            sample = sample_name
        )
}) |> do.call(rbind, args = _)

gc_alldata$sample <- factor(gc_alldata$sample, levels = names(files))

p_gc <- ggplot(gc_alldata, aes(sample, 100 * GC, fill = m1A)) +
    geom_boxplot(outlier.shape = NA) +
    theme_classic() +
    stat_compare_means(aes(label = after_stat(p.signif)), label.y = 55) +
    BoldTheme() +
    coord_cartesian(ylim = c(10, 60)) +
    labs(x = "", y = "GC content in 5'UTR (%)")
ggsave(p_gc, filename = "GC.5UTR.boxplot.pdf", height = 6, width = 20)

st_gc <- compare_means(data = gc_alldata, GC ~ m1A, group.by = "sample")
write.csv(st_gc, "GC.5UTR.boxplot.stat.csv", row.names = FALSE)

# ─── 3. 5'UTR minimum free energy (MFE) ─────────────────────
message("Computing 5'UTR MFE...")

# Write 5'UTR FASTA for RNALfold (run externally if not yet done)
writeXStringSet(all5utr_fa, "all5UTR.fa")
message("  [NOTE] If MFE_RES does not exist, run externally:")
message("  RNALfold -i all5UTR.fa > ", MFE_RES)

xx     <- readLines(MFE_RES)
mfe_df <- data.frame(
    txname = gsub("^>", "", xx[grepl("^>", xx)]),
    value  = as.numeric(gsub(" \\(|\\)", "", xx[grepl("^ \\(", xx)]))
)

mfe_alldata <- lapply(names(all_m1A.info), function(sample_name) {
    mfe_df |>
        left_join(alltrans.info |> dplyr::filter(sample == sample_name),
                  by = "txname") |>
        mutate(
            m1A    = ifelse(is.na(type), "non-m1A", "m1A"),
            type2  = ifelse(is.na(type), "non-m1A", type),
            sample = sample_name
        )
}) |> do.call(rbind, args = _)

mfe_alldata$sample <- factor(mfe_alldata$sample, levels = names(files))

p_mfe <- ggplot(mfe_alldata, aes(sample, value, fill = m1A)) +
    geom_boxplot(outlier.shape = NA) +
    theme_classic() +
    stat_compare_means(aes(label = after_stat(p.signif)), label.y = 10) +
    BoldTheme() +
    coord_cartesian(ylim = c(-100, 20)) +
    labs(x = "", y = "MFE of 5' UTR (kcal/mol)")
ggsave(p_mfe, filename = "EV.boxplot.pdf", height = 6, width = 20)

st_mfe <- compare_means(data = mfe_alldata, value ~ m1A, group.by = "sample")
write.csv(st_mfe, "EV.5UTR.boxplot.stat.csv", row.names = FALSE)

# ─── 4. Per-sample GC & MFE profiles around AUG ─────────────
message("Building per-sample GC and MFE profiles around AUG...")
all_m1A_txnames    <- unique(alltrans.info$txname)
all_m1A_gr         <- atg[atg$txname %in% all_m1A_txnames]
nonm1A_gr          <- sample(atg[!(atg$txname %in% all_m1A_txnames)], BG_SAMPLE_N)

# Sequences centred on AUG ±300 bp
win_size    <- AUG_EXTEND
all_m1A_seq <- getSeq(genome, resize(all_m1A_gr,  width = width(all_m1A_gr)  + win_size, fix = "center"))
nonm1A_seq  <- getSeq(genome, resize(nonm1A_gr,   width = width(nonm1A_gr)   + win_size, fix = "center"))

all_m1A_EV  <- get_EV_steps(all_m1A_seq)
nonm1A_EV   <- get_EV_steps(nonm1A_seq)

x_labels <- as.character(seq(-win_size / 2, win_size / 2, length.out = 7))
x_breaks_gc <- seq(0, win_size, length.out = 7)
x_breaks_ev <- seq(0, ncol(all_m1A_EV), length.out = 7)

pdf("GC.perSample.pdf")
for (sample_name in names(all_m1A.info)) {
    tx_idx    <- match(
        alltrans.info |> dplyr::filter(sample == sample_name) |> pull(txname) |> unique(),
        all_m1A_gr$txname
    )
    m1A_gc    <- gc_content(all_m1A_seq[na.omit(tx_idx)])
    nonm1A_gc <- gc_content(nonm1A_seq)

    gc_data <- data.frame(m1A.GC = m1A_gc, nonm1A.GC = nonm1A_gc) |>
        as.matrix() |> reshape::melt()

    p <- ggplot(gc_data, aes(X1, value, color = X2)) +
        geom_line() +
        theme_classic() +
        BoldTheme() +
        labs(x = "Distance from AUG codon (nts)", y = "Mean GC content (%)",
             color = "Group") +
        scale_x_continuous(labels = x_labels, breaks = x_breaks_gc) +
        theme(legend.position = "bottom") +
        ggtitle(sample_name)
    print(p)
}
dev.off()

pdf("EV.perSample.pdf")
for (sample_name in names(all_m1A.info)) {
    tx_idx  <- match(
        unique(all_m1A.info[[sample_name]]$anno@anno$transcriptId),
        all_m1A_gr$txname
    )
    ev_data <- data.frame(
        m1A.EV    = colMeans(all_m1A_EV[tx_idx, ], na.rm = TRUE),
        nonm1A.EV = colMeans(nonm1A_EV,             na.rm = TRUE)
    ) |>
        as.matrix() |> reshape::melt()

    p <- ggplot(ev_data, aes(X1, value, color = X2)) +
        geom_line() +
        theme_classic() +
        BoldTheme() +
        labs(x = "Distance from AUG codon (nts)", y = "Mean EV") +
        scale_x_continuous(labels = x_labels, breaks = x_breaks_ev) +
        ggtitle(sample_name)
    print(p)
}
dev.off()

message("Done: 05_sequence_features.R")