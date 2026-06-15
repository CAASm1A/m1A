#!/usr/bin/env Rscript
# ============================================================
# 09_functional_analysis.R — Functional analysis of m1A sites
#
# Analyses:
#   1. RNA–protein Spearman correlation (m1A vs non-m1A genes)
#   2. Protein-to-RNA ratio (PTR) comparison
#   3. Cumulative expression (ECDF) and violin plots
#   4. m1A fraction across expression quantile bins
#
# Requires:
#   proj/prepared_data.Rdata    (from 01_prepare_data.R)
#   proj/all_m1A.info.Rdata     (from 02_annotate_peaks.R)
#   proj/allTPM.csv             (from 01_prepare_data.R or featureCounts)
#
# Outputs (in proj/):
#   protein_RNA.spearman.correlation*.pdf / *.csv
#   PTR.boxplot.pdf, PTR.boxplot.stat.csv
#   m1A.expression.ecdf.pdf, m1A.expression.violin.pdf
#   m1A.expression.violin.stat.csv
#   methy_RNAexp2.pdf, methy_RNAexp_allsamples.pdf
#   methy_RNAexp_bins.csv
# ============================================================

suppressPackageStartupMessages({
    library(dplyr)
    library(reshape2)
    library(ggplot2)
    library(ggpubr)
    library(impute)
})

source("00_config.R")
setwd("proj")
load("prepared_data.Rdata")   # sampleinfo, files, files.color
load("all_m1A.info.Rdata")    # all_m1A.info

# ─── Input files ────────────────────────────────────────────
FEAT_COUNT_RES <- "/path/to/all.featureCount.res"
PROTEIN_FILE   <- "/path/to/protein.txt"      # Protein abundance (will be KNN-imputed)
PROTEIN2_FILE  <- "/path/to/protein2.txt"     # Additional protein data (no imputation)

# ─── Sample groups ───────────────────────────────────────────
SS1 <- c("CK", "IAA", "ABA", "ACC", "GA3", "6BA")
# SS2 is derived after protein data is loaded (samples in protein but not in SS1)

# ─── Helper functions ────────────────────────────────────────
`%||%` <- function(a, b) if (!is.null(a)) a else b

name_to_id <- function(sample_name) {
    sampleinfo |> dplyr::filter(ID_name %in% sample_name) |> pull(ID)
}

m1A_gene_ids <- function(sample_name, annotation_pattern = NULL) {
    dat <- all_m1A.info[[name_to_id(sample_name)]]
    if (!is.null(annotation_pattern)) {
        dat <- dat |> dplyr::filter(grepl(annotation_pattern, annotation))
    }
    unique(dat$geneId)
}

# ============================================================
# Compute TPM from featureCounts output
# ============================================================
message("Computing TPM...")
fc         <- read.table(FEAT_COUNT_RES, header = TRUE, check.names = FALSE)
counts_mat <- as.matrix(fc[, 7:ncol(fc)])
rownames(counts_mat) <- fc[, 1]
colnames(counts_mat) <- dirname(gsub("matching/", "", colnames(counts_mat)))

tpm_fn <- function(counts, lengths) {
    rpk <- counts / lengths
    t(t(rpk) * 1e6 / colSums(rpk))
}

allTPM <- tpm_fn(counts_mat, fc$Length)
allTPM <- allTPM[, grepl("Input", colnames(allTPM))]
colnames(allTPM) <- gsub("Input|_", "", colnames(allTPM))

allTPM_m <- reshape2::melt(allTPM) |>
    mutate(ID = gsub("[12]", "", Var2)) |>
    left_join(sampleinfo, by = "ID")

allTPM_m |>
    group_by(Var1, ID_name) |>
    summarise(value = mean(value), .groups = "drop") |>
    reshape2::dcast(Var1 ~ ID_name, value.var = "value") |>
    write.csv("allTPM.csv", row.names = FALSE)

# ============================================================
# Load protein abundance data
# ============================================================
message("Loading protein data...")

load_protein <- function(path, impute_knn = TRUE) {
    dat <- read.delim(path)
    rownames(dat) <- gsub("\\..+", "", dat[, 1])
    mat <- as.matrix(dat[, -1])
    if (impute_knn) mat <- impute::impute.knn(mat)$data
    mat
}

protein_m_all <- rbind(
    reshape2::melt(load_protein(PROTEIN_FILE,  impute_knn = TRUE))  |> mutate(ID_name = gsub("[123]$", "", Var2)),
    reshape2::melt(load_protein(PROTEIN2_FILE, impute_knn = FALSE)) |> mutate(ID_name = gsub("[123]$", "", Var2))
) |> left_join(sampleinfo, by = "ID_name")

SS2 <- sampleinfo$ID_name[
    sampleinfo$ID_name %in% setdiff(unique(protein_m_all$ID_name), SS1)
] |> as.character()

# ============================================================
# 1. RNA–protein Spearman correlation
# ============================================================
message("Computing RNA-protein correlations...")

get_cor <- function(sample_name, region = NULL) {
    anno_pat <- switch(region %||% "All",
        "All"    = NULL,
        "5' UTR" = "5' UTR",
        "3' UTR" = "3' UTR",
        "CDS"    = "Exon|Intron",
        region
    )
    region_label <- if (!is.null(region) && region == "CDS") "CDS" else (region %||% "All")
    m1A_ids      <- m1A_gene_ids(sample_name, anno_pat)
    all_ids      <- unique(allTPM_m$Var1)

    pair_means <- function(gene_set) {
        left_join(
            allTPM_m    |> dplyr::filter(Var1 %in% gene_set, ID_name %in% sample_name) |>
                group_by(Var1) |> summarise(RNA     = mean(value), .groups = "drop"),
            protein_m_all |> dplyr::filter(Var1 %in% gene_set, ID_name %in% sample_name) |>
                group_by(Var1) |> summarise(protein = mean(value), .groups = "drop"),
            by = "Var1"
        )
    }

    m1A_pair     <- pair_means(m1A_ids)
    non_m1A_pair <- pair_means(setdiff(all_ids, m1A_ids))

    data.frame(
        sample     = sample_name,
        m1Cor      = cor(m1A_pair$RNA,     m1A_pair$protein,     method = "spearman", use = "complete.obs"),
        non_m1Cor  = cor(non_m1A_pair$RNA, non_m1A_pair$protein, method = "spearman", use = "complete.obs"),
        region     = region_label,
        m1.size    = nrow(m1A_pair),
        nonm1.size = nrow(non_m1A_pair)
    )
}

run_cor_set <- function(sample_set, suffix = "") {
    cor1 <- lapply(sample_set, get_cor, region = "5' UTR") |> do.call(rbind, args = _)
    cor2 <- lapply(sample_set, get_cor)                    |> do.call(rbind, args = _)
    cor3 <- lapply(sample_set, get_cor, region = "3' UTR") |> do.call(rbind, args = _)
    cor4 <- lapply(sample_set, get_cor, region = "CDS")    |> do.call(rbind, args = _)
    all_cor <- rbind(cor1, cor2, cor3, cor4)

    p1 <- ggplot(
        reshape2::melt(all_cor) |> mutate(sample = factor(sample, levels = sample_set)),
        aes(variable, value, color = sample, group = sample)
    ) +
        geom_point() + geom_line(linetype = "dashed") +
        facet_wrap(~region, scales = "free_y") +
        theme_bw() + BoldTheme() + labs(x = "", y = "Spearman correlation")
    ggsave(p1, filename = paste0("protein_RNA.spearman.correlation", suffix, ".pdf"),
           width = 12, height = 8)

    p2 <- ggplot(
        reshape2::melt(all_cor) |> mutate(sample = factor(sample, levels = sample_set)),
        aes(sample, value, color = variable, shape = region)
    ) +
        geom_point(size = 3) + labs(x = "", y = "Spearman correlation") +
        theme_bw() + BoldTheme()
    ggsave(p2, filename = paste0("protein_RNA.spearman.correlation", suffix, ".points.pdf"),
           width = 7, height = 5)

    write.csv(cor2,               paste0("correlation_all",        suffix, ".csv"), row.names = FALSE)
    write.csv(rbind(cor1,cor3,cor4), paste0("correlation_UTR_CDS", suffix, ".csv"), row.names = FALSE)
    write.csv(all_cor,            paste0("correlation_allregions", suffix, ".csv"), row.names = FALSE)
}

run_cor_set(SS1, suffix = "")
run_cor_set(SS2, suffix = "2")

# ============================================================
# 2. Protein-to-RNA ratio (PTR)
# ============================================================
message("Computing PTR...")

get_PTR <- function(sample_name) {
    m1A_ids <- m1A_gene_ids(sample_name, "5' UTR")

    compute_ptr <- function(gene_set, type_label) {
        paired <- left_join(
            allTPM_m      |> dplyr::filter(Var1 %in% gene_set, ID_name %in% sample_name) |>
                group_by(Var1) |> summarise(RNA     = mean(value), .groups = "drop"),
            protein_m_all |> dplyr::filter(Var1 %in% gene_set, ID_name %in% sample_name) |>
                group_by(Var1) |> summarise(protein = mean(value), .groups = "drop"),
            by = "Var1"
        ) |> na.omit()
        data.frame(
            sample = sample_name, type = type_label, gene = paired$Var1,
            value  = log2(paired$protein + 1) - log2(paired$RNA + 1)
        )
    }

    rbind(
        compute_ptr(m1A_ids,                                  "m1PTR"),
        compute_ptr(setdiff(unique(allTPM_m$Var1), m1A_ids), "non_m1PTR")
    )
}

allPTR1 <- lapply(SS1, get_PTR) |> do.call(rbind, args = _)
allPTR2 <- lapply(SS2, get_PTR) |> do.call(rbind, args = _)

make_ptr_plot <- function(ptr_data, sample_levels, ylim = NULL) {
    p <- ptr_data |>
        mutate(sample = factor(sample, levels = sample_levels)) |>
        ggplot(aes(sample, value, fill = type)) +
        geom_boxplot(outlier.shape = NA) +
        theme_classic() + BoldTheme() +
        stat_compare_means(label = "p.signif") +
        labs(x = "", y = "PTR (log2 protein / RNA)")
    if (!is.null(ylim)) p <- p + coord_cartesian(ylim = ylim)
    p
}

write.csv(
    rbind(
        compare_means(data = allPTR1, value ~ type, group.by = "sample"),
        compare_means(data = allPTR2, value ~ type, group.by = "sample")
    ),
    "PTR.boxplot.stat.csv", row.names = FALSE
)

pdf("PTR.boxplot.pdf", 12, 5)
print(make_ptr_plot(allPTR1, SS1))
print(make_ptr_plot(allPTR2, SS2, ylim = c(-5, 5)))
dev.off()

# ============================================================
# 3. Cumulative expression (ECDF) and violin plots
# ============================================================
message("Plotting expression distributions...")

all_exp <- lapply(unique(as.character(sampleinfo$ID_name)), function(sample_name) {
    m1A_ids <- m1A_gene_ids(sample_name, "5' UTR")
    rbind(
        allTPM_m |> dplyr::filter( Var1 %in%  m1A_ids, ID_name %in% sample_name) |>
            transmute(sample = sample_name, type = "m1A",     value),
        allTPM_m |> dplyr::filter(!Var1 %in%  m1A_ids, ID_name %in% sample_name) |>
            transmute(sample = sample_name, type = "non-m1A", value)
    )
}) |> do.call(rbind, args = _)

ggsave(
    ggplot(all_exp, aes(log2(value), color = type)) +
        stat_ecdf() + facet_wrap(~sample) +
        labs(x = "log2(TPM)", y = "Cumulative fraction") +
        theme_classic() + BoldTheme(),
    filename = "m1A.expression.ecdf.pdf", width = 10, height = 6
)

ggsave(
    ggplot(all_exp, aes(type, log2(value), fill = type)) +
        geom_violin() + facet_wrap(~sample) +
        labs(y = "log2(TPM)", x = "") +
        theme_classic() + BoldTheme(),
    filename = "m1A.expression.violin.pdf", width = 10, height = 6
)

write.csv(
    compare_means(data = all_exp, value ~ type, group.by = "sample"),
    "m1A.expression.violin.stat.csv", row.names = FALSE
)

# ============================================================
# 4. m1A fraction across expression quantile bins
# ============================================================
message("Computing m1A fraction per expression quintile bin...")

rna_wide <- read.csv("allTPM.csv")

compute_bin_data <- function(sample_name) {
    m1A_anno <- all_m1A.info[[name_to_id(sample_name)]] |>
        dplyr::select(geneId, annotation) |> data.frame()

    sel <- data.frame(geneId = rna_wide$Var1,
                      value  = rna_wide[[as.character(sample_name)]]) |>
        left_join(m1A_anno, by = "geneId")

    bks    <- quantile(sel$value, probs = seq(0, 1, length.out = 6))
    bks    <- bks + seq_along(bks) * 1e-12
    bks[1] <- 0

    sel |>
        mutate(
            bins = cut(value, breaks = bks, labels = 1:5, include.lowest = TRUE),
            m1A  = ifelse(is.na(annotation), "no-m1A", "m1A")
        ) |>
        group_by(m1A, bins) |> summarise(n = n(), .groups = "drop") |>
        group_by(bins)       |> mutate(percent = 100 * n / sum(n)) |>
        dplyr::filter(m1A == "m1A") |>
        mutate(sample = sample_name)
}

all_samples <- unique(as.character(sampleinfo$ID_name))

pdf("methy_RNAexp2.pdf")
for (sample_name in all_samples) {
    print(
        ggplot(compute_bin_data(sample_name), aes(bins, percent)) +
            geom_bar(stat = "identity", width = 0.5) +
            theme_classic() + BoldTheme() +
            labs(x = "Expression quintile", y = "% m1A genes per bin") +
            ggtitle(sample_name)
    )
}
dev.off()

all_bin_data <- lapply(all_samples, compute_bin_data) |> do.call(rbind, args = _)

ggsave(
    ggplot(all_bin_data, aes(bins, percent)) +
        geom_bar(stat = "identity", width = 0.5) +
        facet_wrap(~sample) + theme_classic() + BoldTheme() +
        labs(x = "Expression quintile", y = "% m1A genes per bin"),
    filename = "methy_RNAexp_allsamples.pdf", width = 12, height = 8
)

write.csv(all_bin_data, "methy_RNAexp_bins.csv", row.names = FALSE)

message("Done: 09_functional_analysis.R")