#!/usr/bin/env Rscript
# ============================================================
# 06_site_statistics.R — m1A site statistics & overlap analysis
#
# Requires: proj1/prepared_data.Rdata
#           proj1/all_m1A.info.Rdata
#
# Outputs (in proj1/):
#   intensity.pdf           — IP signal intensity violin plot
#   stat_m1A.sites.pdf      — m1A sites per gene (bar chart)
#   m1A.upset.pdf           — UpSet plot of site overlaps across samples
#   exp1.pdf                — TPM expression by annotation region (sample-unique sites)
#   <sample>.genes.m1A.number.csv  — per-gene m1A site counts
# ============================================================

suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(reshape2)
    library(UpSetR)
})

source("00_config.R")
setwd("proj")
load("prepared_data.Rdata")
load("all_m1A.info.Rdata")

# ─── 1. IP signal intensity ───────────────────────────────────
message("Plotting IP intensity...")
inten_raw <- read.csv(INTENSITY, row.names = 1) |> as.matrix()
inten_m   <- inten_raw[, IDS]
colnames(inten_m) <- sampleinfo$ID_name[match(IDS, sampleinfo$ID)]

p_inten <- inten_m |>
    reshape2::melt() |>
    ggplot(aes(Var2, log2(value), fill = Var2)) +
    geom_violin() +
    scale_fill_manual(values = files.color) +
    labs(x = "", y = "log2(intensity)") +
    theme_classic() +
    BoldTheme() +
    NoLegend()
ggsave(p_inten, filename = "intensity.pdf", width = 8, height = 5)

# ─── 2. m1A sites per gene (bar chart) ───────────────────────
message("Computing m1A sites per gene...")
allSites <- lapply(names(all_m1A.info), function(sample_name) {
    all_m1A.info[[sample_name]]$anno@anno |>
        data.frame() |>
        group_by(geneId) |>
        summarise(n = n(), .groups = "drop") |>
        mutate(n = ifelse(n >= 4, "4+", as.character(n))) |>
        group_by(n) |>
        summarise(num = n(), .groups = "drop") |>
        mutate(percent = 100 * num / sum(num), sample = sample_name)
}) |> do.call(rbind, args = _)

m_sites <- allSites |>
    reshape2::acast(n ~ sample, value.var = "percent") |>
    replace(is.na(x = _), 0) |>
    reshape2::melt()

p_sites <- ggplot(
    m_sites |> mutate(sample = factor(Var2, levels = names(files.color))),
    aes(sample, value, fill = Var1)
) +
    geom_bar(stat = "identity", position = position_dodge(), width = 0.6) +
    theme_classic() +
    labs(x = "", y = "Percent of m1A genes", fill = "Sites per gene")
ggsave(p_sites, filename = "stat_m1A.sites.pdf", width = 10, height = 8)

# ─── 3. Per-sample site counts CSV ───────────────────────────
message("Writing per-gene m1A counts...")
lapply(names(all_m1A.info), function(sample_name) {
    all_m1A.info[[sample_name]]$anno@anno |>
        data.frame() |>
        group_by(geneId) |>
        summarise(m1A_sites = n(), .groups = "drop") |>
        write.csv(paste0(sample_name, ".genes.m1A.number.csv"), row.names = FALSE)
})

# ─── 4. UpSet plot ────────────────────────────────────────────
message("Plotting UpSet diagram...")
allsites <- lapply(files, function(bed_path) {
    a <- read.delim(bed_path, header = FALSE)
    paste(a[, 1], a[, 2], a[, 3], sep = "-")
})
names(allsites) <- names(files)

pdf("m1A.upset.pdf", 9, 7)
upset(fromList(allsites), nsets = length(files))
dev.off()

# ─── 5. Sample-unique sites & expression analysis ─────────────
message("Analysing sample-unique m1A sites vs expression...")

sample_unique <- lapply(seq_along(allsites), function(i) {
    others <- unique(unlist(allsites[-i]))
    setdiff(allsites[[i]], others)
})
names(sample_unique) <- names(allsites)

# Load all-sample TPM and reshape
allTPM_raw <- read.csv(TRANS_TPM, row.names = 1) |> as.matrix()
allTPM_m   <- reshape2::melt(allTPM_raw) |>
    mutate(sample = gsub("[0-9]", "", Var2)) |>
    group_by(Var1, sample) |>
    summarise(value = mean(value, na.rm = TRUE), .groups = "drop") |>
    dplyr::filter(grepl("Input", sample)) |>
    mutate(sample = gsub("_|Input", "", sample)) |>
    dplyr::filter(sample %in% sampleinfo$ID[seq_along(IDS)])

allTPM_m$sample <- sampleinfo$ID_name[match(allTPM_m$sample, sampleinfo$ID)]

# Build unique-site expression table with annotation
all_anno_regions <- c("5' UTR", "3' UTR", "Exon", "Intron", "Intergenic")

xxx_data <- lapply(names(all_m1A.info), function(sample_name) {
    all_m1A.info[[sample_name]]$anno@anno |>
        data.frame() |>
        mutate(
            ID         = paste(seqnames, start - 1, end, sep = "-"),
            annotation = gsub(" \\(.*", "", annotation)
        ) |>
        dplyr::filter(ID %in% sample_unique[[sample_name]]) |>
        dplyr::select(ID, annotation, geneId) |>
        left_join(allTPM_m, by = c("geneId" = "Var1")) |>
        mutate(label = sample_name)
}) |> do.call(rbind, args = _)

make_region_plot <- function(data, region) {
    ggplot(
        data |> dplyr::filter(annotation == region),
        aes(x = sample, y = log(value + 1), color = label)
    ) +
        geom_boxplot() +
        facet_wrap(~label, scales = "free_y") +
        theme_classic() +
        BoldTheme() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(x = region, y = "Mean TPM (log)", color = "Sample")
}

pdf("exp1.pdf", 12, 5)
for (region in all_anno_regions) {
    print(make_region_plot(xxx_data, region))
}
dev.off()

message("Done: 06_site_statistics.R")