# N1-methyladenosine mRNA Methylation in *Arabidopsis*

## Abstract

N1-methyladenosine (m1A) is a recently identified mRNA methylation mark that affects gene expression and translation. In the model plant *Arabidopsis thaliana*, its landscape and function remain uncharacterized. Here, we present base-resolution m1A methylomes across diverse *Arabidopsis* tissues and show that m1A is enriched in 5′ UTR and preferentially localizes to highly expressed genes, negatively correlating with mRNA translation.

---

## Table of Contents

- [Pipeline Overview](#pipeline-overview)
- [Repository Structure](#repository-structure)
- [Environment Setup](#environment-setup)
- [Installing Required Libraries](#installing-required-libraries)
- [Clone the Repository](#clone-the-repository)
- [Usage](#usage)
- [File Descriptions](#file-descriptions)
- [Questions and Errors](#questions-and-errors)
- [Contact](#contact)

---

## Repository Structure

```
.
├── 1.align.sh                   # Step 1: QC and alignment
├── 2.callpeak.sh                # Step 2: Peak calling
├── 3.m1Asites.analysis.sh       # Step 3: m1A site identification
├── 4.run_m1A.sh                 # Step 4: Downstream analysis launcher
└── Rscript/
    ├── 00.config.R              # Shared paths, parameters, and theme functions
    ├── 01.prepare.R             # Load genome, TxDb, expression matrices
    ├── 02.annotate_peaks.R      # Annotate m1A sites; compute AUG distances
    ├── 03.genomic_distribution.R# mRNA region distribution and landmark plots
    ├── 04.motif.R               # Sequence motif discovery (MEME)
    ├── 05.GC_MFE.R              # 5′ UTR GC content and minimum free energy
    ├── 06.site_statistics.R     # Site counts, intensity, UpSet, sample-unique sites
    ├── 07.GO_enrichment.R       # GO enrichment and expression heatmap
    ├── 08.PCA.R                 # PCA of IP samples on peak coverage
    ├── 09.functional.R          # RNA–protein correlation, PTR, expression analysis
    ├── callpeak.R               # Per-sample MACS2 wrapper called by 2.callpeak.sh
    └── MFE_UTR.py               # RNALfold sliding-window MFE computation
```

---

## Environment Setup

The pipeline requires the following software. Versions used in development are listed for reproducibility.

| Tool | Version | Purpose |
|------|---------|---------|
| [Trim Galore](https://github.com/FelixKrueger/TrimGalore) | ≥ 0.6 | Read trimming and QC |
| [STAR](https://github.com/alexdobin/STAR) | 2.7.2a | RNA-seq alignment |
| [samtools](http://www.htslib.org/) | ≥ 1.15 | BAM processing |
| [deepTools](https://deeptools.readthedocs.io/) | ≥ 3.5 | bigwig generation and matrix computation |
| [Subread / featureCounts](https://subread.sourceforge.net/) | ≥ 2.0 | Gene-level read quantification |
| [MACS2](https://github.com/macs3-project/MACS) | 2.2.7.1 | Peak calling |
| [sequenza-utils](https://sequenzatools.bitbucket.io/) | ≥ 3.0 | pileup2acgt conversion |
| [MEME Suite](https://meme-suite.org/) | 5.5.2 | Motif discovery |
| [RNAfold / RNALfold](https://www.tbi.univie.ac.at/RNA/) | 2.5.1 | RNA secondary structure MFE |
| R | 4.3.3 | Statistical analysis and visualization |
| Python | 3.12.7 | MFE computation helper script |

A conda environment is recommended:

```bash
conda create -n m1A_pipeline python=3.10
conda activate m1A_pipeline

conda install -c bioconda trim-galore star samtools deeptools subread macs2 sequenza-utils
conda install -c bioconda -c conda-forge viennarna
conda install -c bioconda meme
```

---

## Installing Required Libraries

**R packages:**

```r
# CRAN
install.packages(c("dplyr", "ggplot2", "reshape2", "ggpubr",
                   "pheatmap", "UpSetR", "impute"))

# Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
    "GenomicFeatures", "GenomicRanges", "Biostrings",
    "BSgenome", "ChIPseeker", "clusterProfiler",
    "org.At.tair.db", "IRanges"
))
```

**Python packages** (for `MFE_UTR.py`):

```bash
pip install biopython
```

---

## Clone the Repository

```bash
git clone https://github.com/<your-org>/<your-repo>.git
cd <your-repo>
```

Before running, open `Rscript/00.config.R` and update **all** `/path/to/...` placeholders to match your local file system:

```r
GENOME_FA    <- "/path/to/TAIR10/tair10.noChr.fa"
GTF_PATH     <- "/path/to/TAIR10/genes/genes.gtf"
GENOME_REF   <- "/path/to/TAIR10/star"
...
```

---

## File Descriptions

### Shell scripts

| File | Description |
|------|-------------|
| `1.align.sh` | Trims reads with Trim Galore, aligns with STAR, filters BAMs, generates normalized bigwigs with bamCoverage, and runs featureCounts for gene-level quantification. Accepts a step argument (1/2/3). |
| `2.callpeak.sh` | Calls peaks for each IP/Input sample pair in parallel using MACS2 via `callpeak.R`. |
| `3.m1Asites.analysis.sh` | Runs mpileup on peak regions, converts to per-base counts with `sequenza-utils pileup2acgt`, computes A/T mismatch ratios, and identifies high-confidence m1A sites by comparing IP vs demethylase-treated IPDe samples. |
| `4.run_m1A.sh` | Launcher for the downstream R analysis pipeline. Runs all numbered Rscripts sequentially or a single step by index. |

### R scripts

| File | Description |
|------|-------------|
| `00.config.R` | Central configuration: file paths, sample IDs, analysis parameters, and shared ggplot2 theme functions (`BoldTheme`, `NoLegend`). Source this at the top of every script. |
| `01.prepare.R` | Loads the reference genome, builds TxDb from GTF, reads sample metadata, discovers m1A site BED files, and preprocesses transcript- and gene-level TPM matrices. Saves `prepared_data.Rdata`. |
| `02.annotate_peaks.R` | Annotates m1A sites relative to genomic features using `ChIPseeker::annotatePeak`, computes distance to the nearest AUG codon, and exports per-sample annotation CSVs. Saves `all_m1A.info.Rdata`. |
| `03.genomic_distribution.R` | Plots m1A density across the mRNA architecture (5′ UTR → CDS → 3′ UTR), and distance distributions relative to AUG, stop codon, and TSS, both per sample and merged. |
| `04.motif.R` | Extracts ±N bp sequences around m1A summits and matched background sequences, then runs MEME motif discovery on each sample. |
| `05.GC_MFE.R` | Calculates GC content and RNALfold minimum free energy (MFE) in 5′ UTRs for m1A-modified vs unmodified genes, and generates per-sample positional profiles around AUG codons. |
| `06.site_statistics.R` | Plots IP signal intensity distributions, per-gene m1A site counts, UpSet diagrams of site overlaps across samples, and expression of genes carrying sample-unique sites. |
| `07.GO_enrichment.R` | Runs Gene Ontology enrichment (via `clusterProfiler`) on m1A-modified genes per sample and generates an expression heatmap. |
| `08.PCA.R` | Performs PCA on IP samples using peak coverage profiles from `*.allcoverage.bed` files. |
| `09.functional.R` | Computes RNA–protein Spearman correlations for m1A vs non-m1A gene sets across regions (5′ UTR, CDS, 3′ UTR, all), calculates PTR (protein-to-RNA ratio), plots ECDF and violin expression distributions, and quantifies m1A enrichment across expression quintile bins. |
| `callpeak.R` | Rscript wrapper called by `2.callpeak.sh`. Runs MACS2, processes summit BEDs to ±100 bp windows, and generates ChIPseeker peak annotation and distribution plots. |

### Python script

| File | Description |
|------|-------------|
| `MFE_UTR.py` | Computes sliding-window RNALfold MFE along 5′ UTR sequences. Outputs per-transcript minimum free energy values used in `05.GC_MFE.R`. |

---

## Questions and Errors

If you encounter issues:

Please use the GitHub Issues page and include the error message, the command or script name, and your R/tool version output (`sessionInfo()` for R errors).

---

## Contact

For questions about the analysis pipeline or data, please contact:

liangzhe
[liangzhe@caas.cn]

For questions related to the manuscript, please refer to the corresponding author listed in the paper.