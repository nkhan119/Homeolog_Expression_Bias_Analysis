# 🧬 Homeolog Expression Bias Analysis in Allopolyploids

> Quantify and visualise subgenome expression dominance across tissues in allopolyploid species — comparing homeologous gene copies between the **H**, **St1**, and **St2** subgenomes.

---

## ✨ Features

- Loads pre-averaged TPM expression data for four tissues (Root, Stem, Leaf, Flower)
- Computes pairwise **log₂ expression ratios** for all three subgenome comparisons: H vs St1 · H vs St2 · St1 vs St2
- Classifies each gene triplet as **H-dominant**, **St1-dominant**, **St2-dominant**, or **Balanced**
- Generates publication-ready **histogram plots** with bias statistics and mean lines annotated per tissue
- Colour-coded gradient: 🟢 St2/St1-biased → ⚪ Unbiased → 🔵 H-biased
- Outputs both **PDF** and **JPG** at 400 dpi

---

## Repository Structure

```
homeolog-expression-bias/
├── homeolog_expression_bias.R    # Main analysis script
├── README.md

> **Data Availability** — The TSV expression files used in this study are available upon reasonable request. Please contact [nadeem.khan@inrs.ca](mailto:nadeem.khan@inrs.ca).

---

## Requirements

**R** ≥ 4.1 with the following packages:

```r
install.packages(c("dplyr", "tidyr", "ggplot2", "patchwork", "ggrepel", "viridis"))

# From CRAN or Bioconductor
install.packages("ggalluvial")
install.packages("ggtern")
```

---

## Input Format

Four tab-separated TSV files — one per tissue — each with these columns:

| Column | Description |
|--------|-------------|
| `chrH` | Gene ID on the H subgenome |
| `chrSt1` | Homeolog ID on the St1 subgenome |
| `chrSt2` | Homeolog ID on the St2 subgenome |
| `TPM1` | Expression of H (pre-averaged across replicates) |
| `TPM2` | Expression of St1 |
| `TPM3` | Expression of St2 |

Example row:
```
chrH              chrSt1            chrSt2            TPM1   TPM2   TPM3
EhH_g00001.1      EhSt1_g00001.1    EhSt2_g00001.1    12.4   4.1    3.9
```

---

## Usage

**1. Edit the CONFIG section** at the top of the script to point to your input files:

```r
INPUT_FILES <- list(
  Root   = "Root_homologous_expression.tsv",
  Stem   = "Stem_homologous_expression.tsv",
  Leaf   = "Leaf_homologous_expression.tsv",
  Flower = "Flower_homologous_expression.tsv"
)
OUTPUT_DIR <- "output"
```

**2. Run the script:**

```bash
Rscript homeolog_expression_bias.R
```

Or from within R:

```r
source("homeolog_expression_bias.R")
```

---

## Output

All files are saved to the `output/` directory:

| File | Description |
|------|-------------|
| `histogram_H_vs_St1_across_tissues.pdf/.jpg` | Log₂(H/St1) per tissue |
| `histogram_H_vs_St2_across_tissues.pdf/.jpg` | Log₂(H/St2) per tissue |
| `histogram_St1_vs_St2_across_tissues.pdf/.jpg` | Log₂(St1/St2) per tissue |

Each plot shows:
- Histogram of log₂ ratios coloured by bias direction
- Dashed red lines at ±1 (2-fold threshold)
- Annotated mean lines with gene counts per bias category

---

## Methods

**Dominance classification** (fold-change threshold = 2×, i.e. |log₂ ratio| > 1):

| Class | Criterion |
|-------|-----------|
| `H_dominant` | H > St1 and H > St2 by ≥ 2× |
| `St1_dominant` | St1 > H and St1 > St2 by ≥ 2× |
| `St2_dominant` | St2 > H and St2 > St1 by ≥ 2× |
| `Balanced` | No subgenome exceeds 2× over others |
| `LowExpr` | All subgenomes below 25th percentile of non-zero TPM |

A **pseudocount** (1/10 of the minimum non-zero TPM) is added before log transformation to handle zero-expression values.

---

## Citation

If you use this script in your research, please cite:

> Ji Y, Chaudhary R, Khan N, Perumal S, Wang Z, Moghaddam Moghanloo L, Hucl P, Biligetu B, Sharpe AG, Jin L (in
submission) Haplotype-resolved genome assemblies of hybrid wheatgrass and bluebunch wheatgrass reveal the
stepwise polyploid origin and biased subgenome dominance.

---

## Contact

**Nadeem Khan, PhD**  
Bioinformatician — INRS–Centre Armand-Frappier Santé-Biotechnologie, Laval, QC, Canada  
nkhan119@uottawa.ca  
[@nkhan119](https://github.com/nkhan119)
