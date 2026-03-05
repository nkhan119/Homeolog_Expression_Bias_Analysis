#!/usr/bin/env Rscript
# ═══════════════════════════════════════════════════════════════════════════════
#  homeolog_expression_bias.R
#
#  Homeolog Expression Bias Analysis in Allopolyploids
#  Subgenomes: H (hybrid) · St1 (progenitor 1) · St2 (progenitor 2)
#
#  Input  : Per-tissue TSV files with pre-averaged TPM values across replicates
#           Columns: chrH, chrSt1, chrSt2, TPM1 (=H), TPM2 (=St1), TPM3 (=St2)
#  Output : PDF/JPG histograms of log2 expression ratios across tissues
#
#  Usage  : Rscript homeolog_expression_bias.R
#           (adjust file paths in the CONFIG section below)
#
#  Author : Nadeem Khan
#  Affil  : INRS–Centre Armand-Frappier Santé-Biotechnologie, Laval, QC, Canada
#  Date   : 2025–2026
# ═══════════════════════════════════════════════════════════════════════════════

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(ggalluvial)
  library(ggtern)
  library(viridis)
  library(ggrepel)
})

# ─────────────────────────────────────────────────────────────────────────────
#  CONFIG  —  edit paths here
# ─────────────────────────────────────────────────────────────────────────────

INPUT_FILES <- list(
  Root   = "Root_homologous_expression.tsv",
  Stem   = "Stem_homologous_expression.tsv",
  Leaf   = "Leaf_homologous_expression.tsv",
  Flower = "Flower_homologous_expression.tsv"
)

OUTPUT_DIR  <- "output"          # all plots saved here
FC_THRESH   <- log2(2)           # fold-change threshold for dominance (default: 2x)

# ─────────────────────────────────────────────────────────────────────────────
#  Setup
# ─────────────────────────────────────────────────────────────────────────────

dir.create(OUTPUT_DIR, showWarnings = FALSE)

message("══════════════════════════════════════════════════")
message("  Homeolog Expression Bias Analysis")
message("  Subgenomes: H | St1 | St2")
message("══════════════════════════════════════════════════\n")
message("NOTE: TPM1 = H, TPM2 = St1, TPM3 = St2")
message("      TPM values are pre-averaged across replicates\n")

# ─────────────────────────────────────────────────────────────────────────────
#  1. Data Loading & Validation
# ─────────────────────────────────────────────────────────────────────────────

message("[1/4] Loading input files...")

load_tissue <- function(path, name) {
  if (!file.exists(path)) stop("File not found: ", path)
  df <- read.delim(path, header = TRUE, sep = "\t")
  required <- c("chrH", "chrSt1", "chrSt2", "TPM1", "TPM2", "TPM3")
  missing  <- setdiff(required, colnames(df))
  if (length(missing) > 0) stop("Missing columns in ", name, ": ", paste(missing, collapse = ", "))
  message("  ✔  ", name, " — ", nrow(df), " gene triplets loaded")
  df
}

raw_data <- mapply(load_tissue, INPUT_FILES, names(INPUT_FILES), SIMPLIFY = FALSE)

# Use Root as reference gene order
H_gene_list   <- raw_data$Root$chrH
St1_gene_list <- raw_data$Root$chrSt1
St2_gene_list <- raw_data$Root$chrSt2

# Check gene order consistency across tissues
check_gene_consistency <- function(df, name) {
  if (!all(df$chrH   == H_gene_list))   warning("chrH gene order mismatch in ",   name)
  if (!all(df$chrSt1 == St1_gene_list)) warning("chrSt1 gene order mismatch in ", name)
  if (!all(df$chrSt2 == St2_gene_list)) warning("chrSt2 gene order mismatch in ", name)
}
invisible(mapply(check_gene_consistency,
                 raw_data[c("Stem","Leaf","Flower")],
                 c("Stem","Leaf","Flower")))

# Extract TPM matrices
prepare_expression <- function(df, name) {
  data.frame(
    H_expr   = df$TPM1,
    St1_expr = df$TPM2,
    St2_expr = df$TPM3
  )
}
expression_data <- mapply(prepare_expression, raw_data, names(raw_data), SIMPLIFY = FALSE)

# Retain only rows with complete data across all tissues
valid_rows      <- Reduce("&", lapply(expression_data, complete.cases))
n_removed       <- sum(!valid_rows)
expression_data <- lapply(expression_data, `[`, valid_rows, )
H_gene_list     <- H_gene_list[valid_rows]
St1_gene_list   <- St1_gene_list[valid_rows]
St2_gene_list   <- St2_gene_list[valid_rows]

message("\n  Gene triplets retained : ", sum(valid_rows))
if (n_removed > 0) message("  Rows removed (NA)      : ", n_removed)

# Global expression thresholds
all_tpm         <- unlist(lapply(expression_data, unlist))
non_zero        <- all_tpm[all_tpm > 0]
low_expr_thresh <- quantile(non_zero, 0.25)
pseudocount     <- min(non_zero) / 10

message(sprintf("  Low-expression threshold (25th pctile): %.4f", low_expr_thresh))
message(sprintf("  Pseudocount                           : %.6f\n", pseudocount))

# ─────────────────────────────────────────────────────────────────────────────
#  2. Log2 Ratio Calculation
# ─────────────────────────────────────────────────────────────────────────────

message("[2/4] Calculating log2 expression ratios...")

calculate_homoeolog_ratios <- function(data_list, pseudo) {
  lapply(names(data_list), function(tissue) {
    df <- data_list[[tissue]]
    data.frame(
      GeneID     = H_gene_list,
      H_vs_St1   = log2(pmax(df$H_expr,   pseudo) / pmax(df$St1_expr, pseudo)),
      H_vs_St2   = log2(pmax(df$H_expr,   pseudo) / pmax(df$St2_expr, pseudo)),
      St1_vs_St2 = log2(pmax(df$St1_expr, pseudo) / pmax(df$St2_expr, pseudo)),
      H_expr     = df$H_expr,
      St1_expr   = df$St1_expr,
      St2_expr   = df$St2_expr
    )
  }) |> setNames(names(data_list))
}

ratios <- calculate_homoeolog_ratios(expression_data, pseudocount)
message("  ✔  Ratios computed for all tissues\n")

# ─────────────────────────────────────────────────────────────────────────────
#  3. Dominance Classification
# ─────────────────────────────────────────────────────────────────────────────

message("[3/4] Classifying homeolog dominance...")

classify_homoeolog_dominance <- function(ratios_list, low_thresh, fc_thresh = FC_THRESH) {
  lapply(ratios_list, function(df) {
    dom <- character(nrow(df))
    for (i in seq_len(nrow(df))) {
      h    <- df$H_expr[i]
      s1   <- df$St1_expr[i]
      s2   <- df$St2_expr[i]
      r_h1 <- df$H_vs_St1[i]
      r_h2 <- df$H_vs_St2[i]
      r_12 <- df$St1_vs_St2[i]

      if (any(!is.finite(c(r_h1, r_h2, r_12)))) { dom[i] <- "Invalid";      next }
      if (max(h, s1, s2) < low_thresh)            { dom[i] <- "LowExpr";     next }

      if      (abs(r_h1) > fc_thresh && abs(r_h2) > fc_thresh && r_h1 > 0 && r_h2 > 0) {
        dom[i] <- "H_dominant"
      } else if (abs(r_h1) > fc_thresh && abs(r_12) > fc_thresh && r_h1 < 0 && r_12 > 0) {
        dom[i] <- "St1_dominant"
      } else if (abs(r_h2) > fc_thresh && abs(r_12) > fc_thresh && r_h2 < 0 && r_12 < 0) {
        dom[i] <- "St2_dominant"
      } else {
        dom[i] <- "Balanced"
      }
    }
    df$Dominance <- dom
    df
  })
}

dominance <- classify_homoeolog_dominance(ratios, low_expr_thresh)

# Print dominance summary
for (tissue in names(dominance)) {
  tbl <- table(dominance[[tissue]]$Dominance)
  message(sprintf("  %-8s  %s", tissue,
                  paste(names(tbl), tbl, sep = "=", collapse = " | ")))
}
message()

# ─────────────────────────────────────────────────────────────────────────────
#  4. Visualisation — Log2 Ratio Histograms
# ─────────────────────────────────────────────────────────────────────────────

message("[4/4] Generating plots...")

generate_combined_histogram <- function(ratios_list, subA, subB, pseudo) {

  expr_col_A <- paste0(subA, "_expr")
  expr_col_B <- paste0(subB, "_expr")

  plot_data <- bind_rows(lapply(names(ratios_list), function(tissue) {
    df <- ratios_list[[tissue]]
    ra <- log2(pmax(df[[expr_col_A]], pseudo) / pmax(df[[expr_col_B]], pseudo))
    valid <- is.finite(ra) & ra >= -20 & ra <= 20
    data.frame(Tissue = tissue, Log2_Ratio = ra[valid])
  }))

  # Per-tissue bias statistics
  bias_stats <- plot_data |>
    group_by(Tissue) |>
    summarise(
      Unbiased    = sum(Log2_Ratio >= -1 & Log2_Ratio <= 1),
      SubA_biased = sum(Log2_Ratio >  1),
      SubB_biased = sum(Log2_Ratio < -1),
      Mean_Unbiased = mean(Log2_Ratio[Log2_Ratio >= -1 & Log2_Ratio <= 1], na.rm = TRUE),
      Mean_SubA     = mean(Log2_Ratio[Log2_Ratio >  1],                     na.rm = TRUE),
      Mean_SubB     = mean(Log2_Ratio[Log2_Ratio < -1],                     na.rm = TRUE),
      .groups = "drop"
    ) |>
    pivot_longer(starts_with("Mean_"), names_to = "BiasType", values_to = "Mean") |>
    mutate(
      Bias = case_when(
        BiasType == "Mean_Unbiased" ~ "Unbiased",
        BiasType == "Mean_SubA"     ~ paste(subA, "biased"),
        BiasType == "Mean_SubB"     ~ paste(subB, "biased")
      ),
      Count = case_when(
        Bias == "Unbiased"              ~ Unbiased,
        Bias == paste(subA, "biased")   ~ SubA_biased,
        Bias == paste(subB, "biased")   ~ SubB_biased
      )
    ) |>
    filter(!is.na(Mean)) |>
    group_by(Tissue) |>
    mutate(
      Y_Position = max(
        hist(plot_data$Log2_Ratio[plot_data$Tissue == Tissue[1]],
             breaks = 200, plot = FALSE)$counts
      ) * 0.9
    ) |>
    ungroup()

  color_map <- setNames(
    c("#1f77b4", "#d62728", "#2ca02c"),
    c("Unbiased", paste(subA, "biased"), paste(subB, "biased"))
  )

  ggplot(plot_data, aes(x = Log2_Ratio)) +
    geom_histogram(aes(fill = after_stat(x)), bins = 200,
                   alpha = 0.7, color = "black", linewidth = 0.15) +
    scale_fill_gradientn(
      colors = c("#2ca02c", "#2ca02c", "#c0c0c0", "#1f77b4", "#1f77b4"),
      values = scales::rescale(c(-20, -1, 0, 1, 20)),
      name   = bquote(log[2] ~ "ratio")
    ) +
    geom_vline(xintercept = c(-1, 1),
               linetype = "dashed", color = "red", linewidth = 0.8) +
    geom_vline(data = bias_stats,
               aes(xintercept = Mean, color = Bias),
               linetype = "dotdash", linewidth = 0.9) +
    scale_color_manual(values = color_map, guide = "none") +
    geom_label_repel(
      data          = bias_stats,
      aes(x = Mean, y = Y_Position, color = Bias,
          label = sprintf("%s\nN = %d | mean = %.2f", Bias, Count, Mean)),
      size          = 3.0,
      box.padding   = 0.6,
      point.padding = 0.4,
      segment.color = "grey40",
      segment.size  = 0.4,
      direction     = "y",
      force         = 2.0,
      max.overlaps  = 20
    ) +
    facet_wrap(~ Tissue, ncol = 2) +
    scale_x_continuous(breaks = seq(-20, 20, 5), limits = c(-20, 20)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.35))) +
    theme_minimal(base_size = 13) +
    theme(
      panel.grid      = element_blank(),
      axis.line       = element_line(color = "black", linewidth = 0.5),
      strip.text      = element_text(face = "bold", size = 12),
      legend.position = "right",
      panel.spacing   = unit(1.2, "lines"),
      plot.margin     = margin(10, 15, 10, 10)
    ) +
    labs(
      x = bquote(log[2] ~ "(" * .(subA) / .(subB) * ")"),
      y = "Gene count"
    )
}

# Generate plots for all three pairwise comparisons
comparisons <- list(c("H", "St1"), c("H", "St2"), c("St1", "St2"))

for (comp in comparisons) {
  subA <- comp[1]; subB <- comp[2]
  message(sprintf("  Plotting %s vs %s ...", subA, subB))

  p <- generate_combined_histogram(ratios, subA, subB, pseudocount)

  base <- file.path(OUTPUT_DIR,
                    sprintf("histogram_%s_vs_%s_across_tissues", subA, subB))
  ggsave(paste0(base, ".pdf"), p, width = 10, height = 6, dpi = 400)
  ggsave(paste0(base, ".jpg"), p, width = 10, height = 6, dpi = 400)
  message(sprintf("  ✔  Saved: %s.pdf / .jpg", base))
}

message("\n══════════════════════════════════════════════════")
message("  ✅  Analysis complete! Output in: ./", OUTPUT_DIR)
message("══════════════════════════════════════════════════\n")
