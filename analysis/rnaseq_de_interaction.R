###############################################################################
# Load required libraries
###############################################################################
library(tximport)
library(readr)
library(dplyr)
library(tidyr)
library(tibble)
library(edgeR)
library(limma)
library(org.Hs.eg.db)
library(ggplot2)
library(ggrepel)
library(ComplexHeatmap)
library(patchwork)

###############################################################################
# Setup paths and directories
###############################################################################
salmon_dir        <- "results/star_salmon"
sample_sheet_file <- "data/samplesheet_Tomo.csv"
tx2gene_file      <- "results/star_salmon/tx2gene.tsv"

dir.create("results_dge", showWarnings = FALSE)
dir.create("plots_dge",   showWarnings = FALSE)

###############################################################################
# Read and process sample sheet
###############################################################################
sample_sheet <- read_csv(sample_sheet_file) %>%
  mutate(
    replicate      = factor(replicate),
    # Ensure gene_target and DNAPKi_treated have appropriate reference levels
    # so that NTC_FALSE becomes the reference in the single-condition analysis,
    # and so that GAL4_FALSE is the reference in the interaction analysis if desired.
    gene_target    = factor(gene_target),
    DNAPKi_treated = factor(DNAPKi_treated, levels = c("FALSE","TRUE")),
    condition      = factor(paste0(gene_target, "_", DNAPKi_treated))
  ) %>%
  # Relevel condition so that NTC_FALSE is the reference for single-comparison design
  mutate(condition = relevel(condition, ref = "NTC_FALSE"))

###############################################################################
# Import Salmon counts
###############################################################################
all_dirs    <- list.dirs(salmon_dir, full.names = TRUE, recursive = FALSE)
sample_dirs <- all_dirs[basename(all_dirs) %in% sample_sheet$sample_name]
sample_sheet<- sample_sheet %>% filter(sample_name %in% basename(sample_dirs))
files       <- file.path(sample_dirs, "quant.sf")
names(files)<- basename(sample_dirs)

tx2gene <- read_tsv(tx2gene_file, col_names = c("transcript_id","gene_id"))
txi     <- tximport(files, type = "salmon", tx2gene = tx2gene, dropInfReps = TRUE)

###############################################################################
# Create DGEList and filter low expression genes
###############################################################################
dge <- DGEList(counts = txi$counts, samples = sample_sheet)

keep_expr <- filterByExpr(dge, group = dge$samples$condition)
dge <- dge[keep_expr, keep.lib.sizes=FALSE]
dge <- calcNormFactors(dge)

###############################################################################
# Helper function to extract DE results with gene symbols
###############################################################################
extract_results <- function(fit_obj, coef_name) {
  topTable(fit_obj, coef = coef_name, number = Inf, sort.by = "P") %>%
    rownames_to_column("gene_id") %>%
    mutate(
      gene_symbol = mapIds(
        org.Hs.eg.db,
        keys    = sub("\\..*$", "", gene_id),
        keytype = "ENSEMBL",
        column  = "SYMBOL"
      )
    ) %>%
    dplyr::select(gene_id, gene_symbol, logFC, AveExpr, t, P.Value, adj.P.Val, B)
}

###############################################################################
# Volcano plot function
###############################################################################
create_volcano_plot <- function(df, title_str) {
  df <- df %>%
    mutate(signif_cat = case_when(
      adj.P.Val < 0.05 & logFC >  1  ~ "Up",
      adj.P.Val < 0.05 & logFC < -1  ~ "Down",
      TRUE                           ~ "NotSig"
    ))
  
  # Highlight a few interesting genes
  lowest_pval    <- df %>% arrange(P.Value) %>% head(10)
  largest_fc_pos <- df %>% arrange(desc(logFC)) %>% head(5)
  largest_fc_neg <- df %>% arrange(logFC) %>% head(5)
  highlight_genes <- bind_rows(lowest_pval, largest_fc_pos, largest_fc_neg) %>%
    distinct(gene_id, .keep_all = TRUE)
  
  p <- ggplot(df, aes(x = logFC, y = -log10(P.Value), color = signif_cat)) +
    geom_point(alpha = 0.6) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
    geom_vline(xintercept = c(-1, 1),      linetype = "dashed", color = "grey40") +
    scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NotSig" = "grey60")) +
    labs(title = title_str, x = "log2 Fold Change", y = "-log10 P-value") +
    theme_bw() +
    theme(
      legend.position = "bottom",
      plot.title      = element_text(hjust = 0.5)
    )
  
  p + geom_text_repel(
    data = highlight_genes,
    aes(label = gene_symbol),
    size = 3,
    max.overlaps = Inf
  )
}

###############################################################################
# 1) DIFFERENTIAL EXPRESSION: Single-condition design (as before)
#    ~ replicate + condition
###############################################################################
design <- model.matrix(~ replicate + condition, data = dge$samples)

v   <- voom(dge, design, plot = FALSE)
fit <- lmFit(v, design)
fit <- eBayes(fit, robust = TRUE)

# Single comparisons vs NTC_FALSE
all_conditions <- setdiff(levels(dge$samples$condition), "NTC_FALSE")
single_comparison_results <- list()

for (cond in all_conditions) {
  coef_name <- paste0("condition", cond)
  res       <- extract_results(fit, coef_name)
  write_csv(res, paste0("results_dge/DE_", cond, "_vs_NTC_FALSE.csv"))
  single_comparison_results[[cond]] <- res
}

# Volcano plots for single comparisons
single_volcano_list <- list()
for (cond in names(single_comparison_results)) {
  plt_title <- paste0(cond, " vs NTC_FALSE")
  single_volcano_list[[cond]] <- create_volcano_plot(single_comparison_results[[cond]], plt_title)
}
combined_single_volcano <- wrap_plots(single_volcano_list, ncol = 3)

ggsave(
  filename = "plots_dge/single_comparisons_volcano.png",
  plot     = combined_single_volcano,
  width    = 15,
  height   = 15,
  dpi      = 300
)

# Contrasts for DNAPKi effects within each gene target
contr_matrix <- makeContrasts(
  GAL4_DNAPKi_effect  = conditionGAL4_TRUE  - conditionGAL4_FALSE,
  PRDX1_DNAPKi_effect = conditionPRDX1_TRUE - conditionPRDX1_FALSE,
  levels = design
)
contr_fit <- contrasts.fit(fit, contr_matrix)
contr_fit <- eBayes(contr_fit, robust = TRUE)

double_comparison_results <- list(
  GAL4_DNAPKi_effect  = extract_results(contr_fit, "GAL4_DNAPKi_effect"),
  PRDX1_DNAPKi_effect = extract_results(contr_fit, "PRDX1_DNAPKi_effect")
)

for (nm in names(double_comparison_results)) {
  write_csv(double_comparison_results[[nm]], paste0("results_dge/DE_", nm, ".csv"))
}

# Volcano plots for DNAPKi contrasts
double_volcano_list <- list()
for (nm in names(double_comparison_results)) {
  double_volcano_list[[nm]] <- create_volcano_plot(
    double_comparison_results[[nm]],
    paste0(nm, " (DNAPKi effect)")
  )
}
combined_double_volcano <- wrap_plots(double_volcano_list, ncol = 2)

ggsave(
  filename = "plots_dge/double_comparisons_volcano.png",
  plot     = combined_double_volcano,
  width    = 10,
  height   = 5,
  dpi      = 300
)

###############################################################################
# 2) DIFFERENTIAL EXPRESSION WITH INTERACTION TERMS
#    ~ replicate + gene_target * DNAPKi_treated
#
# This design estimates:
# - Intercept = expression in [gene_target = reference, DNAPKi_treated = reference]
# - Main effect of gene_target (PRDX1 vs GAL4)
# - Main effect of DNAPKi_treated (TRUE vs FALSE)
# - Interaction: difference in difference (PRDX1:TRUE vs. everything else)
###############################################################################
# Make sure factor levels are set so that:
# gene_target = "GAL4" is reference
# DNAPKi_treated = "FALSE" is reference
# If not, you can do:
# dge$samples$gene_target    <- relevel(dge$samples$gene_target, "GAL4")
# dge$samples$DNAPKi_treated <- relevel(dge$samples$DNAPKi_treated, "FALSE")

design_interaction <- model.matrix(~ replicate + gene_target * DNAPKi_treated, data = dge$samples)

v_int   <- voom(dge, design_interaction, plot = FALSE)
fit_int <- lmFit(v_int, design_interaction)
fit_int <- eBayes(fit_int, robust = TRUE)

# Extract DE results for each coefficient:
# 1) gene_targetPRDX1              : main effect of PRDX1 vs GAL4 (in absence of DNAPKi)
# 2) DNAPKi_treatedTRUE            : main effect of DNAPKi (TRUE vs FALSE) in GAL4 background
# 3) gene_targetPRDX1:DNAPKi_treatedTRUE : interaction effect
interaction_results <- list(
  gene_targetPRDX1                 = extract_results(fit_int, "gene_targetPRDX1"),
  DNAPKi_treatedTRUE               = extract_results(fit_int, "DNAPKi_treatedTRUE"),
  gene_targetPRDX1_DNAPKi_treatedTRUE = extract_results(fit_int, "gene_targetPRDX1:DNAPKi_treatedTRUE")
)

# Write out CSV results for each coefficient
write_csv(interaction_results$gene_targetPRDX1,
          "results_dge/DE_interaction_gene_targetPRDX1.csv")
write_csv(interaction_results$DNAPKi_treatedTRUE,
          "results_dge/DE_interaction_DNAPKiTRUE.csv")
write_csv(interaction_results$gene_targetPRDX1_DNAPKi_treatedTRUE,
          "results_dge/DE_interaction_PRDX1xDNAPKi.csv")

# Create and save volcano plots for each interaction term
interaction_volcano_list <- list()
interaction_volcano_list[["gene_targetPRDX1"]] <- create_volcano_plot(
  interaction_results[["gene_targetPRDX1"]], 
  "Main effect: PRDX1 vs GAL4 (no DNAPKi)"
)
interaction_volcano_list[["DNAPKi_treatedTRUE"]] <- create_volcano_plot(
  interaction_results[["DNAPKi_treatedTRUE"]], 
  "Main effect: DNAPKi (TRUE vs FALSE) in GAL4"
)
interaction_volcano_list[["gene_targetPRDX1:DNAPKi_treatedTRUE"]] <- create_volcano_plot(
  interaction_results[["gene_targetPRDX1_DNAPKi_treatedTRUE"]], 
  "Interaction: (PRDX1:DNAPKi) - [GAL4:no_DNAPKi]"
)

combined_interaction_volcano <- wrap_plots(interaction_volcano_list, ncol = 3)

ggsave(
  filename = "plots_dge/interaction_terms_volcano.png",
  plot     = combined_interaction_volcano,
  width    = 18,
  height   = 6,
  dpi      = 300
)

###############################################################################
# End of script
###############################################################################
message("Analysis complete. Results in 'results_dge/', plots in 'plots_dge/'.")
