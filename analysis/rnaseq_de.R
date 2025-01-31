# Load required libraries
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

# Setup paths and directories
salmon_dir        <- "results/star_salmon"
sample_sheet_file <- "data/samplesheet_Tomo.csv"
tx2gene_file      <- "results/star_salmon/tx2gene.tsv"

dir.create("results_dge", showWarnings = FALSE)
dir.create("plots_dge", showWarnings = FALSE)

# Read and process sample sheet
sample_sheet <- read_csv(sample_sheet_file) %>%
  mutate(
    replicate      = factor(replicate),
    gene_target    = factor(gene_target),
    DNAPKi_treated = factor(DNAPKi_treated, levels = c("FALSE","TRUE")),
    condition      = factor(paste0(gene_target, "_", DNAPKi_treated))
  ) %>%
  mutate(condition = relevel(condition, ref = "NTC_FALSE"))

# Import Salmon counts
all_dirs    <- list.dirs(salmon_dir, full.names = TRUE, recursive = FALSE)
sample_dirs <- all_dirs[basename(all_dirs) %in% sample_sheet$sample_name]
sample_sheet<- sample_sheet %>% filter(sample_name %in% basename(sample_dirs))
files       <- file.path(sample_dirs, "quant.sf")
names(files)<- basename(sample_dirs)

tx2gene <- read_tsv(tx2gene_file, col_names = c("transcript_id","gene_id"))
txi     <- tximport(files, type = "salmon", tx2gene = tx2gene, dropInfReps = TRUE)

stopifnot(all(colnames(txi$counts) == sample_sheet$sample_name))

# Create DGEList and filter low expression genes
dge <- DGEList(counts = txi$counts, samples = sample_sheet)
keep_expr <- filterByExpr(dge, group = dge$samples$condition)
dge <- dge[keep_expr, keep.lib.sizes=FALSE]
dge <- calcNormFactors(dge)

# Design matrix: ~ replicate + condition
# - Accounts for replicate effects
# - Tests each condition against NTC_FALSE (reference level)
design <- model.matrix(~ replicate + condition, data = dge$samples)
v   <- voom(dge, design, plot = TRUE)
fit <- lmFit(v, design)
fit <- eBayes(fit, robust = TRUE)

# Helper function to extract DE results with gene symbols
extract_results <- function(fit_obj, coef_name) {
  topTable(fit_obj, coef = coef_name, number = Inf, sort.by = "P") %>%
    rownames_to_column("gene_id") %>%
    mutate(
      gene_symbol = mapIds(
        org.Hs.eg.db,
        keys   = sub("\\..*$", "", gene_id),
        keytype= "ENSEMBL",
        column = "SYMBOL"
      )
    ) %>%
    dplyr::select(gene_id, gene_symbol, logFC, AveExpr, t, P.Value, adj.P.Val, B)
}

# Single comparisons vs NTC_FALSE
all_conditions <- setdiff(levels(dge$samples$condition), "NTC_FALSE")
single_comparison_results <- list()
for (cond in all_conditions) {
  coef_name <- paste0("condition", cond)
  res       <- extract_results(fit, coef_name)
  write_csv(res, paste0("results_dge/DE_", cond, "_vs_NTC_FALSE.csv"))
  single_comparison_results[[cond]] <- res
}

# Contrasts for DNAPKi effects within each gene target
# Tests: (Target_TRUE - Target_FALSE) for each target
contr_matrix <- makeContrasts(
  GAL4_DNAPKi_effect  = conditionGAL4_TRUE - conditionGAL4_FALSE,
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

# Volcano plot function
create_volcano_plot <- function(df, title_str) {
  df <- df %>%
    mutate(signif_cat = case_when(
      adj.P.Val < 0.05 & logFC >  1  ~ "Up",
      adj.P.Val < 0.05 & logFC < -1  ~ "Down",
      TRUE                           ~ "NotSig"
    ))
  
  lowest_pval <- df %>% arrange(P.Value) %>% head(10)
  largest_fc_pos <- df %>% arrange(desc(logFC)) %>% head(5)
  largest_fc_neg <- df %>% arrange(logFC) %>% head(5)
  
  highlight_genes <- bind_rows(lowest_pval, largest_fc_pos, largest_fc_neg) %>%
    distinct(gene_id, .keep_all = TRUE)
  
  p <- ggplot(df, aes(x = logFC, y = -log10(P.Value), color = signif_cat)) +
    geom_point(alpha = 0.6) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40") +
    scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NotSig" = "grey60")) +
    labs(title = title_str, x = "log2 Fold Change", y = "-log10 P-value") +
    theme_bw() +
    theme(
      legend.position = "bottom",
      plot.title     = element_text(hjust = 0.5)
    )
  
  p + geom_text_repel(
    data = highlight_genes,
    aes(label = gene_symbol),
    size = 3,
    max.overlaps = Inf
  )
}

# Generate volcano plots
single_volcano_list <- list()
for (cond in names(single_comparison_results)) {
  plt_title <- paste0(cond, " vs NTC_FALSE")
  single_volcano_list[[cond]] <- create_volcano_plot(single_comparison_results[[cond]], plt_title)
}
combined_single_volcano <- wrap_plots(single_volcano_list, ncol = 3)
ggsave("plots_dge/single_comparisons_volcano.png", 
       combined_single_volcano, width = 15, height = 15, dpi = 300)

double_volcano_list <- list()
for (nm in names(double_comparison_results)) {
  double_volcano_list[[nm]] <- create_volcano_plot(
    double_comparison_results[[nm]], 
    paste0(nm, " (DNAPKi effect)")
  )
}
combined_double_volcano <- wrap_plots(double_volcano_list, ncol = 2)
ggsave("plots_dge/double_comparisons_volcano.png", 
       combined_double_volcano, width = 10, height = 5, dpi = 300)

message("Analysis complete. Results in 'results_dge/', plots in 'plots_dge/'.")