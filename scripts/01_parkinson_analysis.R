# =============================================================================
# Parkinson's Disease — RNA-Seq Exploratory Analysis & Biomarker Discovery
# =============================================================================
# Author      : Ibrahim (Bembo)
# Course      : Multiomics Data Analysis — Bioinformatics Diploma
# Date        : April 2025
# Description : Identification of most variable genes as potential biomarkers
#               for Parkinson's disease using microarray expression data.
#
# Pipeline:
#   1. Load & validate data
#   2. Quality check — Histogram & Density Plot
#   3. Dimensionality reduction — PCA (2D & 3D)
#   4. Most variable genes — rowVars()
#   5. Heatmaps — Absolute expression & Z-score
# =============================================================================


# ── 0. Install packages (run once) ───────────────────────────────────────────
# install.packages(c("ggplot2", "plotly", "matrixStats", "htmlwidgets"))
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("ComplexHeatmap")


# ── 1. Load libraries ─────────────────────────────────────────────────────────
library(ggplot2)
library(plotly)
library(matrixStats)
library(ComplexHeatmap)
library(circlize)
library(htmlwidgets)


# ── 2. Set paths ──────────────────────────────────────────────────────────────
data_dir   <- file.path(dirname(getwd()), "data")
output_dir <- file.path(dirname(getwd()), "results", "plots")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)


# ── 3. Load data ──────────────────────────────────────────────────────────────
message("[ 1/5 ] Loading data...")

exp_raw <- read.table(
  file.path(data_dir, "Parkinson_exp.txt"),
  header      = TRUE,
  sep         = "\t",
  row.names   = 1,
  check.names = FALSE
)

pheno <- read.table(
  file.path(data_dir, "Parkinson_phenotable.txt"),
  header    = TRUE,
  sep       = "\t",
  row.names = 1
)

# Validate
stopifnot(
  "Sample names mismatch between expression and pheno table!" =
    all(colnames(exp_raw) == rownames(pheno))
)

message(sprintf("  Expression matrix : %d genes x %d samples", nrow(exp_raw), ncol(exp_raw)))
message(sprintf("  Phenotype table   : %d samples x %d variables", nrow(pheno), ncol(pheno)))
message(sprintf("  Groups            : %s", paste(table(pheno$sample.type), collapse = " Ctrl / "), "PD"))


# ── 4. QC — Histogram & Density Plot ─────────────────────────────────────────
message("[ 2/5 ] Quality check — Histogram & Density Plot...")

exp_values <- as.vector(as.matrix(exp_raw))

# 4a. Histogram
pdf(file.path(output_dir, "01_Histogram_Expression.pdf"), width = 8, height = 6)
hist(exp_values,
     breaks = 100,
     col    = "steelblue",
     border = "white",
     main   = "Distribution of Expression Values — Parkinson Dataset",
     xlab   = "Expression Value (log2 scale)",
     ylab   = "Frequency")
dev.off()

# 4b. Density plot (per sample)
exp_long        <- stack(as.data.frame(exp_raw))
colnames(exp_long) <- c("expression", "sample")
exp_long$group  <- pheno[exp_long$sample, "sample.type"]

p_density <- ggplot(exp_long, aes(x = expression, color = sample, linetype = group)) +
  geom_density(linewidth = 0.6) +
  theme_minimal(base_size = 12) +
  labs(
    title    = "Density Plot of Expression Values per Sample",
    subtitle = "Parkinson's Disease (PD) vs Controls (Ctrl)",
    x        = "Expression Value (log2 scale)",
    y        = "Density",
    color    = "Sample",
    linetype = "Group"
  ) +
  theme(legend.position = "right")

pdf(file.path(output_dir, "02_Density_Plot_per_Sample.pdf"), width = 11, height = 6)
print(p_density)
dev.off()


# ── 5. PCA — Dimensionality Reduction ─────────────────────────────────────────
message("[ 3/5 ] PCA — 2D & 3D...")

pca_result <- prcomp(t(exp_raw), scale. = TRUE)
variance   <- summary(pca_result)$importance[2, 1:3] * 100

# 5a. PCA 2D
pca_df         <- as.data.frame(pca_result$x[, 1:2])
pca_df$sample  <- rownames(pca_df)
pca_df$group   <- pheno[rownames(pca_df), "sample.type"]
pca_df$gender  <- pheno[rownames(pca_df), "gender"]

p_pca2d <- ggplot(pca_df, aes(x = PC1, y = PC2,
                               color = group, shape = gender, label = sample)) +
  geom_point(size = 4, alpha = 0.85) +
  geom_text(vjust = -0.9, size = 2.8, show.legend = FALSE) +
  stat_ellipse(aes(group = group), linetype = "dashed", linewidth = 0.4) +
  theme_minimal(base_size = 12) +
  labs(
    title    = "PCA — Parkinson's Disease vs Control",
    subtitle = paste0("PC1 explains ", round(variance[1], 1),
                      "% | PC2 explains ", round(variance[2], 1), "% of variance"),
    x        = paste0("PC1 (", round(variance[1], 1), "%)"),
    y        = paste0("PC2 (", round(variance[2], 1), "%)"),
    color    = "Group",
    shape    = "Gender"
  ) +
  scale_color_manual(values = c("Ctrl" = "steelblue", "PD" = "tomato")) +
  theme(plot.subtitle = element_text(size = 10, color = "gray50"))

pdf(file.path(output_dir, "03_PCA_2D.pdf"), width = 9, height = 7)
print(p_pca2d)
dev.off()

# 5b. PCA 3D
pca_3d         <- as.data.frame(pca_result$x[, 1:3])
pca_3d$group   <- pheno[rownames(pca_3d), "sample.type"]
pca_3d$gender  <- pheno[rownames(pca_3d), "gender"]
pca_3d$sample  <- rownames(pca_3d)

p_pca3d <- plot_ly(
  data   = pca_3d,
  x = ~PC1, y = ~PC2, z = ~PC3,
  color  = ~group,
  symbol = ~gender,
  colors = c("Ctrl" = "steelblue", "PD" = "tomato"),
  text   = ~sample,
  type   = "scatter3d",
  mode   = "markers",
  marker = list(size = 6)
) %>%
  layout(
    title = "PCA 3D — Parkinson's Disease vs Control",
    scene = list(
      xaxis = list(title = paste0("PC1 (", round(variance[1], 1), "%)")),
      yaxis = list(title = paste0("PC2 (", round(variance[2], 1), "%)")),
      zaxis = list(title = paste0("PC3 (", round(variance[3], 1), "%)"))
    )
  )

saveWidget(p_pca3d,
           file          = file.path(output_dir, "04_PCA_3D.html"),
           selfcontained = TRUE)


# ── 6. Most Variable Genes ────────────────────────────────────────────────────
message("[ 4/5 ] Calculating row variances & selecting top 100 genes...")

gene_vars   <- rowVars(as.matrix(exp_raw))
names(gene_vars) <- rownames(exp_raw)

top100_idx  <- order(gene_vars, decreasing = TRUE)[1:100]
top100_expr <- exp_raw[top100_idx, ]

# Save top 100 gene list
top100_table <- data.frame(
  probe_gene = names(sort(gene_vars, decreasing = TRUE)[1:100]),
  variance   = sort(gene_vars, decreasing = TRUE)[1:100]
)
write.csv(top100_table,
          file      = file.path(output_dir, "../top100_variable_genes.csv"),
          row.names = FALSE)

message("  Top 5 most variable genes:")
print(head(top100_table, 5))


# ── 7. Heatmaps ───────────────────────────────────────────────────────────────
message("[ 5/5 ] Generating heatmaps...")

col_anno <- HeatmapAnnotation(
  Group  = pheno$sample.type,
  Gender = pheno$gender,
  col    = list(
    Group  = c("Ctrl" = "steelblue", "PD" = "tomato"),
    Gender = c("Female" = "#F4A0B5", "Male" = "#A0C4F4")
  ),
  annotation_name_gp = gpar(fontsize = 9)
)

# 7a. Absolute expression
pdf(file.path(output_dir, "05_Heatmap_Absolute_Expression.pdf"),
    width = 12, height = 10)
draw(
  Heatmap(
    as.matrix(top100_expr),
    name                 = "Expression",
    top_annotation       = col_anno,
    show_row_names       = FALSE,
    show_column_names    = TRUE,
    cluster_rows         = TRUE,
    cluster_columns      = TRUE,
    col                  = colorRamp2(c(6, 9, 14), c("blue", "white", "red")),
    column_names_gp      = gpar(fontsize = 7),
    column_title         = "Top 100 Most Variable Genes — Absolute Expression",
    column_title_gp      = gpar(fontsize = 12, fontface = "bold"),
    heatmap_legend_param = list(title = "Expression\n(log2)")
  )
)
dev.off()

# 7b. Z-score (relative expression)
top100_zscore <- t(scale(t(as.matrix(top100_expr))))

pdf(file.path(output_dir, "06_Heatmap_Zscore.pdf"),
    width = 12, height = 10)
draw(
  Heatmap(
    top100_zscore,
    name                 = "Z-score",
    top_annotation       = col_anno,
    show_row_names       = FALSE,
    show_column_names    = TRUE,
    cluster_rows         = TRUE,
    cluster_columns      = TRUE,
    col                  = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
    column_names_gp      = gpar(fontsize = 7),
    column_title         = "Top 100 Most Variable Genes — Z-score (Relative Expression)",
    column_title_gp      = gpar(fontsize = 12, fontface = "bold"),
    heatmap_legend_param = list(title = "Z-score")
  )
)
dev.off()


# ── 8. Done ───────────────────────────────────────────────────────────────────
message("\n✅ Analysis complete! Output files saved to:")
message("   ", normalizePath(output_dir))
message("\nFiles generated:")
for (f in list.files(output_dir)) message("   • ", f)
message("\nTop 100 variable genes table saved to results/top100_variable_genes.csv\n")
