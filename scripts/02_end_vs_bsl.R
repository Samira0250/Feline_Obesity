## =============================================================================
## Adipose Tissue DGE Visualization & GO Enrichment
## Comparison: END vs BSL (direct, not sequential)
## Depots: Abdominal (AB), Subcutaneous (SQ)
## Design: Paired repeated measures (~Animal + Timepoint, 8 cats)
##
## Inputs (place in data/ directory):
##   data/List of DEGs for supplementary.xlsx
##
## Outputs (written to plots2/):
##   volcano_END_vs_BSL.pdf/png
##   heatmaps_END_vs_BSL.pdf
##   DEG_END_vs_BSL_tables.xlsx
##   go_dotplots_END_vs_BSL.pdf
##   GO_results_END_vs_BSL.xlsx
##   individual_volcanos/  individual_heatmaps/  individual_GO/
##
## Note: Run 01_dge_volcanos_heatmaps_GO.R first — output directories are
##       created there. run_go() is redefined here so this script can also
##       be run independently.
## =============================================================================

# Set your working directory to the repo root, or update this path:
# setwd("path/to/adipose-dge")

## -----------------------------------------------------------------------------
## Libraries
## -----------------------------------------------------------------------------

library(openxlsx)
library(tidyverse)
library(ggrepel)
library(patchwork)
library(grid)
library(pheatmap)
library(RColorBrewer)
library(clusterProfiler)
library(org.Hs.eg.db)

## -----------------------------------------------------------------------------
## Output directories
## -----------------------------------------------------------------------------

dir.create("plots2",                      showWarnings = FALSE)
dir.create("plots2/individual_volcanos",  showWarnings = FALSE)
dir.create("plots2/individual_heatmaps",  showWarnings = FALSE)
dir.create("plots2/individual_GO",        showWarnings = FALSE)

## -----------------------------------------------------------------------------
## Sample IDs by depot and timepoint
## -----------------------------------------------------------------------------

ab_bsl_samples <- c("SL364592","SL364598","SL364604","SL364610",
                    "SL368255","SL368261","SL368267","SL368273")
ab_end_samples <- c("SL364596","SL364602","SL364608","SL364614",
                    "SL368259","SL368265","SL368271","SL368277")
sq_bsl_samples <- c("SL364593","SL364599","SL364605","SL364611",
                    "SL368256","SL368262","SL368268","SL368274")
sq_end_samples <- c("SL364597","SL364603","SL364609","SL364615",
                    "SL368260","SL368266","SL368272","SL368278")

## -----------------------------------------------------------------------------
## Color palette for timepoint annotations
## -----------------------------------------------------------------------------

tp_colors <- c("Baseline" = "#66C2A5", "Middle" = "#FC8D62", "End" = "#8DA0CB")

## -----------------------------------------------------------------------------
## 1. LOAD AND ANNOTATE END vs BSL SHEETS
## -----------------------------------------------------------------------------

load_end_vs_bsl <- function(sheet_name) {
  df <- read.xlsx("data/List of DEGs for supplementary.xlsx", sheet = sheet_name)
  
  # Flip LFC: sheets use BSL as reference; positive LFC = upregulated at END
  df$log2FoldChange <- -df$log2FoldChange
  
  # Gene label: use gene symbol; suppress Ensembl IDs from appearing on plots
  df$label <- ifelse(
    !is.na(df$gene_name) & df$gene_name != "" & !grepl("^ENSFCAG", df$gene_name),
    df$gene_name, NA
  )
  
  # Classify DEGs: |log2FC| > 1, padj < 0.05
  df$diffexpressed <- "NO"
  df$diffexpressed[df$log2FoldChange >  1 & !is.na(df$padj) & df$padj < 0.05] <- "UP"
  df$diffexpressed[df$log2FoldChange < -1 & !is.na(df$padj) & df$padj < 0.05] <- "DOWN"
  df$diffexpressed <- factor(df$diffexpressed, levels = c("DOWN", "NO", "UP"))
  
  df
}

evb_data <- list(
  AB_END_BSL = load_end_vs_bsl("Group_1vsGroup_3"),
  SQ_END_BSL = load_end_vs_bsl("Group_4vsGroup_6")
)

# Print DEG counts as a quick sanity check
for (depot in names(evb_data)) {
  cat(depot, "— Up:", sum(evb_data[[depot]]$diffexpressed == "UP",   na.rm = TRUE),
      "| Down:",      sum(evb_data[[depot]]$diffexpressed == "DOWN", na.rm = TRUE), "\n")
}

## -----------------------------------------------------------------------------
## 2. VOLCANO PLOT FUNCTION
## Labels top n significant named genes (no Ensembl IDs)
## x-axis fixed to -8 / +8 with squishing for outliers
## -----------------------------------------------------------------------------

plot_volcano <- function(df, title, n_label = 20) {
  
  # Select top n DEGs to label, ranked by adjusted p-value
  df$delabel    <- NA
  de_idx        <- which(df$diffexpressed != "NO" & !is.na(df$label))
  de_idx_sorted <- de_idx[order(df$padj[de_idx])]
  top_idx       <- de_idx_sorted[seq_len(min(n_label, length(de_idx_sorted)))]
  df$delabel[top_idx] <- df$label[top_idx]
  
  mycolors <- c("DOWN" = "blue", "NO" = "black", "UP" = "red")
  myshapes <- c("DOWN" = 17,     "NO" = 16,       "UP" = 17)
  
  ggplot(df, aes(x = log2FoldChange, y = -log10(padj),
                 col = diffexpressed, shape = diffexpressed, label = delabel)) +
    geom_point(alpha = 0.5, size = 1) +
    geom_text_repel(size = 2.5, max.overlaps = 25,
                    box.padding = 0.3, na.rm = TRUE, show.legend = FALSE) +
    scale_color_manual(values = mycolors, name = "diffexpressed") +
    scale_shape_manual(values = myshapes, name = "diffexpressed") +
    scale_x_continuous(limits = c(-8, 8), oob = scales::squish) +
    geom_vline(xintercept = c(-1, 1), col = "red") +
    geom_hline(yintercept = -log10(0.05), col = "red") +
    labs(title = title, x = "log2FoldChange", y = "-log10(padj)") +
    theme_minimal() +
    theme(
      panel.grid.minor = element_blank(),
      axis.line        = element_blank(),
      legend.position  = "right",
      legend.title     = element_text(size = 8),
      legend.text      = element_text(size = 7),
      plot.title       = element_text(face = "bold", size = 11)
    )
}

## -----------------------------------------------------------------------------
## 3. GENERATE VOLCANO PLOTS — individual PNGs + combined side-by-side grid
## -----------------------------------------------------------------------------

p_ab_evb <- plot_volcano(evb_data$AB_END_BSL, "Abdominal — END vs BSL")
p_sq_evb <- plot_volcano(evb_data$SQ_END_BSL, "Subcutaneous — END vs BSL")

ggsave("plots2/individual_volcanos/volcano_AB_END_vs_BSL.png",
       p_ab_evb, width = 7, height = 5, dpi = 300)
ggsave("plots2/individual_volcanos/volcano_SQ_END_vs_BSL.png",
       p_sq_evb, width = 7, height = 5, dpi = 300)
cat("Saved individual volcanos: AB_END_vs_BSL, SQ_END_vs_BSL\n")

# Combined side-by-side grid with depot and comparison labels
blank     <- wrap_elements(full = textGrob("", gp = gpar(fontsize = 1)))
abd_title <- wrap_elements(full = textGrob("Abdominal Adipose Tissue",
                                           gp = gpar(fontface = "bold", fontsize = 12)))
sub_title <- wrap_elements(full = textGrob("Subcutaneous Adipose Tissue",
                                           gp = gpar(fontface = "bold", fontsize = 12)))
row_evb   <- wrap_elements(full = textGrob("END vs BSL", rot = 90,
                                           gp = gpar(fontface = "bold", fontsize = 11)))

volcano_evb_grid <-
  (blank   | abd_title | sub_title) /
  (row_evb | p_ab_evb  | p_sq_evb) +
  plot_layout(heights = c(0.1, 1), widths = c(0.07, 1, 1))

ggsave("plots2/volcano_END_vs_BSL.pdf", volcano_evb_grid, width = 14, height = 6, dpi = 300)
ggsave("plots2/volcano_END_vs_BSL.png", volcano_evb_grid, width = 14, height = 6, dpi = 300)
cat("Saved: volcano_END_vs_BSL\n")

## -----------------------------------------------------------------------------
## 4. HEATMAPS — top 50 DEGs, Z-scored pseudo counts
## -----------------------------------------------------------------------------

heatmap_evb_configs <- list(
  list(
    label      = "AB END vs BSL",
    de_df      = evb_data$AB_END_BSL,
    expr_sheet = "Group_1vsGroup_3",
    grp1_pat   = ab_bsl_samples,
    grp2_pat   = ab_end_samples,
    grp1_name  = "Baseline", grp2_name = "End"
  ),
  list(
    label      = "SQ END vs BSL",
    de_df      = evb_data$SQ_END_BSL,
    expr_sheet = "Group_4vsGroup_6",
    grp1_pat   = sq_bsl_samples,
    grp2_pat   = sq_end_samples,
    grp1_name  = "Baseline", grp2_name = "End"
  )
)

# Helper: build and draw a single heatmap given a config and a top-gene data frame
draw_heatmap <- function(cfg, top50) {
  
  expr_df     <- read.xlsx("data/List of DEGs for supplementary.xlsx", sheet = cfg$expr_sheet)
  pseudo_cols <- grep("\\.pseudo$", colnames(expr_df), value = TRUE)
  expr_mat    <- as.matrix(expr_df[, pseudo_cols])
  rownames(expr_mat) <- expr_df$Geneid
  colnames(expr_mat) <- gsub("\\.pseudo$", "", colnames(expr_mat))
  
  grp1_cols <- intersect(cfg$grp1_pat, colnames(expr_mat))
  grp2_cols <- intersect(cfg$grp2_pat, colnames(expr_mat))
  mat       <- expr_mat[top50$Geneid, c(grp1_cols, grp2_cols), drop = FALSE]
  rownames(mat) <- top50$label
  mat_z     <- t(scale(t(mat)))   # Z-score across samples per gene
  
  ann_col <- data.frame(
    Timepoint = c(rep(cfg$grp1_name, length(grp1_cols)),
                  rep(cfg$grp2_name, length(grp2_cols))),
    row.names = c(grp1_cols, grp2_cols)
  )
  ann_colours <- list(Timepoint = tp_colors[unique(ann_col$Timepoint)])
  
  pheatmap(mat_z,
           annotation_col    = ann_col,
           annotation_colors = ann_colours,
           color             = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
           cluster_rows      = TRUE,  cluster_cols  = FALSE,
           show_rownames     = TRUE,  show_colnames = TRUE,
           fontsize_row = 7,  fontsize_col = 8,
           main         = cfg$label, border_color = NA,
           cutree_rows  = 2)
}

# Save all heatmaps to a single PDF and individual PNGs
pdf("plots2/heatmaps_END_vs_BSL.pdf", width = 10, height = 12)
for (cfg in heatmap_evb_configs) {
  
  top50 <- cfg$de_df %>%
    filter(diffexpressed != "NO") %>%
    arrange(padj) %>%
    slice_head(n = 50) %>%
    mutate(label = ifelse(!is.na(gene_name) & gene_name != "", gene_name, Geneid))
  
  if (nrow(top50) < 2) { message("Skipping heatmap (too few DEGs): ", cfg$label); next }
  
  draw_heatmap(cfg, top50)
  
  # Also save individual PNG
  fname <- paste0("plots2/individual_heatmaps/heatmap_", gsub(" ", "_", cfg$label), ".png")
  png(fname, width = 10, height = 12, units = "in", res = 300)
  draw_heatmap(cfg, top50)
  dev.off()
  cat("Saved heatmap:", cfg$label, "\n")
}
dev.off()
cat("Saved: heatmaps_END_vs_BSL.pdf\n")

## -----------------------------------------------------------------------------
## 5. UP/DOWN DEG TABLES → Excel
## -----------------------------------------------------------------------------

out_cols <- c("Geneid", "gene_name", "Gene.description",
              "baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj")

wb_evb <- createWorkbook()

# Summary tab: gene counts per depot
summary_evb <- data.frame(
  Comparison = names(evb_data),
  Up   = sapply(evb_data, function(df) sum(df$diffexpressed == "UP",   na.rm = TRUE)),
  Down = sapply(evb_data, function(df) sum(df$diffexpressed == "DOWN", na.rm = TRUE))
)
addWorksheet(wb_evb, "Summary")
writeData(wb_evb, "Summary", summary_evb)

# One sheet per UP/DOWN subset
for (depot in names(evb_data)) {
  df <- evb_data[[depot]]
  
  up_df <- df[!is.na(df$diffexpressed) & df$diffexpressed == "UP",   out_cols]
  dn_df <- df[!is.na(df$diffexpressed) & df$diffexpressed == "DOWN", out_cols]
  
  addWorksheet(wb_evb, paste0(depot, "_UP"))
  writeData(wb_evb,   paste0(depot, "_UP"),   up_df)
  addWorksheet(wb_evb, paste0(depot, "_DOWN"))
  writeData(wb_evb,   paste0(depot, "_DOWN"), dn_df)
}

saveWorkbook(wb_evb, "plots2/DEG_END_vs_BSL_tables.xlsx", overwrite = TRUE)
cat("Saved: DEG_END_vs_BSL_tables.xlsx\n")

## -----------------------------------------------------------------------------
## 6. GO BIOLOGICAL PROCESS ENRICHMENT
## Uses human gene symbols (org.Hs.eg.db) as proxy for feline gene annotation
## -----------------------------------------------------------------------------

run_go <- function(gene_symbols, universe_symbols) {
  sig <- toupper(gene_symbols[!is.na(gene_symbols) & gene_symbols != ""])
  uni <- toupper(universe_symbols[!is.na(universe_symbols) & universe_symbols != ""])
  
  if (length(sig) < 5) return(NULL)
  
  eg_sig <- bitr(sig, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  eg_uni <- bitr(uni, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  
  if (nrow(eg_sig) == 0) return(NULL)
  
  ego <- enrichGO(
    gene          = eg_sig$ENTREZID,
    universe      = eg_uni$ENTREZID,
    OrgDb         = org.Hs.eg.db,
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.20,
    readable      = TRUE
  )
  
  if (is.null(ego) || nrow(ego) == 0) return(NULL)
  ego
}

go_evb <- list()
for (depot in names(evb_data)) {
  cat("Running GO enrichment:", depot, "\n")
  
  df       <- evb_data[[depot]]
  universe <- df$gene_name
  up_genes <- df$gene_name[!is.na(df$diffexpressed) & df$diffexpressed == "UP"]
  dn_genes <- df$gene_name[!is.na(df$diffexpressed) & df$diffexpressed == "DOWN"]
  
  go_evb[[paste0(depot, "_UP")]]   <- run_go(up_genes, universe)
  go_evb[[paste0(depot, "_DOWN")]] <- run_go(dn_genes, universe)
  go_evb[[paste0(depot, "_both")]] <- run_go(c(up_genes, dn_genes), universe)
}

# Save dotplots to a single PDF
pdf("plots2/go_dotplots_END_vs_BSL.pdf", width = 10, height = 8)
for (nm in names(go_evb)) {
  if (!is.null(go_evb[[nm]]) && nrow(go_evb[[nm]]) > 0) {
    print(
      dotplot(go_evb[[nm]], showCategory = 20, font.size = 9,
              title = gsub("_", " ", nm)) +
        theme(plot.title = element_text(face = "bold"))
    )
  }
}
dev.off()

# Save individual PNG per comparison
for (nm in names(go_evb)) {
  if (!is.null(go_evb[[nm]]) && nrow(go_evb[[nm]]) > 0) {
    p_go <- dotplot(go_evb[[nm]], showCategory = 20, font.size = 9,
                    title = gsub("_", " ", nm)) +
      theme(plot.title = element_text(face = "bold"))
    ggsave(paste0("plots2/individual_GO/GO_", nm, ".png"),
           p_go, width = 10, height = 8, dpi = 300)
    cat("Saved GO PNG:", nm, "\n")
  }
}

# Save GO results to Excel
go_out_cols <- c("ID", "Description", "GeneRatio", "BgRatio",
                 "pvalue", "p.adjust", "qvalue", "Count", "geneID")
wb_go_evb <- createWorkbook()
for (nm in names(go_evb)) {
  if (!is.null(go_evb[[nm]]) && nrow(go_evb[[nm]]) > 0) {
    sheet_name <- substr(paste0("GO_", nm), 1, 31)
    addWorksheet(wb_go_evb, sheet_name)
    writeData(wb_go_evb, sheet_name,
              as.data.frame(go_evb[[nm]]) %>%
                dplyr::select(all_of(go_out_cols)) %>%
                arrange(p.adjust))
  }
}
saveWorkbook(wb_go_evb, "plots2/GO_results_END_vs_BSL.xlsx", overwrite = TRUE)
cat("Saved: go_dotplots_END_vs_BSL.pdf + GO_results_END_vs_BSL.xlsx\n")

cat("\nEND vs BSL analysis complete. All outputs written to plots2/\n")
