## =============================================================================
## 01_dge_volcanos_heatmaps_GO.R
##
## Adipose Tissue DGE Visualization & GO Enrichment
## Comparisons: MID vs BSL, END vs MID
## Depots: Abdominal (AB), Subcutaneous (SQ)
## Design: Paired repeated measures (~Animal + Timepoint, 8 cats)
##
## Inputs (place in data/ directory):
##   data/List of DEGs for supplementary.xlsx
##
## Outputs (written to plots2/):
##   volcano_2x2.pdf/png
##   heatmaps.pdf
##   DEG_up_down_tables.xlsx
##   go_dotplots.pdf
##   GO_results.xlsx
##   individual_volcanos/  individual_heatmaps/  individual_GO/
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

dir.create("plots2",                          showWarnings = FALSE)
dir.create("plots2/individual_volcanos",      showWarnings = FALSE)
dir.create("plots2/individual_heatmaps",      showWarnings = FALSE)
dir.create("plots2/individual_GO",            showWarnings = FALSE)

## -----------------------------------------------------------------------------
## 1. LOAD THE 4 FULL DEG SHEETS
## -----------------------------------------------------------------------------

# Sheet names in the supplementary Excel file
sheets <- list(
  AB_BSL_MID = " AB BSL-MID full",
  AB_MID_END = "AB MID-END full",
  SQ_BSL_MID = "SQ BSL-MID full",
  SQ_MID_END = "SQ MID-END full"
)

full_data <- lapply(sheets, function(sheet_name) {
  df <- read.xlsx("data/List of DEGs for supplementary.xlsx", sheet = sheet_name)
  
  # Flip LFC: sheets use BSL as reference; we want the later timepoint as reference
  # so positive LFC = upregulated at MID or END
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
})

## -----------------------------------------------------------------------------
## 2. VOLCANO PLOT FUNCTION
## Labels top n significant named genes (no Ensembl IDs)
## x-axis fixed to -8 / +8 with squishing for outliers
## -----------------------------------------------------------------------------

plot_volcano <- function(df, title, n_label = 20) {
  
  # Select top n DEGs to label, ranked by adjusted p-value
  df$delabel   <- NA
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
## 3. GENERATE VOLCANO PLOTS — individual PNGs + combined 2x2 grid
## -----------------------------------------------------------------------------

volcano_titles <- list(
  AB_BSL_MID = "Abdominal — MID vs BSL",
  AB_MID_END = "Abdominal — END vs MID",
  SQ_BSL_MID = "Subcutaneous — MID vs BSL",
  SQ_MID_END = "Subcutaneous — END vs MID"
)

# Save individual PNGs and collect plot objects
p_volc <- list()
for (comparison in names(full_data)) {
  p <- plot_volcano(full_data[[comparison]], title = volcano_titles[[comparison]])
  p_volc[[comparison]] <- p
  ggsave(
    paste0("plots2/individual_volcanos/volcano_", comparison, ".png"),
    p, width = 7, height = 5, dpi = 300
  )
  cat("Saved individual volcano:", comparison, "\n")
}

# Build 2x2 combined grid with depot and comparison labels
blank     <- wrap_elements(full = textGrob("", gp = gpar(fontsize = 1)))
abd_title <- wrap_elements(full = textGrob("Abdominal Adipose Tissue",
                                           gp = gpar(fontface = "bold", fontsize = 12)))
sub_title <- wrap_elements(full = textGrob("Subcutaneous Adipose Tissue",
                                           gp = gpar(fontface = "bold", fontsize = 12)))
row_bm    <- wrap_elements(full = textGrob("MID vs BSL", rot = 90,
                                           gp = gpar(fontface = "bold", fontsize = 11)))
row_me    <- wrap_elements(full = textGrob("END vs MID", rot = 90,
                                           gp = gpar(fontface = "bold", fontsize = 11)))

volcano_grid <-
  (blank  | abd_title            | sub_title)            /
  (row_bm | p_volc$AB_BSL_MID   | p_volc$SQ_BSL_MID)   /
  (row_me | p_volc$AB_MID_END   | p_volc$SQ_MID_END)   +
  plot_layout(heights = c(0.07, 1, 1), widths = c(0.07, 1, 1))

ggsave("plots2/volcano_2x2.pdf", volcano_grid, width = 14, height = 10, dpi = 300)
ggsave("plots2/volcano_2x2.png", volcano_grid, width = 14, height = 10, dpi = 300)
cat("Saved: volcano_2x2\n")

## -----------------------------------------------------------------------------
## 4. HEATMAPS — top 50 DEGs, Z-scored pseudo counts
## -----------------------------------------------------------------------------

tp_colors <- c("Baseline" = "#66C2A5", "Middle" = "#FC8D62", "End" = "#8DA0CB")

# Each config defines a comparison: which DEG sheet, which expression sheet,
# and which sample IDs belong to each group
heatmap_configs <- list(
  list(
    label      = "AB MID vs BSL",
    de_sheet   = " AB BSL-MID full",
    expr_sheet = "Group_1vsGroup_3",
    grp1_pat   = c("SL364592","SL364598","SL364604","SL364610",
                   "SL368255","SL368261","SL368267","SL368273"),
    grp2_pat   = c("SL364596","SL364602","SL364608","SL364614",
                   "SL368259","SL368265","SL368271","SL368277"),
    grp1_name  = "Baseline", grp2_name = "Middle"
  ),
  list(
    label      = "AB END vs MID",
    de_sheet   = "AB MID-END full",
    expr_sheet = "Group_1vsGroup_3",
    grp1_pat   = c("SL364596","SL364602","SL364608","SL364614",
                   "SL368259","SL368265","SL368271","SL368277"),
    grp2_pat   = c("SL364592","SL364598","SL364604","SL364610",
                   "SL368255","SL368261","SL368267","SL368273"),
    grp1_name  = "Middle", grp2_name = "End"
  ),
  list(
    label      = "SQ MID vs BSL",
    de_sheet   = "SQ BSL-MID full",
    expr_sheet = "Group_4vsGroup_6",
    grp1_pat   = c("SL364593","SL364599","SL364605","SL364611",
                   "SL368256","SL368262","SL368268","SL368274"),
    grp2_pat   = c("SL364597","SL364603","SL364609","SL364615",
                   "SL368260","SL368266","SL368272","SL368278"),
    grp1_name  = "Baseline", grp2_name = "Middle"
  ),
  list(
    label      = "SQ END vs MID",
    de_sheet   = "SQ MID-END full",
    expr_sheet = "Group_4vsGroup_6",
    grp1_pat   = c("SL364597","SL364603","SL364609","SL364615",
                   "SL368260","SL368266","SL368272","SL368278"),
    grp2_pat   = c("SL364593","SL364599","SL364605","SL364611",
                   "SL368256","SL368262","SL368268","SL368274"),
    grp1_name  = "Middle", grp2_name = "End"
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
pdf("plots2/heatmaps.pdf", width = 10, height = 12)
for (cfg in heatmap_configs) {
  
  de_df <- read.xlsx("data/List of DEGs for supplementary.xlsx", sheet = cfg$de_sheet)
  de_df$log2FoldChange <- -de_df$log2FoldChange
  de_df$label <- ifelse(!is.na(de_df$gene_name) & de_df$gene_name != "",
                        de_df$gene_name, de_df$Geneid)
  de_df$diffexpressed <- "NO"
  de_df$diffexpressed[de_df$log2FoldChange >  1 & !is.na(de_df$padj) & de_df$padj < 0.05] <- "UP"
  de_df$diffexpressed[de_df$log2FoldChange < -1 & !is.na(de_df$padj) & de_df$padj < 0.05] <- "DOWN"
  
  top50 <- de_df %>%
    filter(diffexpressed != "NO") %>%
    arrange(padj) %>%
    slice_head(n = 50)
  
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
cat("Saved: heatmaps.pdf\n")

## -----------------------------------------------------------------------------
## 5. UP/DOWN DEG TABLES → Excel
## -----------------------------------------------------------------------------

up_down_sheets <- list(
  AB_MID_vs_BSL_UP   = "AB B-M up",
  AB_MID_vs_BSL_DOWN = "AB B-M down",
  AB_END_vs_MID_UP   = "AB M-E up",
  AB_END_vs_MID_DOWN = "AB M-E down",
  SQ_MID_vs_BSL_UP   = "SQ B-M up",
  SQ_MID_vs_BSL_DOWN = "SQ B-M down",
  SQ_END_vs_MID_UP   = "SQ M-E up",
  SQ_END_vs_MID_DOWN = "SQ M-E down"
)

out_cols <- c("Geneid", "gene_name", "Gene.description",
              "baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj")

wb <- createWorkbook()

# Summary tab: gene counts per comparison
summary_rows <- lapply(names(up_down_sheets), function(nm) {
  df <- read.xlsx("data/List of DEGs for supplementary.xlsx", sheet = up_down_sheets[[nm]])
  data.frame(Comparison = nm, N_genes = nrow(df))
})
addWorksheet(wb, "Summary")
writeData(wb, "Summary", do.call(rbind, summary_rows))

# One sheet per UP/DOWN subset
for (nm in names(up_down_sheets)) {
  df <- read.xlsx("data/List of DEGs for supplementary.xlsx", sheet = up_down_sheets[[nm]])
  sheet_name <- substr(nm, 1, 31)
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, df[, out_cols])
}

saveWorkbook(wb, "plots2/DEG_up_down_tables.xlsx", overwrite = TRUE)
cat("Saved: DEG_up_down_tables.xlsx\n")

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

go_configs <- list(
  list(comparison = "AB_MID_vs_BSL", full_sheet = " AB BSL-MID full",
       up_sheet = "AB B-M up", down_sheet = "AB B-M down"),
  list(comparison = "AB_END_vs_MID", full_sheet = "AB MID-END full",
       up_sheet = "AB M-E up", down_sheet = "AB M-E down"),
  list(comparison = "SQ_MID_vs_BSL", full_sheet = "SQ BSL-MID full",
       up_sheet = "SQ B-M up", down_sheet = "SQ B-M down"),
  list(comparison = "SQ_END_vs_MID", full_sheet = "SQ MID-END full",
       up_sheet = "SQ M-E up", down_sheet = "SQ M-E down")
)

go_results <- list()
for (cfg in go_configs) {
  cat("Running GO enrichment:", cfg$comparison, "\n")
  
  full_df  <- read.xlsx("data/List of DEGs for supplementary.xlsx", sheet = cfg$full_sheet)
  up_df    <- read.xlsx("data/List of DEGs for supplementary.xlsx", sheet = cfg$up_sheet)
  dn_df    <- read.xlsx("data/List of DEGs for supplementary.xlsx", sheet = cfg$down_sheet)
  universe <- full_df$gene_name
  
  go_results[[paste0(cfg$comparison, "_UP")]]   <- run_go(up_df$gene_name, universe)
  go_results[[paste0(cfg$comparison, "_DOWN")]] <- run_go(dn_df$gene_name, universe)
  go_results[[paste0(cfg$comparison, "_both")]] <- run_go(
    c(up_df$gene_name, dn_df$gene_name), universe
  )
}

# Save dotplots to a single PDF
pdf("plots2/go_dotplots.pdf", width = 10, height = 8)
for (nm in names(go_results)) {
  if (!is.null(go_results[[nm]]) && nrow(go_results[[nm]]) > 0) {
    print(
      dotplot(go_results[[nm]], showCategory = 20, font.size = 9,
              title = gsub("_", " ", nm)) +
        theme(plot.title = element_text(face = "bold"))
    )
  }
}
dev.off()

# Save individual PNG per comparison
for (nm in names(go_results)) {
  if (!is.null(go_results[[nm]]) && nrow(go_results[[nm]]) > 0) {
    p_go <- dotplot(go_results[[nm]], showCategory = 20, font.size = 9,
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
wb_go <- createWorkbook()
for (nm in names(go_results)) {
  if (!is.null(go_results[[nm]]) && nrow(go_results[[nm]]) > 0) {
    sheet_name <- substr(paste0("GO_", nm), 1, 31)
    addWorksheet(wb_go, sheet_name)
    writeData(wb_go, sheet_name,
              as.data.frame(go_results[[nm]]) %>%
                dplyr::select(all_of(go_out_cols)) %>%
                arrange(p.adjust))
  }
}
saveWorkbook(wb_go, "plots2/GO_results.xlsx", overwrite = TRUE)
cat("Saved: go_dotplots.pdf + GO_results.xlsx\n")

cat("\nAll outputs written to plots2/\n")
