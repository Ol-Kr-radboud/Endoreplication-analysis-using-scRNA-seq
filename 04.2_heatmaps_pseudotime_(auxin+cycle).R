library(Seurat)
library(tradeSeq)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(ggplot2)
library(openxlsx)
library(grid)

#open RDS file
path <- "C:/Users/Олеся/Desktop/graduation/data_sets"
seurat_obj <- readRDS(paste0(path, "/shahan_with_sling&ploidy_0706.rds"))


#markers for aux and cycle 
Aux_data <- read.xlsx(paste0(path,"/hormones_gene_annotation.xlsx"), sheet = "Auxin")
cycle_data <- read.xlsx(paste0(path,"/Cell_cycle_related.xlsx"))

#since gene lists are large they were separated into two per set
aux_marker1 <-Aux_data[Aux_data$Type %in% c("synthesis", "auxin response"),]
aux_marker2 <-Aux_data[Aux_data$Type %in% c("conjugation/deconjugation", "transport"),]
cycle_marker1 <-cycle_data[cycle_data$gene_groups %in% c("CDK", "CDK_related", "CYC"),]
cycle_marker2 <-cycle_data[!(cycle_data$gene_groups %in% c("CDK", "CDK_related", "CYC")),]

genes_to_plot_aux1 <- as.character(aux_marker1$TAIR_ID)
genes_to_plot_aux2 <- as.character(aux_marker2$TAIR_ID)
genes_to_plot_cycle1 <- as.character(cycle_marker1$TAIR_ID)
genes_to_plot_cycle2 <- as.character(cycle_marker2$TAIR_ID)


DefaultAssay(seurat_obj) <- "RNA"
seurat_obj <- NormalizeData(seurat_obj)
exprs_mat <- GetAssayData(seurat_obj, slot = "data")


##rerun following part until end per each lineage ##
#------------------------------------------------------

# Extract normalized expression matrix (genes x cells)
exprs_mat <- GetAssayData(seurat_obj, slot = "data", assay = "RNA")

# Extract pseudotime values from metadata
pt <- seurat_obj$sling_lineage9
names(pt) <- colnames(seurat_obj)  # ensure names match column names of exprs_mat

# Filter to valid pseudotime cells
valid_cells <- names(pt)[!is.na(pt)]
pt <- pt[valid_cells]
exprs_mat <- exprs_mat[, valid_cells]

# -------------------------
plot_pseudotime_heatmap <- function(genes_to_plot, exprs_mat, pt, gene_labels, title = "Heatmap") {
  genes_present <- intersect(genes_to_plot, rownames(exprs_mat))
  exprs_sub <- exprs_mat[genes_present, , drop = FALSE]
  
  order_cells <- order(pt)
  exprs_ordered <- exprs_sub[, order_cells, drop = FALSE]
  pt_ordered <- pt[order_cells]
  
  exprs_scaled <- t(apply(exprs_ordered, 1, function(x) {
    x <- pmax(x, 0)
    if(max(x) == 0) return(x)
    (x - min(x)) / (max(x) - min(x))
  }))
  
  gene_names <- sapply(rownames(exprs_scaled), function(gene) {
    if (!is.na(gene_labels[gene]) && gene_labels[gene] != "") {
      return(gene_labels[gene])
    } else {
      return(gene)
    }
  })
  rownames(exprs_scaled) <- gene_names
  
  heatmap_colors <- colorRampPalette(c("#ece7f2", "#a6bddb", "#2b8cbe"))(100)
  pt_col_fun <- colorRamp2(range(pt_ordered), c("#fef0d9","#e34a33"))
  
  # Calculate cell indices for pseudotime breaks at steps of 10
  
  # Pseudotime color bar annotation (first track)
  pt_anno <- HeatmapAnnotation(
    Pseudotime = pt_ordered,
    col = list(Pseudotime = pt_col_fun),
    annotation_name_side = "left",
    height = unit(8, "mm")
  )
  # Create pseudotime color bar (continuous)
  pt_col_fun <- colorRamp2(range(pt_ordered), c("#fef0d9","#e34a33"))
  column_ha <- HeatmapAnnotation(
    Pseudotime = pt_ordered,
    col = list(Pseudotime = pt_col_fun),
    annotation_name_side = "left"
  )
  
  # Draw heatmap
  ht<- Heatmap(
    exprs_scaled,
    name = "Expression",
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    show_column_names = FALSE,
    column_title = title,
    row_title = "Genes",
    col = heatmap_colors,
    top_annotation = column_ha
  )
  return(ht)
}

add_pseudotime_ticks <- function(pt_ordered, step = 20, heatmap_name = "Expression") {
  breaks <- seq(0, max(pt_ordered), by = step)
  n_cols <- length(pt_ordered)
  
  decorate_heatmap_body(heatmap_name, {
    for (b in breaks) {
      idx <- which.min(abs(pt_ordered - b))
      x_pos <- idx / n_cols  # normalize index to [0,1]
      
      grid.lines(x = unit(c(x_pos, x_pos), "npc"),
                 y = unit(c(0, 1), "npc"),
                 gp = gpar(col = "black", lty = 2, lwd = 1))
    }
  })
}

gene_labels_aux1 <- setNames(as.character(aux_marker1$gene_labels), aux_marker1$TAIR_ID)
gene_labels_aux2 <- setNames(as.character(aux_marker2$gene_labels), aux_marker2$TAIR_ID)
gene_labels_cycle1 <- setNames(as.character(cycle_marker1$gene_labels), cycle_marker1$TAIR_ID)
gene_labels_cycle2 <- setNames(as.character(cycle_marker2$gene_labels), cycle_marker2$TAIR_ID)


heatmap_aux1 <- plot_pseudotime_heatmap(genes_to_plot_aux1, exprs_mat, pt, gene_labels_aux1, "Auxin related Markers for lineage9, graph 1")
heatmap_aux2 <- plot_pseudotime_heatmap(genes_to_plot_aux2, exprs_mat, pt, gene_labels_aux2, "Auxin related Markers for lineage9, graph 2")
heatmap_cycle1 <- plot_pseudotime_heatmap(genes_to_plot_cycle1, exprs_mat, pt, gene_labels_cycle1, "Cell cylcle related Markers for lineage9, grahp 1")
heatmap_cycle2 <- plot_pseudotime_heatmap(genes_to_plot_cycle2, exprs_mat, pt, gene_labels_cycle2, "Cell cylcle related Markers for lineage9, grahp 2")


pt_ordered <- pt[order(pt)]

draw(heatmap_aux1)
add_pseudotime_ticks(pt_ordered)

draw(heatmap_aux2)
add_pseudotime_ticks(pt_ordered)

draw(heatmap_cycle1)
add_pseudotime_ticks(pt_ordered)

draw(heatmap_cycle2)
add_pseudotime_ticks(pt_ordered)

setwd("C:/Users/Олеся/Desktop/graduation/data_sets//shahan_sling_0306/HMs/")
png("heatmap_aux1_lineage9.png", width = 3600, height = 4500, res = 300)  # size and resolution adjustable
draw(heatmap_aux1)
add_pseudotime_ticks(pt_ordered)
dev.off()

setwd("C:/Users/Олеся/Desktop/graduation/data_sets//shahan_sling_0306/HMs/")
png("heatmap_aux2_lineage9.png", width = 3600, height = 4500, res = 300)  # size and resolution adjustable
draw(heatmap_aux2)
add_pseudotime_ticks(pt_ordered)
dev.off()

png("heatmap_cycle1_lineage9.png", width = 3600, height = 4500, res = 300)
draw(heatmap_cycle1)
add_pseudotime_ticks(pt_ordered)
dev.off()

png("heatmap_cycle2_lineage9.png", width = 3600, height = 4500, res = 300)
draw(heatmap_cycle2)
add_pseudotime_ticks(pt_ordered)
dev.off()
