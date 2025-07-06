library(Seurat)
library(tradeSeq)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(ggplot2)
library(openxlsx)
library(grid)


#load data
path <- "C:/Users/Олеся/Desktop/graduation/data_sets"
seurat_obj <- readRDS(paste0(path, "/shahan_with_sling&ploidy_0706.rds"))
m<-seurat_obj@meta.data

#DimPlot(shahan, group.by = "celltype.anno.crude" )
#Select one of the subsets
#seurat_obj <- subset(shahan, subset = celltype.anno.crude %in% "Atrichoblast")
#seurat_obj <- subset(shahan, subset = celltype.anno.crude %in% "Trichoblast")


#markers for AT and t
marker_data <- read.xlsx(paste0(path,"/AT_T_markers.xlsx"))
AT_marker <- marker_data[marker_data$Cell_type == "Atrichoblast", ]
T_marker <- marker_data[marker_data$Cell_type == "Trichoblast", ]
Ep_marker <- marker_data[marker_data$Cell_type == "Epidermis", ]


genes_to_plot_AT <-as.character(AT_marker$TAIR_ID)
genes_to_plot_T <-as.character(T_marker$TAIR_ID)
genes_to_plot_Ep <-as.character(Ep_marker$TAIR_ID)



DefaultAssay(seurat_obj) <- "RNA"
seurat_obj <- NormalizeData(seurat_obj)
exprs_mat <- GetAssayData(seurat_obj, slot = "data")


##rerun following part until end per each lineage ##
#------------------------------------------------------
# Extract normalized expression matrix (genes x cells)
exprs_mat <- GetAssayData(seurat_obj, slot = "data", assay = "RNA")

# Extract pseudotime values from metadata
pt <- seurat_obj$sling_lineage1
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

gene_labels_AT <- setNames(as.character(AT_marker$gene_labels), AT_marker$TAIR_ID)
gene_labels_T  <- setNames(as.character(T_marker$gene_labels), T_marker$TAIR_ID)
gene_labels_Ep  <- setNames(as.character(Ep_marker$gene_labels), Ep_marker$TAIR_ID)

heatmap_AT <- plot_pseudotime_heatmap(genes_to_plot_AT, exprs_mat, pt, gene_labels_AT, "Atrichoblast Markers for lineage1")
heatmap_T  <- plot_pseudotime_heatmap(genes_to_plot_T, exprs_mat, pt, gene_labels_T, "Trichoblast Marker for lineage1")
heatmap_Ep  <- plot_pseudotime_heatmap(genes_to_plot_Ep, exprs_mat, pt, gene_labels_Ep, "Epidermis Marker for lineage1")

pt_ordered <- pt[order(pt)]

draw(heatmap_AT)
add_pseudotime_ticks(pt_ordered)

draw(heatmap_T)
add_pseudotime_ticks(pt_ordered)

draw(heatmap_Ep)
add_pseudotime_ticks(pt_ordered)

setwd("C:/Users/Олеся/Desktop/graduation/data_sets//shahan_sling_0306/HMs/")
png("heatmap_all_ATmarkers_lineage1.png", width = 1200, height = 900, res = 150)  # size and resolution adjustable
draw(heatmap_AT)
add_pseudotime_ticks(pt_ordered)
dev.off()

png("heatmap_all_Tmarkers_lineage1.png", width = 1200, height = 900, res = 150)
draw(heatmap_T)
add_pseudotime_ticks(pt_ordered)
dev.off()

png("heatmap_all_Epidermis_lineage1.png", width = 1200, height = 900, res = 150)
draw(heatmap_Ep)
add_pseudotime_ticks(pt_ordered)
dev.off()
