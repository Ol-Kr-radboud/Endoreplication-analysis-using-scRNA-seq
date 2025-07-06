library(tidyverse)
library(Seurat)
library(here)
library(readxl)

set.seed(100)
path <- here()
#in_path <- paste0(path, "/in")
in_path <- "E:/R_images/RootCellAtlas/Shahan2020"

out_path <- paste0(path, "/out")
seurat_obj <- readRDS(file = paste0(in_path, "/GSE152766_Epidermis_LRC_Atlas_Shahan.rds"))

DimPlot(seurat_obj)

meta.data <- seurat_obj@meta.data

#exploring original annotation ----

class(seurat_obj@meta.data$celltype.ID)
if (is.list(seurat_obj@meta.data$celltype.ID)) {
  seurat_obj@meta.data$celltype.ID <- as.character(unlist(seurat_obj@meta.data$celltype.ID))
}
levels(seurat_obj@meta.data$celltype.ID)
seurat_obj@meta.data$celltype.ID <- factor(seurat_obj@meta.data$celltype.ID)
table(seurat_obj@meta.data$celltype.ID)

table(seurat_obj@meta.data$final.anno.P)

DimPlot(seurat_obj, group.by = "celltype.ID")

DimPlot(
  object = seurat_obj,
  reduction = "umap",       # or "tsne", "pca"
  group.by = "celltype.ID",
  label = TRUE,             # optionally label clusters
  label.size = 4
) + ggtitle("UMAP colored by Celltype.ID")

DimPlot(
  object = seurat_obj,
  reduction = "umap",       # or "tsne", "pca"
  group.by = "final.anno.P",
  label = TRUE,             # optionally label clusters
  label.size = 4
) + ggtitle("UMAP colored by Celltype.ID")

keep_clusters <- c("hair cells", "LRC & non-hair cells", "non-hair cells", "QC cells")

# Subset Seurat object
subset_obj_small <- subset(seurat_obj, subset = celltype.ID %in% keep_clusters)

DimPlot(subset_obj_small, group.by = "celltype.ID")
rm(subset_obj_small)
#conclusion: original annotation is far from perfect



#Subsetting the data ----
DimPlot(seurat_obj, label = TRUE,label.size = 8)

levels(Idents(seurat_obj))

keep_clusters <- c("16", "8", "18", "15", "11", "12", "17", "9", "13","0", "2", "6")

# Subset Seurat object
# subset_obj_small <- subset(seurat_obj, idents = keep_clusters)
# DimPlot(subset_obj_small, label = TRUE,label.size = 8)

seurat_obj$cluster <- Idents(seurat_obj)
unique(seurat_obj$cluster)

cells_to_keep <- rownames(seurat_obj@meta.data[seurat_obj$cluster %in% keep_clusters, ])

# Subset using cells vector (avoids re-evaluating large matrices)
seurat_obj <- subset(seurat_obj, cells = cells_to_keep)
DimPlot(seurat_obj, label = TRUE,label.size = 8)


#Data reclustering:

seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)
seurat_obj <- FindClusters(seurat_obj, resolution = 1.5)  # adjust as needed
seurat_obj <- RunUMAP(seurat_obj, dims = 1:20)
p2<- DimPlot(seurat_obj, reduction = "umap", label = TRUE, label.size = 5)
p2
ggsave(paste0(out_path, "/umap_epidermis_1.pdf"), plot = p2, width = 10, height = 12)


#Checking marker genes
markers_df<- read_excel("E:/R_images/RootCellAtlas/Markers_13.12.2021.xlsx")

gene_ids <- markers_df$GeneID
gene_names <- markers_df$GeneName

# Filter genes that exist in Seurat object
valid_genes <- gene_ids %in% rownames(seurat_obj)
gene_ids <- gene_ids[valid_genes]
gene_names <- gene_names[valid_genes]

# Create mapping: "GeneID (GeneName)"
gene_labels <- paste0(gene_ids, " (", gene_names, ")")
names(gene_labels) <- gene_ids

# DoHeatmap + custom y-axis labels
p <- DoHeatmap(seurat_obj, features = gene_ids) +
  scale_y_discrete(labels = gene_labels[gene_ids]) +
  ggtitle("Marker genes in epidermis")

p

ggsave(paste0(out_path, "/heatmap_marker_genes_2.pdf"), plot = p, width = 10, height = 12)

#vizualizing cell cycle markers
markers_df<- read.csv("E:/Datasets/!Datasets/Cell_cycle_related.csv")
gene_ids <- markers_df$Cell.cycle.genes
gene_names <- markers_df$X.1

valid_genes <- gene_ids %in% rownames(seurat_obj)
gene_ids <- gene_ids[valid_genes]
gene_names <- gene_names[valid_genes]
gene_labels <- paste0(gene_ids, " (", gene_names, ")")
names(gene_labels) <- gene_ids

# DoHeatmap + custom y-axis labels
p <- DoHeatmap(seurat_obj, features = gene_ids) +
  scale_y_discrete(labels = gene_labels[gene_ids]) +
  ggtitle("cell cycle genes in epidermis")

p

ggsave(paste0(out_path, "/heatmap_cell_cycle_genes.pdf"), plot = p, width = 10, height = 12)

##vizualizing LRP markers
markers_df<- read_excel("E:/R_images/RootCellAtlas/Markers_LRP.xlsx")
gene_ids <- markers_df$Tair.ID
gene_names <- markers_df$Gene

valid_genes <- gene_ids %in% rownames(seurat_obj)
gene_ids <- gene_ids[valid_genes]
gene_names <- gene_names[valid_genes]
gene_labels <- paste0(gene_ids, " (", gene_names, ")")
names(gene_labels) <- gene_ids

# DoHeatmap + custom y-axis labels
p <- DoHeatmap(seurat_obj, features = gene_ids) +
  scale_y_discrete(labels = gene_labels[gene_ids]) +
  ggtitle("LRP genes in epidermis")

#no LRP marker genes are present

#Analysis suggest that most of the clusters are in place, but clusters 34 (FEZ positive dividing LRC) and 31 (FEZ positive differentiating LRC) need to be removed
exclude_clusters <- c("31", "34")

# Identify cells to keep
cells_to_keep <- WhichCells(seurat_obj, idents = setdiff(levels(Idents(seurat_obj)), exclude_clusters))

# Subset the object safely
seurat_obj <- subset(seurat_obj, cells = cells_to_keep)
DimPlot(seurat_obj, label = TRUE,label.size = 8)

#Second round of reclustering

seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)
seurat_obj <- FindClusters(seurat_obj, resolution = 1.5)  # adjust as needed
seurat_obj <- RunUMAP(seurat_obj, dims = 1:20)
DimPlot(seurat_obj, reduction = "umap", label = TRUE, label.size = 5)

# redo heatmaps for marker genes

saveRDS(seurat_obj, file = "E:/R_images/RootCellAtlas/Shahan2020/Shahan_epidermis_020625.rds")

DimPlot(seurat_obj, split.by = "orig.ident", label = TRUE)
p <- DimPlot(seurat_obj, split.by = "orig.ident", label = TRUE)

ggsave(paste0(out_path, "/DimPlot_different_samples.pdf"), plot = p, width = 20, height = 12)


