library(tradeSeq)
library(scater)
library(slingshot)
library(SingleCellExperiment)
library(Seurat)
library(tidyverse)
library(ggplot2)
library(CytoTRACE)

#open AT_T data set
path <- "C:/Users/Олеся/Desktop/graduation/data_sets/shahan_sling_0306"
AT_T_data <- readRDS("C:/Users/Олеся/Desktop/graduation/data_sets/shahan_reclustered_030625.rds")

dp <- DimPlot(AT_T_data, label=TRUE, , label.size = 5)
dp
#ggsave(filename = paste0(path, "/R_plots/umap_epidermis_4.06.png"), plot = dp, width = 4, height = 3.5, dpi = 150)

colors <- c('#fef0d9','#fdd49e','#fc8d59','#b30000')

# perform cytotrace
mat <- GetAssayData(object = AT_T_data, assay = "RNA", slot = "counts")

as_matrix <- function(mat){
  tmp <- matrix(data=0L, nrow = mat@Dim[1], ncol = mat@Dim[2])
  row_pos <- mat@i+1
  col_pos <- findInterval(seq(mat@x)-1,mat@p[-1])+1
  val <- mat@x
  
  for (i in seq_along(val)){
    tmp[row_pos[i],col_pos[i]] <- val[i]
  }
  row.names(tmp) <- mat@Dimnames[[1]]
  colnames(tmp) <- mat@Dimnames[[2]]
  return(tmp)
}

mat2 <- as_matrix(mat)

cyto_ra <- CytoTRACE(mat2,ncores = 1,enableFast = T)

cyto_ra$CytoTRACE <- 1-cyto_ra$CytoTRACE

AT_T_data <- AddMetaData(
  object = AT_T_data,
  metadata = cyto_ra$CytoTRACE,
  col.name = "differentiation_level"
)
meta<-AT_T_data@meta.data
FeaturePlot(AT_T_data, features = c('differentiation_level'))+
  scale_color_gradientn(colors = colors)+theme(text=element_text(size=10))+NoAxes()

ggsave(filename = paste0(path, "/R_plots/cytotrace_0606.png"), plot = last_plot(), width = 4, height = 3.5, dpi = 150)

# slingshot for full set together for multiple curves

sce_T.sling <- as.SingleCellExperiment(AT_T_data)

# Run slingshot to infer multiple trajectories
# Specify the starting cluster 
sce_T.sling <- slingshot(sce_T.sling,
                         cluster = sce_T.sling$seurat_clusters,
                         start.clus = "18",
                         reducedDim = 'PCA') 

# Check the inferred pseudotime for multiple lineages
head(slingPseudotime(sce_T.sling))

# Embed the curves onto your reduced dimension (UMAP in this case)
embedded_T <- embedCurves(sce_T.sling, "UMAP") # Assuming your UMAP is stored with this name

# Get the UMAP coordinates
umap_df <- as.data.frame(reducedDim(sce_T.sling, "UMAP"))

# Start the plot with the UMAP embedding (simple version)
p <- ggplot(umap_df, aes(x = umap_1, y = umap_2)) +
  geom_point(size = 0.5, color = "lightgray") +
  theme_bw()

# Add the Slingshot trajectories as separate paths
for (i in seq_along(slingCurves(embedded_T))) {
  curve_data <- as.data.frame(slingCurves(embedded_T)[[i]]$s[slingCurves(embedded_T)[[i]]$ord, ])
  if (nrow(curve_data) > 1) {
    p <- p + geom_path(data = curve_data, aes(x = umap_1, y = umap_2),
                       linewidth = 1.2, color = "black")
  }
}

print(p)
ggsave( filename = paste0(path, "/R_plots/sling_multi(18)_06.06.png"),
        plot = p,
        width = 4,
        height = 3.5,
        dpi = 150)

# selecting trajectories
slingLineages(sce_T.sling) # cheaking lineage names and clusters
pseudotime_lineage1 <- slingPseudotime(sce_T.sling)[, "Lineage1"]
pseudotime_lineage2 <- slingPseudotime(sce_T.sling)[, "Lineage2"]
pseudotime_lineage3 <- slingPseudotime(sce_T.sling)[, "Lineage3"]
pseudotime_lineage4 <- slingPseudotime(sce_T.sling)[, "Lineage4"]
pseudotime_lineage5 <- slingPseudotime(sce_T.sling)[, "Lineage5"]
pseudotime_lineage6 <- slingPseudotime(sce_T.sling)[, "Lineage6"]
pseudotime_lineage7 <- slingPseudotime(sce_T.sling)[, "Lineage7"]
pseudotime_lineage8 <- slingPseudotime(sce_T.sling)[, "Lineage8"]
pseudotime_lineage9 <- slingPseudotime(sce_T.sling)[, "Lineage9"]

# saving lineages to AT_T_data
AT_T_data <- AddMetaData(AT_T_data, metadata = pseudotime_lineage1, col.name = "sling_lineage1")
AT_T_data <- AddMetaData(AT_T_data, metadata = pseudotime_lineage2, col.name = "sling_lineage2")
AT_T_data <- AddMetaData(AT_T_data, metadata = pseudotime_lineage3, col.name = "sling_lineage3")
AT_T_data <- AddMetaData(AT_T_data, metadata = pseudotime_lineage4, col.name = "sling_lineage4")
AT_T_data <- AddMetaData(AT_T_data, metadata = pseudotime_lineage5, col.name = "sling_lineage5")
AT_T_data <- AddMetaData(AT_T_data, metadata = pseudotime_lineage6, col.name = "sling_lineage6")
AT_T_data <- AddMetaData(AT_T_data, metadata = pseudotime_lineage7, col.name = "sling_lineage7")
AT_T_data <- AddMetaData(AT_T_data, metadata = pseudotime_lineage8, col.name = "sling_lineage8")
AT_T_data <- AddMetaData(AT_T_data, metadata = pseudotime_lineage9, col.name = "sling_lineage9")

#combining lineages together(optional)
#AT_T_data$Combined_lineage_56 <- ifelse(!is.na(AT_T_data$sling_lineage5_T),
#                                        AT_T_data$sling_lineage5_T,
#                                        AT_T_data$sling_lineage6_T)

#AT_T_data$Combined_lineage_356 <- ifelse(!is.na(AT_T_data$sling_lineage3_T),
#                                         AT_T_data$sling_lineage3_T,
#                                         ifelse(!is.na(AT_T_data$sling_lineage5_T),
#                                                AT_T_data$sling_lineage5_T,
#                                                AT_T_data$sling_lineage6_T))

#AT_T_data$Combined_lineage_37 <- ifelse(!is.na(AT_T_data$sling_lineage3_T),
#                                        AT_T_data$sling_lineage3_T,
#                                        AT_T_data$sling_lineage7_AT)
meta_shahan <- AT_T_data@meta.data

#saving the rds file (it was given SCT related error, so i tried to clean it)
AT_T_clean <- CreateSeuratObject(counts = GetAssayData(AT_T_data, assay = "RNA", slot = "counts"))

AT_T_clean@meta.data <- AT_T_data@meta.data
AT_T_clean@reductions$umap <- AT_T_data@reductions$umap
AT_T_clean@active.ident <- AT_T_data@active.ident
saveRDS(AT_T_clean, file = paste0(path, "/shahan_with_sling_0606.rds"))


# plot saved pseudotimes
umap_df <- as.data.frame(Embeddings(AT_T_data, "umap"))
umap_df$pseudotime <- AT_T_data$sling_lineage9 # select a lineage to plot

ggplot(umap_df, aes(x = umap_1, y = umap_2, color = pseudotime)) +
  geom_point(size = 1.2) +
  scale_color_viridis_c() +
  theme_minimal() +
  ggtitle("UMAP with Lineage 9 Pseudotime")

ggsave( filename = paste0(path, "/R_plots/sling_lineage_9_21.05.png"),
       width = 4,
       height = 3.5,
       dpi = 150)



