library(dplyr)
library(Seurat)
library(patchwork)
library(tidyverse)
library(SeuratObject)
library(here)
library(readxl)

path <- "C:/Users/Олеся/Desktop/graduation/data_sets"
seurat_obj <- readRDS(paste0(path, "/shahan_with_sling&ploidy_0706.rds"))
M<- seurat_obj@meta.data
DefaultAssay(AT1) <- "RNA"
seurat_obj<- NormalizeData(seurat_obj)

#select the cells you are intrested in (per lineage)
cells_to_keepAT2 <- rownames(seurat_obj@meta.data)[!is.na(seurat_obj@meta.data$sling_lineage1)]
cells_to_keepT2 <- rownames(seurat_obj@meta.data)[!is.na(seurat_obj@meta.data$sling_lineage2)]
cells_to_keepAT1 <- rownames(seurat_obj@meta.data)[!is.na(seurat_obj@meta.data$sling_lineage4)]
cells_to_keepT1 <- rownames(seurat_obj@meta.data)[!is.na(seurat_obj@meta.data$sling_lineage6)]
# Subset the Seurat object using those cells
AT2 <- subset(seurat_obj, cells = cells_to_keepAT2)
T2 <- subset(seurat_obj, cells = cells_to_keepT2)
AT1 <- subset(seurat_obj, cells = cells_to_keepAT1)
T1 <- subset(seurat_obj, cells = cells_to_keepT1)

# Remove unwanted metadata columns
AT2@meta.data <- AT2@meta.data[, !colnames(AT2@meta.data) %in% c(
  "sling_lineage2", "sling_lineage3", "sling_lineage4",
  "sling_lineage5", "sling_lineage6", "sling_lineage7",
  "sling_lineage8", "sling_lineage9"
)]

T2@meta.data <- T2@meta.data[, !colnames(T2@meta.data) %in% c(
  "sling_lineage1", "sling_lineage3", "sling_lineage4",
  "sling_lineage5", "sling_lineage6", "sling_lineage7",
  "sling_lineage8", "sling_lineage9"
)]
AT1@meta.data <- AT1@meta.data[, !colnames(AT1@meta.data) %in% c(
  "sling_lineage2", "sling_lineage3", "sling_lineage1",
  "sling_lineage5", "sling_lineage6", "sling_lineage7",
  "sling_lineage8", "sling_lineage9"
)]
T1@meta.data <- T1@meta.data[, !colnames(T1@meta.data) %in% c(
  "sling_lineage2", "sling_lineage3", "sling_lineage4",
  "sling_lineage5", "sling_lineage1", "sling_lineage7",
  "sling_lineage8", "sling_lineage9"
)]
# Rename 'sling_lineage1' to 'sling_pseudo'
colnames(AT2@meta.data)[colnames(AT2@meta.data) == "sling_lineage1"] <- "sling_pseudo"
colnames(T2@meta.data)[colnames(T2@meta.data) == "sling_lineage2"] <- "sling_pseudo"
colnames(AT1@meta.data)[colnames(AT1@meta.data) == "sling_lineage4"] <- "sling_pseudo"
colnames(T1@meta.data)[colnames(T1@meta.data) == "sling_lineage6"] <- "sling_pseudo"

# Define breakpoints and labels
breaksAT2 <- c(-Inf, 40, 60, 80, 100, 110, Inf)
labelsAT2 <- c("Atrichoblast2_m1", "Atrichoblast2_m2", "Atrichoblast2_t",
            "Atrichoblast2_e1", "Atrichoblast2_e2", "Atrichoblast2_d")
breaksT2 <- c(-Inf, 30, 60, 90, 110, 120, Inf)
labelsT2 <- c("Trichoblast2_m1", "Trichoblast2_m2", "Trichoblast2_t",
               "Trichoblast2_e1", "Trichoblast2_e2", "Trichoblast2_d")

breaksAT1 <- c(-Inf, 20, 60, 80, 100, 120, Inf)
labelsAT1 <- c("Atrichoblast1_m1", "Atrichoblast1_m2", "Atrichoblast1_t",
               "Atrichoblast1_e1", "Atrichoblast1_e2", "Atrichoblast1_d")

breaksT1 <- c(-Inf, 20, 60, 80, 90, 110, Inf)
labelsT1 <- c("Trichoblast1_m1", "Trichoblast1_m2", "Trichoblast1_t",
              "Trichoblast1_e1", "Trichoblast1_e2", "Trichoblast1_d")


# Add new annotation column
AT2$Cell_type<- cut(AT2@meta.data$sling_pseudo, breaks = breaksAT2, labels = labelsAT2, right = TRUE)
view(AT2@meta.data)
T2$Cell_type<- cut(T2@meta.data$sling_pseudo, breaks = breaksT2, labels = labelsT2, right = TRUE)
view(T2@meta.data)
AT1$Cell_type<- cut(AT1@meta.data$sling_pseudo, breaks = breaksAT1, labels = labelsAT1, right = TRUE)
view(AT1@meta.data)
T1$Cell_type<- cut(T1@meta.data$sling_pseudo, breaks = breaksT1, labels = labelsT1, right = TRUE)
view(T1@meta.data)

saveRDS(AT2, file = "C:/Users/Олеся/Desktop/graduation/data_sets/realistic_layout/AT2.rds")
saveRDS(T2, file = "C:/Users/Олеся/Desktop/graduation/data_sets/realistic_layout/T2.rds")
saveRDS(AT1, file = "C:/Users/Олеся/Desktop/graduation/data_sets/realistic_layout/AT1.rds")
saveRDS(T1, file = "C:/Users/Олеся/Desktop/graduation/data_sets/realistic_layout/T1.rds")
