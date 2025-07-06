library(Matrix)
library(DropletUtils)
library(ggplot2)
library(scales)
library(Seurat)
library(tidyverse)

# Load atlas
setwd("C:/Users/Олеся/Desktop/graduation/data_sets/")

rc.integrated <- readRDS("Shahan_epidermis_020625.rds")
# Load reference expression profiles from Bhosale et al. 2018, The Plant Cell 
load(file="endo_exp.RD")

# Reference expression profiles for ploidy
endo_exp[1:10,]

# Extract matrix of SCTransformed expression value (due to a lagre size they been split)
rc1 <- as.matrix(rc.integrated@assays$SCT@data[,1:10000])
rc2 <- as.matrix(rc.integrated@assays$SCT@data[,10001:20000])
rc3 <- as.matrix(rc.integrated@assays$SCT@data[,20001:25054])

rc <- cbind(rc1, rc2)
rc <- cbind(rc, rc3)
str(rc)

# Merge the reference expression profile with the normalized expression matrix of our sample  

merge.rownames <- function (x,y){
  dat <- merge(x = x, y = y, by = "row.names")
  rownames(dat) <- dat$Row.names
  dat <- dat[,-1]
  return(dat)
}

ploidy <- Reduce(merge.rownames, list(endo_exp,rc))

# Prepare customized label name (optional)
ploidy_label=c("2C", "4C", "8C", "16C")

ploidy[,1:10]

#Calculating the correlation coefficient of each cell to each reference expression profile and annotate the cell as the label that it has the highest correlation coefficient with.  
ploidy_stat <- suppressWarnings(sapply(5:ncol(ploidy), function(i) sapply(1:4, function(j) cor.test(ploidy[,i],ploidy[,j],method = "pearson")[c(3,4)])))
ploidy_cor <- ploidy_stat[seq(2,nrow(ploidy_stat),2),]
ploidy_pvalue <- ploidy_stat[seq(1,nrow(ploidy_stat)-1,2),]
ploidy_max <- sapply(1:(ncol(ploidy)-4), function(i) max(as.numeric(ploidy_cor[,i])))
ploidy_ident <- sapply(1:(ncol(ploidy)-4), function(i) ploidy_label[which(as.numeric(ploidy_cor[,i])==max(as.numeric(ploidy_cor[,i])))])
ploidy_maxp <- sapply(1:(ncol(ploidy)-4), function(i) as.numeric(ploidy_pvalue[,i])[which(as.numeric(ploidy_cor[,i])==max(as.numeric(ploidy_cor[,i])))])
names(ploidy_max) <- ploidy_ident


#Store the annotation, correlation coefficient and the p-value in Seurat object
rc.integrated@meta.data$ploidy.ID <- as.character(ploidy_ident)


# In case there is cell with insufficient information for annotation, label them as "unknown"
rc.integrated@meta.data$ploidy.ID[which(rc.integrated@meta.data$ploidy.ID=='character(0)')]="unknown"

view(rc.integrated@meta.data)

options(repr.plot.width=10, repr.plot.height=8)
order <- c("2C","4C","8C","16C","unknown")
palette <- c("#DCEDC8","#42B3D5","#FDEA6F","#CF4F29","#cccccc")

DimPlot(rc.integrated, group.by="ploidy.ID", cols=palette)

saveRDS(rc.integrated, "shahan_ploidy.ID")


#applying ploidy.ID to the set with calculated slingshot 
metadata<-rc.integrated@meta.data["ploidy.ID"]
#load the set
path <- "C:/Users/Олеся/Desktop/graduation/data_sets/shahan_sling_0306"
seurat_obj<- readRDS("C:/Users/Олеся/Desktop/graduation/data_sets/shahan_reclustered_030625.rds")

#find intersect between sets
common_cells <- intersect(Cells(seurat_obj), rownames(metadata))

# Subset metadata and new Seurat object to only matching cells
metadata <- metadata[common_cells, , drop = FALSE]
seurat_obj <- subset(new_seurat_obj, cells = common_cells)

# Add metadata to new Seurat object
seurat_obj <- AddMetaData(seurat_obj, metadata = metadata$ploidy.ID, col.name = "ploidy.ID")
DimPlot(seurat_obj, group.by="ploidy.ID", cols=palette)
#save the file
saveRDS(seurat_obj, file = "shahan_with_sling&ploidy_0706.rds")
