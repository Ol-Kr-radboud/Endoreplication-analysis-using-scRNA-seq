library(xml2)
library(RCurl)
library(ggplot2)
library(Seurat)
library(dplyr)
library(tidyr)
library(ggpubr)
library(devtools)
if (!requireNamespace("ggPlantmap", quietly = TRUE)) {
  devtools::install_github("leonardojo/ggPlantmap")
} else {
  message("ggPlantmap is already installed.")
}
library(ggPlantmap)
library(gridExtra)
library(patchwork)
library(readr)
library(readxl)
library(cowplot)
library(stringr)
library(RColorBrewer)
library(openxlsx)

# setwd()
setwd("C:/Users/Олеся/Desktop/graduation/data_sets/realistic_layout/")

# Load object
AT1 <- readRDS("AT1.rds")
AT2 <- readRDS("AT2.rds")
T1 <- readRDS("T1.rds")
T2 <- readRDS("T2.rds")



# Function to calculate a point on a cubic Bézier curve
calculate_bezier_point <- function(t, P0, P1, P2, P3) {
  (1 - t)^3 * P0 + 3 * (1 - t)^2 * t * P1 + 3 * (1 - t) * t^2 * P2 + t^3 * P3
}

# Creating new ggPlantmap data frame
new.ggPlantmap = data.frame(matrix( 
  vector(), 0, 7, dimnames=list(c(), c("ROI.name", "Cell_type", "Colour", "ROI.id","point","x","y"))), 
  stringsAsFactors=F)

# load in svg_file
svg_file <- readLines("Epidermis_Olesia.svg", warn = FALSE)
svg_content <- paste(svg_file, collapse = "\n")

# Parse the SVG content
svg <- read_xml(svg_content)

# Find all path elements in the SVG
ns <- c(svg = "http://www.w3.org/2000/svg")
paths <- xml_find_all(svg, "//svg:path", ns)
paths <- paths[grepl("[zZ]", xml_attr(paths, "d"))]

# give ID
for (i in seq_along(paths)) {
  xml_set_attr(paths[[i]], "id", paste0("poly", i))
}

# take styles
style_node <- xml_find_first(svg, "//svg:style", ns)

# Extract the CSS text
css_text <- xml_text(style_node)

# Remove newlines and tabs for easier parsing
css_text <- gsub("[\r\n\t]", "", css_text)

# Split by class definitions
# This assumes CSS like ".st0{fill:#FFFFFF;stroke:#010101;}.st1{...}"
rules <- str_split(css_text, "\\.st")[[1]]
rules <- rules[nzchar(rules)]

# Parse into a named list: names = st0, st1, ...
css_list <- lapply(rules, function(x) {
  parts <- str_match(x, "(\\d+)\\{([^}]*)\\}")
  if (!is.na(parts[1])) {
    class_num <- paste0("st", parts[2])
    style_val <- parts[3]
    return(setNames(style_val, class_num))
  }
  return(NULL)
})

# Combine list into named vector
css_list <- unlist(css_list)

for (path_node in paths) {
  # Get class attribute (e.g., "st0")
  class_val <- xml_attr(path_node, "class")
  
  # Find style string for that class
  style_val <- css_list[class_val]
  
  # If style found, add style attribute
  if (!is.na(style_val)) {
    xml_set_attr(path_node, "style", style_val)
  }
}

# Loop through each path and extract information
for (path in paths) {
  
  fill_value <- xml_attr(path, "style")
  fill_value <- sub(".*fill:([^;]*);.*", "\\1", fill_value)
  
  if (fill_value!='#231f20' && !is.na(fill_value)) {
    
    path_id<-xml_attr(path, "id")
    
    # Extract and add each point
    path_d <- xml_attr(path, "d")
    coordinates <- str_extract_all(path_d, "[A-Za-z]|-?\\d*\\.?\\d+|-?\\d+\\.?\\d*")[[1]]
    
    x<-0
    y<-0
    l<-1
    k<-1
    point<-0
    for (i in seq(1, length(coordinates), by = 1)) {
      # Update coordinates using relative values
      
      if ((coordinates[i])=='h') {
        h<-1
        l<-0
        v<-0
        c<-0
        k<-1
      } else if ((coordinates[i])=='v') {
        h<-0
        l<-0
        v<-1
        c<-0
        k<-1
      } else if ((coordinates[i])=='m' | (coordinates[i])=='l') {
        h<-0
        l<-1
        v<-0
        c<-0
        k<-1
      } else if ((coordinates[i])=='M' | (coordinates[i])=='L') {
        h<-0
        l<-2
        v<-0
        c<-0
        k<-1
      } else if ((coordinates[i])=='H') {
        h<-2
        l<-0
        v<-0
        c<-0
        k<-1
      } else if ((coordinates[i])=='V') {
        h<-0
        l<-0
        v<-2
        c<-0
        k<-1
      } else if ((coordinates[i])=='C') {
        h<-0
        l<-0
        v<-0
        c<-2
        k<-1
      } else if ((coordinates[i])=='c') {
        h<-0
        l<-0
        v<-0
        c<-1
        k<-1
      } else if ((coordinates[i])!='l'&(coordinates[i])!='h'&(coordinates[i])!='v' & l==1 & k==1) {
        x <- x + as.numeric(coordinates[i])
        y <- y + as.numeric(coordinates[i + 1])
        point<-point+1
        new.ggPlantmap[nrow(new.ggPlantmap)+1,]<-c(path_id,NA,fill_value,path_id,point,x,y)
        k<-2
      } else if ((coordinates[i])!='l'&(coordinates[i])!='h'&(coordinates[i])!='v' & l==2 & k==1) {
        x <- as.numeric(coordinates[i])
        y <- as.numeric(coordinates[i + 1])
        point<-point+1
        new.ggPlantmap[nrow(new.ggPlantmap)+1,]<-c(path_id,NA,fill_value,path_id,point,x,y)
        k<-2
      } else if (k>1) {
        k<-k-1
      } else if ((coordinates[i])!='l'&(coordinates[i])!='h'&(coordinates[i])!='v' & h==1) {
        x <- x + as.numeric(coordinates[i])
        point<-point+1
        new.ggPlantmap[nrow(new.ggPlantmap)+1,]<-c(path_id,NA,fill_value,path_id,point,x,y)
      } else if ((coordinates[i])!='l'&(coordinates[i])!='h'&(coordinates[i])!='v' & h==2) {
        x <- as.numeric(coordinates[i])
        point<-point+1
        new.ggPlantmap[nrow(new.ggPlantmap)+1,]<-c(path_id,NA,fill_value,path_id,point,x,y)
      } else if ((coordinates[i])!='l'&(coordinates[i])!='h'&(coordinates[i])!='v' & v==1) {
        y <- y + as.numeric(coordinates[i])
        point<-point+1
        new.ggPlantmap[nrow(new.ggPlantmap)+1,]<-c(path_id,NA,fill_value,path_id,point,x,y)
      } else if ((coordinates[i])!='l'&(coordinates[i])!='h'&(coordinates[i])!='v' & v==2) {
        y <- as.numeric(coordinates[i])
        point<-point+1
        new.ggPlantmap[nrow(new.ggPlantmap)+1,]<-c(path_id,NA,fill_value,path_id,point,x,y)
      } else if ((coordinates[i])!='l'&(coordinates[i])!='h'&(coordinates[i])!='v' & c==2 & k==1) {
        # Absolute cubic Bézier curve
        P0_x <- x
        P0_y <- y
        P1_x <- as.numeric(coordinates[i]) # First control point
        P1_y <- as.numeric(coordinates[i + 1])
        P2_x <- as.numeric(coordinates[i + 2]) # Second control point
        P2_y <- as.numeric(coordinates[i + 3])
        P3_x <- as.numeric(coordinates[i + 4]) # Endpoint
        P3_y <- as.numeric(coordinates[i + 5])
        
        # Sample points along the curve
        num_samples <- 10  # Number of points to approximate the curve
        for (t in seq(0, 1, length.out = num_samples)) {
          sampled_x <- calculate_bezier_point(t, P0_x, P1_x, P2_x, P3_x)
          sampled_y <- calculate_bezier_point(t, P0_y, P1_y, P2_y, P3_y)
          point <- point + 1
          new.ggPlantmap[nrow(new.ggPlantmap) + 1, ] <- c(path_id,NA,fill_value,path_id,point,sampled_x,sampled_y)
        }
        
        # Update the current position to the endpoint of the curve
        x <- P3_x
        y <- P3_y
        
        k<-6
      }
      
      else if ((coordinates[i])!='l'&(coordinates[i])!='h'&(coordinates[i])!='v' & c==1 & k==1) {
        # Relative cubic Bézier curve
        P0_x <- x
        P0_y <- y
        P1_x <- x + as.numeric(coordinates[i])  # First control point
        P1_y <- y + as.numeric(coordinates[i + 1])
        P2_x <- x + as.numeric(coordinates[i + 2])  # Second control point
        P2_y <- y + as.numeric(coordinates[i + 3])
        P3_x <- x + as.numeric(coordinates[i + 4])  # Endpoint
        P3_y <- y + as.numeric(coordinates[i + 5])
        
        # Sample points along the curve
        num_samples <- 10  # Number of points to approximate the curve
        for (t in seq(0, 1, length.out = num_samples)) {
          sampled_x <- calculate_bezier_point(t, P0_x, P1_x, P2_x, P3_x)
          sampled_y <- calculate_bezier_point(t, P0_y, P1_y, P2_y, P3_y)
          point <- point + 1
          new.ggPlantmap[nrow(new.ggPlantmap) + 1, ] <- c(path_id,NA,fill_value,path_id,point,sampled_x,sampled_y)
        }
        
        # Update the current position to the endpoint of the curve
        x <- P3_x
        y <- P3_y
        
        k<-6
      }
      
    }
  }
}

# as character
{
  new.ggPlantmap$ROI.name = as.character(new.ggPlantmap$ROI.name) #character()
  new.ggPlantmap$Cell_type = as.character(new.ggPlantmap$Cell_type) #character()
  new.ggPlantmap$Colour = as.character(new.ggPlantmap$Colour) #character()
  new.ggPlantmap$ROI.id = as.character(new.ggPlantmap$ROI.id) #integer()
  new.ggPlantmap$point = as.numeric(new.ggPlantmap$point) #integer()
  new.ggPlantmap$x = as.numeric(new.ggPlantmap$x) #double()
  new.ggPlantmap$y = as.numeric(new.ggPlantmap$y) #double()
}

# ??
new.ggPlantmap$y<-new.ggPlantmap$y+2*(mean(max(new.ggPlantmap$y, na.rm = TRUE),min(new.ggPlantmap$y, na.rm = TRUE))-new.ggPlantmap$y)

# draw plot
ggPlantmap.plot(data = new.ggPlantmap, layer = Colour)  


# extract colour mapping
unique(new.ggPlantmap$Colour)

# assign cell types to colours
new.ggPlantmap$Cell_type[new.ggPlantmap$Colour=="#FCD3C1" ]<-"Atrichoblast1-m1"
new.ggPlantmap$Cell_type[new.ggPlantmap$Colour=="#FFF6D4" ]<-"Atrichoblast1-m2"
new.ggPlantmap$Cell_type[new.ggPlantmap$Colour=="#FFEFAE"]<-"Atrichoblast1-m2"
new.ggPlantmap$Cell_type[new.ggPlantmap$Colour=="#FFE35B"]<-"Atrichoblast1-t"
new.ggPlantmap$Cell_type[new.ggPlantmap$Colour=="#FFDE17" ]<-"Atrichoblast1-e1"
new.ggPlantmap$Cell_type[new.ggPlantmap$Colour=="#CAB453"]<-"Atrichoblast1-e2"
new.ggPlantmap$Cell_type[new.ggPlantmap$Colour=="#958A5B"]<-"Atrichoblast1-d"

new.ggPlantmap$Cell_type[new.ggPlantmap$Colour=="#F9AA8F"]<-"Atrichoblast2-m1"
new.ggPlantmap$Cell_type[new.ggPlantmap$Colour=="#FFFFFF"]<-"Atrichoblast2-m2"
new.ggPlantmap$Cell_type[new.ggPlantmap$Colour=="#DBD3D2"]<-"Atrichoblast2-m2"
new.ggPlantmap$Cell_type[new.ggPlantmap$Colour=="#EBF5EA"]<-"Atrichoblast2-t"
new.ggPlantmap$Cell_type[new.ggPlantmap$Colour=="#D1E9D1"]<-"Atrichoblast2-e1"
new.ggPlantmap$Cell_type[new.ggPlantmap$Colour=="#ADBCAF"]<-"Atrichoblast2-e2"
new.ggPlantmap$Cell_type[new.ggPlantmap$Colour=="#868F89"]<-"Atrichoblast2-d"

new.ggPlantmap$Cell_type[new.ggPlantmap$Colour=="#D3E5D1"]<-"Trichoblast1-m1"
new.ggPlantmap$Cell_type[new.ggPlantmap$Colour=="#ACD0AA"]<-"Trichoblast1-m2"
new.ggPlantmap$Cell_type[new.ggPlantmap$Colour=="#84BE88"]<-"Trichoblast1-m2"
new.ggPlantmap$Cell_type[new.ggPlantmap$Colour=="#59AE66"]<-"Trichoblast1-t"
new.ggPlantmap$Cell_type[new.ggPlantmap$Colour=="#00A14B"]<-"Trichoblast1-e1"
new.ggPlantmap$Cell_type[new.ggPlantmap$Colour=="#4D8E57"]<-"Trichoblast1-e2"
new.ggPlantmap$Cell_type[new.ggPlantmap$Colour=="#56765A"]<-"Trichoblast1-d"

new.ggPlantmap$Cell_type[new.ggPlantmap$Colour=="#798B91"]<-"Trichoblast2-m1"
new.ggPlantmap$Cell_type[new.ggPlantmap$Colour=="#CCCADB"]<-"Trichoblast2-m2"
new.ggPlantmap$Cell_type[new.ggPlantmap$Colour=="#A2A1BF"]<-"Trichoblast2-m2"
new.ggPlantmap$Cell_type[new.ggPlantmap$Colour=="#908FB3"]<-"Trichoblast2-t"
new.ggPlantmap$Cell_type[new.ggPlantmap$Colour=="#85859C"]<-"Trichoblast2-e1"
new.ggPlantmap$Cell_type[new.ggPlantmap$Colour=="#72727F"]<-"Trichoblast2-e2"
new.ggPlantmap$Cell_type[new.ggPlantmap$Colour=="#56565C"]<-"Trichoblast2-d"


# save new.ggplantmap
write.csv(new.ggPlantmap, "new_ggPlantmap_working.csv")

# set default assay
DefaultAssay(AT1) <- "RNA"
DefaultAssay(AT2) <- "RNA"
DefaultAssay(T1) <- "RNA"
DefaultAssay(T2) <- "RNA"

#set idents
Idents(AT1) <- "Cell_type"
Idents(AT2) <- "Cell_type"
Idents(T1) <- "Cell_type"
Idents(T2) <- "Cell_type"


# Average expression AT1
avg_exp <- AverageExpression(AT1, assays = "RNA")
colnames(avg_exp$RNA)
avg_exp<-avg_exp$RNA
colnames(avg_exp)
avg_exp1<-t(avg_exp)
avg_exp1<-as.data.frame((avg_exp1))
row.names(avg_exp1)
avg_exp1<-cbind.data.frame(row.names(avg_exp1), avg_exp1)
colnames(avg_exp1)[1] <- "Cell_type"

# Average expression AT2
avg_exp <- AverageExpression(AT2, assays = "RNA")
colnames(avg_exp$RNA)
avg_exp<-avg_exp$RNA
colnames(avg_exp)
avg_exp2<-t(avg_exp)
avg_exp2<-as.data.frame((avg_exp2))
row.names(avg_exp2)
avg_exp2<-cbind.data.frame(row.names(avg_exp2), avg_exp2)
colnames(avg_exp2)[1] <- "Cell_type"

# Average expression T1
avg_exp <- AverageExpression(T1, assays = "RNA")
colnames(avg_exp$RNA)
avg_exp<-avg_exp$RNA
colnames(avg_exp)
avg_exp3<-t(avg_exp)
avg_exp3<-as.data.frame((avg_exp3))
row.names(avg_exp3)
avg_exp3<-cbind.data.frame(row.names(avg_exp3), avg_exp3)
colnames(avg_exp3)[1] <- "Cell_type"

#Average expression T2
avg_exp <- AverageExpression(T2, assays = "RNA")
colnames(avg_exp$RNA)
avg_exp<-avg_exp$RNA
colnames(avg_exp)
avg_exp4<-t(avg_exp)
avg_exp4<-as.data.frame((avg_exp4))
row.names(avg_exp4)
avg_exp4<-cbind.data.frame(row.names(avg_exp4), avg_exp4)
colnames(avg_exp4)[1] <- "Cell_type"

#Combining DFs together

all_genes <- unique(unlist(lapply(list(avg_exp1, avg_exp2, avg_exp3, avg_exp4), colnames)))

pad_genes <- function(df, all_genes) {
  missing_genes <- setdiff(all_genes, colnames(df))
  df[missing_genes] <- NA  # Add missing genes as NA
  return(df[all_genes])    # Reorder columns to match all_genes
}

df1_p <- pad_genes(avg_exp1, all_genes)
df2_p <- pad_genes(avg_exp2, all_genes)
df3_p <- pad_genes(avg_exp3, all_genes)
df4_p <- pad_genes(avg_exp4, all_genes)

merged_avg <- rbind(df1_p, df2_p, df3_p, df4_p)

# Merge

new.ggPlantmap_atlas <- ggPlantmap.merge(new.ggPlantmap, merged_avg, "Cell_type", "Cell_type")
saveRDS(new.ggPlantmap_atlas, file = "C:/Users/Олеся/Desktop/graduation/data_sets/realistic_layout/epidermis_map.rds")


new.ggPlantmap_atlas<- readRDS("C:/Users/Олеся/Desktop/graduation/data_sets/realistic_layout/epidermis_map.rds")
#some plots
markers<- read.xlsx("C:/Users/Олеся/Desktop/graduation/data_sets/Cell_cycle_related.xlsx")
gene_labels <- setNames(as.character(markers$gene_labels), markers$TAIR_ID)
genes <- markers$TAIR_ID

# change function
ggPlantmap.heatmap_edit <- function(map.quant, value.quant = value.quant, show.legend = TRUE, linewidth = 0.5) {
  ggPlantmap <- map.quant
  ar <- (max(ggPlantmap$y) - min(ggPlantmap$y)) / (max(ggPlantmap$x) - min(ggPlantmap$x))
  
  ggplot2::ggplot(ggPlantmap, aes(x = x, y = y)) +
    geom_polygon(
      aes(group = ROI.id, fill = .data[[value.quant]]),
      colour = "black",
      linewidth = linewidth,
      show.legend = show.legend
    ) +
    theme_void() +
    theme(
      panel.grid = element_blank(),
      legend.position = "right",
      aspect.ratio = ar
    )
}

setwd("C:/Users/Олеся/Desktop/graduation/data_sets/realistic_layout/")
# plot loop
for (gene_to_plot in genes) {
  # Skip if gene not in dataset
  if (!gene_to_plot %in% colnames(new.ggPlantmap_atlas)) {
    message("Skipping ", gene_to_plot, ": not found in expression dataset.")
    next
  }
  
  clean_atlas <- new.ggPlantmap_atlas %>% 
    dplyr::filter(is.finite(x), is.finite(y))
  
  # Lookup the label, fallback to ID if not found
  label <- if (!is.na(gene_labels[gene_to_plot]) && gene_labels[gene_to_plot] != "") {
    gene_labels[gene_to_plot]
  } else {
    gene_to_plot
  }
  print(label)
  p <- ggPlantmap.heatmap_edit(clean_atlas, gene_to_plot) +
    scale_fill_distiller(palette = "YlOrRd", direction = 1, na.value = "white") +
    ggtitle(label)
  
  ggsave(filename = paste0("plot/Epidermis_", label, "_heatmap.png"), plot = p, width = 2.7, height = 3.6, dpi = 300, bg = "white")
}
