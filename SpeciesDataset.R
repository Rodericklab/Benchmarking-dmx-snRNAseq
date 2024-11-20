library(BUSpaRse)
library(tidyverse)
library(DropletUtils)
library(Matrix)
library(Seurat)
library(ggplot2)
library(liger)
library(Matrix)
library(rliger)
library(Seurat)
library(SeuratWrappers)
library(hdf5r)
library(rhdf5)
library(readr)
library(DoubletFinder)
library(Azimuth)
library(SeuratData)
library(dplyr)
library(magrittr)
library(patchwork)
library(viridis)
library(scCustomize)
library(qs)
library(ggpubr)
library(wesanderson)
library(paletteer)

#### 01 read in data and demultiplx human and sheep ####
cell_bender_mat <- Read_CellBender_h5_Mat(file_name = "data/human_sheep/cellbender/map_kallisto/raw_feature_bc_hush_cellbender.h5")
dim(cell_bender_mat)

## remove empty droplets
tot_counts <- Matrix::colSums(cell_bender_mat)
summary(tot_counts)

## Compute barcode rank from Dropletutils
bc_rank <- barcodeRanks(cell_bender_mat)

qplot(bc_rank$total, bc_rank$rank, geom = "line") +
  geom_vline(xintercept = metadata(bc_rank)$knee, color = "darkblue", linetype = 2) +
  geom_vline(xintercept = metadata(bc_rank)$inflection, color = "darkred", linetype = 2) +
  annotate("text", y = 1000, x = 1.5 * c(metadata(bc_rank)$knee, metadata(bc_rank)$inflection),
           label = c("knee", "inflection"), color = c("darkblue", "darkred")) +
  scale_x_log10() +
  scale_y_log10() +
  labs(y = "Barcode rank", x = "Total UMI count") +
  # theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = "black", fill=NA, size=1.5))

## filter the cells
cell_bender_mat <- cell_bender_mat[, tot_counts > metadata(bc_rank)$inflection]
dim(cell_bender_mat)

gene_species <- ifelse(str_detect(rownames(cell_bender_mat), "^ENSOARG"), "sheep", "human")
sheep_inds <- gene_species == "sheep"
human_inds <- gene_species == "human"

# mark cells as sheep or human
cell_species <- tibble(n_sheep_umi = Matrix::colSums(cell_bender_mat[sheep_inds,]),
                       n_human_umi = Matrix::colSums(cell_bender_mat[human_inds,]),
                       tot_umi = Matrix::colSums(cell_bender_mat),
                       prop_sheep = n_sheep_umi / tot_umi,
                       prop_human = n_human_umi / tot_umi)

cell_species<- cell_species %>% 
  mutate(Species = case_when(
    prop_sheep > 0.7 ~ "Sheep",
    prop_human > 0.7 ~ "Human",
    TRUE ~ "Doublets"
  ))
table(cell_species$Species)

ggscatter(cell_species, x = "n_human_umi", y = "n_sheep_umi",
          palette = c("#777772", "darkred", "darkblue"), shape = 20, size = 1, 
          color = "Species"
)

cell_species$barcode <- paste0(colnames(cell_bender_mat), "-1")

#### 02 process data mapped on human ####
map_human <- Read10X("data/human_sheep/cellranger/map_human/filtered_feature_bc_matrix/")

seurat_human <- CreateSeuratObject(map_human, project = "Map_human", min.cells = 3, min.features = 200)
seurat_human[["percent.mt"]] <- PercentageFeatureSet(seurat_human, pattern = "^MT-")
seurat_human <- subset(seurat_human, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Visualize QC metrics as a violin plot
VlnPlot(seurat_human, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

## add dmx info for human
seurat_human@meta.data$barcode <- rownames(seurat_human@meta.data)
meta_human <- seurat_human@meta.data

kalisto_dmx <- read_csv("results/hush/map_kallisto/kalisto_dmx.csv")
meta_human_dmx <- left_join(meta_human, kalisto_dmx, by = "barcode")
sum(is.na(meta_human_dmx$Species))
seurat_human <- AddMetaData(seurat_human, meta_human_dmx$Species, col.name = "Species")

## remove unknown cells
Idents(seurat_human) <- "Species"
seurat_human_subset <- subset(x = seurat_human, idents = c ("Human", "Doublets"))
seurat_human_subset$dmx_doublet <- ifelse(seurat_human_subset$Species == "Human", "Singlets", "Doublets")

## normalize data
seurat_human_subset <- seurat_human_subset %>% 
  NormalizeData(verbose = FALSE) %>% 
  ScaleData(verbose = FALSE) %>% 
  FindVariableFeatures(verbose = FALSE) %>%
  RunPCA(verbose = FALSE, npcs = 30) %>%
  RunTSNE(dims = 1:20, check_duplicates = FALSE) %>%
  RunUMAP(dims = 1:20, check_duplicates = FALSE) %>%
  FindNeighbors(dims = 1:10) %>%
  FindClusters(resolution = 0.5)

## dmx doublets
DimPlot(seurat_human_subset, reduction = "umap", pt.size = 0.5, group.by = "dmx_doublet", cols = c("red", "black"))

## identify the heterotypic doubelts
## pK Identification (no ground-truth) 
sweep.res.list_human <- paramSweep(seurat_human_subset, PCs = 1:10, sct = FALSE)
sweep.stats_human <- summarizeSweep(sweep.res.list_human, GT = FALSE)
bcmvn_human <- find.pK(sweep.stats_human)

## Homotypic Doublet Proportion Estimate
homotypic.prop <- modelHomotypic(seurat_human_subset@meta.data$seurat_clusters)           
nExp_poi <- round(0.075*nrow(seurat_human_subset@meta.data))  
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies
seurat_human_subset <- doubletFinder(seurat_human_subset, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
DimPlot(seurat_human_subset, reduction = "umap", pt.size = 0.5, group.by = "DF.classifications_0.25_0.09_223", cols = c("red", "black"))

## all doublets
seurat_human_subset@meta.data$All_doublets <- ifelse(seurat_human_subset@meta.data$dmx_doublet == "Doublets" | 
                                                   seurat_human_subset@meta.data$DF.classifications_0.25_0.09_223 == "Doublet", "Doublets", "Singlets")
DimPlot(seurat_human_subset, reduction = "umap", pt.size = 0.5, group.by = "All_doublets", cols = c("red", "black"))

## remove doublets
Idents(seurat_human_subset) <- "All_doublets"
seurat_human_subset_clean <- subset(x = seurat_human_subset, idents = c("Singlets"))

seurat_human_subset_clean <- seurat_human_subset_clean %>% 
  NormalizeData(verbose = FALSE) %>% 
  ScaleData(verbose = FALSE) %>% 
  FindVariableFeatures(verbose = FALSE) %>%
  RunPCA(verbose = FALSE, npcs = 30) %>%
  RunTSNE(dims = 1:20, check_duplicates = FALSE) %>%
  RunUMAP(dims = 1:20, check_duplicates = FALSE) %>%
  FindNeighbors(dims = 1:10) %>%
  FindClusters(resolution = 0.5)

DimPlot(seurat_human_subset_clean, reduction = "umap", pt.size = 0.5, group.by = "seurat_clusters")

## annotate cells 
seurat_human_subset_clean <- RunAzimuth(seurat_human_subset_clean, reference = "data/reference/heart/")
color_scheme <- readRDS("results/application/GC151435/processed/color_scheme.rds")
DimPlot(seurat_human_subset_clean, reduction = "umap", group.by = "predicted.celltype.l1", raster =FALSE, cols = color_scheme)

### cell composition
plt_cell_count_human <- seurat_human_subset_clean@meta.data[,c("Species", "predicted.celltype.l1")]
colnames(plt_cell_count_human)<-c("Sample", "Celltype")

ggplot(plt_cell_count_human,aes(x=Sample,fill=Celltype))+
  geom_bar(position="fill")+
  scale_fill_manual(values = color_scheme) +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())


#### 03 process data mapped on sheep ####
map_sheep <- Read10X("data/human_sheep/cellranger/map_sheep/filtered_feature_bc_matrix/")

seurat_sheep <- CreateSeuratObject(map_sheep, project = "Map_sheep", min.cells = 3, min.features = 200)
seurat_sheep[["percent.mt"]] <- PercentageFeatureSet(seurat_sheep, pattern = "^MT-")
seurat_sheep <- subset(seurat_sheep, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)

# Visualize QC metrics as a violin plot
VlnPlot(seurat_sheep, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

## add dmx info
seurat_sheep@meta.data$barcode <- rownames(seurat_sheep@meta.data)
meta_sheep <- seurat_sheep@meta.data

kalisto_dmx <- read_csv("results/hush/map_kallisto/kalisto_dmx.csv")

meta_sheep_dmx <- left_join(meta_sheep, kalisto_dmx, by = "barcode")
sum(is.na(meta_sheep_dmx$Species))

seurat_sheep <- AddMetaData(seurat_sheep, meta_sheep_dmx$Species, col.name = "Species")

Idents(seurat_sheep) <- "Species"
seurat_sheep_subset <- subset(x = seurat_sheep, idents = c ("Sheep", "Doublets"))
seurat_sheep_subset$dmx_doublet <- ifelse(seurat_sheep_subset$Species == "Sheep", "Singlets", "Doublets")

## normalize data
seurat_sheep_subset <- seurat_sheep_subset %>% 
  NormalizeData(verbose = FALSE) %>% 
  ScaleData(verbose = FALSE) %>% 
  FindVariableFeatures(verbose = FALSE) %>%
  RunPCA(verbose = FALSE, npcs = 30) %>%
  RunTSNE(dims = 1:20, check_duplicates = FALSE) %>%
  RunUMAP(dims = 1:20, check_duplicates = FALSE) %>%
  FindNeighbors(dims = 1:10) %>%
  FindClusters(resolution = 0.5)

DimPlot(seurat_sheep_subset, reduction = "umap", pt.size = 0.5, group.by = "dmx_doublet", cols = c("red", "black"))

## identify heterotypic doubelts

## pK Identification (no ground-truth) 
sweep.res.list_sheep <- paramSweep(seurat_sheep_subset, PCs = 1:10, sct = FALSE)
sweep.stats_sheep <- summarizeSweep(sweep.res.list_sheep, GT = FALSE)
bcmvn_sheep <- find.pK(sweep.stats_sheep)

## Homotypic Doublet Proportion Estimate 
homotypic.prop <- modelHomotypic(seurat_sheep_subset@meta.data$seurat_clusters)          
nExp_poi <- round(0.075*nrow(seurat_sheep_subset@meta.data))  
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies 
seurat_sheep_subset <- doubletFinder(seurat_sheep_subset, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
DimPlot(seurat_sheep_subset, reduction = "umap", pt.size = 0.5, group.by = "DF.classifications_0.25_0.09_687", cols = c("red", "black"))

## combine all doublets
seurat_sheep_subset@meta.data$All_doublets <- ifelse(seurat_sheep_subset@meta.data$dmx_doublet == "Doublets" | 
                                                       seurat_sheep_subset@meta.data$DF.classifications_0.25_0.09_687 == "Doublet", "Doublets", "Singlets")
DimPlot(seurat_sheep_subset, reduction = "umap", pt.size = 0.5, group.by = "All_doublets", cols = c("red", "black"))

## subset to singlets
Idents(seurat_sheep_subset) <- "All_doublets"
seurat_sheep_subset_clean <- subset(x = seurat_sheep_subset, idents = c("Singlets"))

seurat_sheep_subset_clean <- seurat_sheep_subset_clean %>% 
  NormalizeData(verbose = FALSE) %>% 
  ScaleData(verbose = FALSE) %>% 
  FindVariableFeatures(verbose = FALSE) %>%
  RunPCA(verbose = FALSE, npcs = 30) %>%
  RunTSNE(dims = 1:20, check_duplicates = FALSE) %>%
  RunUMAP(dims = 1:20, check_duplicates = FALSE) %>%
  FindNeighbors(dims = 1:10) %>%
  FindClusters(resolution = 0.5)

## annotate cells
color_scheme <- readRDS("results/application/GC151435/processed/color_scheme.rds")
seurat_sheep_subset_clean <- RunAzimuth(seurat_sheep_subset_clean, reference = "data/reference/heart/")
DimPlot(seurat_sheep_subset_clean, reduction = "ref.umap", pt.size = 0.5, group.by = "seurat_clusters")
DimPlot(seurat_sheep_subset_clean, reduction = "ref.umap", group.by = "predicted.celltype.l1", raster =FALSE, cols = color_scheme)

## cell composition
plt_cell_count_sheep <- seurat_sheep_subset_clean@meta.data[,c("Species", "predicted.celltype.l1")]
colnames(plt_cell_count_sheep)<-c("Sample", "Celltype")

ggplot(plt_cell_count_sheep,aes(x=Sample,fill=Celltype))+
  geom_bar(position="fill")+
  scale_fill_manual(values = color_scheme) +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())

#### 04 plot doublet number ####
table(seurat_sheep_subset@meta.data$dmx_doublet, seurat_sheep_subset@meta.data$DF.classifications_0.25_0.09_687)
table(seurat_human_subset@meta.data$dmx_doublet, seurat_human_subset@meta.data$DF.classifications_0.25_0.09_223)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
color <- gg_color_hue(8); color

data <- data.frame(
  Group = c("Sheep", "Sheep", "Sheep", "Sheep", 
            "Human", "Human","Human","Human"),
  Type = c("Singlets", "Doublets_Dmultiplexing", "Doublets_DoubletFinder", "Doublets_Both", 
               "Singlets", "Doublets_Dmultiplexing", "Doublets_DoubletFinder", "Doublets_Both"),
  Value = c(8085, 385, 649, 38, 
            2361, 386, 186, 37)
)
data <- data %>%
  group_by(Group) %>%
  mutate(Total = sum(Value),        
         Percentage = Value / Total * 100)  
ggplot(data, aes(x = Group, y = Percentage, fill = Type)) + 
  geom_bar(stat = "identity") +
  labs(title = "Doublet Porportion", x = "", y = "Percentage (%)") + 
  scale_fill_manual(values = color) +
  theme(axis.title.y = element_text(size = 16),  # Make the y-axis title bigger and bold
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        plot.title = element_text(size = 16, hjust = 0.5))  +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.background = element_rect(fill = "white", colour = "white")) 
