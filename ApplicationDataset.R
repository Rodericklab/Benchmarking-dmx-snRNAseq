library(dplyr)
library(Seurat)
library(ggplot2)
library(liger)
library(rliger)
library(Matrix)
library(rliger)
library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(hdf5r)
library(rhdf5)
library(dplyr)
library(readr)
library(DoubletFinder)
library(Seurat)
library(Azimuth)
library(SeuratData)
library(patchwork)
library(scCustomize)
library(wesanderson)
library(paletteer)

#### 01 read in the data #### 
NucSeq.data <- Read10X_h5("data/application/GC151435/filtered_feature_bc_matrix.h5")
NucSeq <- CreateSeuratObject(
  counts = NucSeq.data
)

#### 02 read in the vireo results ####
dmx_vireo <- read_tsv("data/application/vireo/donor_ids.tsv") %>% dplyr::select(cell, donor_id, best_singlet)

NucSeq$donor_id <- dmx_vireo$donor_id
NucSeq$best_singlet <- dmx_vireo$best_singlet

NucSeq$sample <- ifelse(NucSeq$best_singlet == "donor0", "AF65",
                        ifelse(NucSeq$best_singlet == "donor1", "AF25",
                               ifelse(NucSeq$best_singlet == "donor2", "AF10",
                                      ifelse(NucSeq$best_singlet == "donor3", "AF64",
                                                    NucSeq$best_singlet))))
NucSeq$dmx_singlets <- ifelse(NucSeq$donor_id %in% c("donor0" ,"donor1", "donor2", "donor3"), TRUE, FALSE)

#### 03 add QC info ####
NucSeq[["percent.mt"]] <- PercentageFeatureSet(NucSeq, pattern = "^MT-")
NucSeq <- subset(NucSeq, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Visualize QC metrics as a violin plot
VlnPlot(NucSeq, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#### 04 identify doublets ####
## Pre-process Seurat object (standard) --------------------------------------------------------------------------------------
NucSeq <- NormalizeData(NucSeq)
NucSeq <- FindVariableFeatures(NucSeq, selection.method = "vst", nfeatures = 2000)
NucSeq <- ScaleData(NucSeq)
NucSeq <- RunPCA(NucSeq)
NucSeq <- FindNeighbors(NucSeq, dims = 1:10)
NucSeq <- FindClusters(NucSeq, resolution = 0.5)
NucSeq <- RunUMAP(NucSeq, dims = 1:10)

## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.list_Nuc <- paramSweep(NucSeq, PCs = 1:10, sct = FALSE)
sweep.stats_Nuc <- summarizeSweep(sweep.res.list_Nuc, GT = FALSE)
bcmvn_Nuc <- find.pK(sweep.stats_Nuc)

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
homotypic.prop <- modelHomotypic(NucSeq@meta.data$seurat_clusters)           
nExp_poi <- round(0.075*nrow(NucSeq@meta.data))  
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
NucSeq <- doubletFinder(NucSeq, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

DimPlot(NucSeq, reduction = "umap", group.by = "dmx_singlets", cols = c("red", "black"), raster =FALSE)
DimPlot(NucSeq, reduction = "umap", group.by = "DF.classifications_0.25_0.09_813", cols = c("red", "black"), raster =FALSE)
DimPlot(NucSeq, reduction = "umap", group.by = "sample", raster =FALSE)
DimPlot(NucSeq, reduction = "umap", group.by = "seurat_clusters", raster =FALSE)

#### 05 merge all of doublets ####
NucSeq@meta.data$All_doublets <- ifelse(NucSeq@meta.data$dmx_singlets == "FALSE" | 
                                          NucSeq@meta.data$DF.classifications_0.25_0.09_813 == "Doublet", "Doublets", "Singlets")
NucSeq@meta.data$both <- ifelse(NucSeq@meta.data$dmx_singlets == "FALSE" & 
                                          NucSeq@meta.data$DF.classifications_0.25_0.09_813 == "Doublet", "Doublets", "Singlets")
DimPlot(NucSeq, reduction = "umap", pt.size = 0.5, group.by = "All_doublets", cols = c("red", "black"))

table(NucSeq@meta.data$sample, NucSeq@meta.data$dmx_singlets)
table(NucSeq@meta.data$sample, NucSeq@meta.data$DF.classifications_0.25_0.09_813)
table(NucSeq@meta.data$sample, NucSeq@meta.data$All_doublets)
table(NucSeq@meta.data$sample, NucSeq@meta.data$both)

#### 06 annotating the cells with doublets ####  
DimPlot(NucSeq_dblts, reduction = "umap", group.by = "sample", raster =FALSE)
DimPlot(NucSeq_dblts, reduction = "umap", group.by = "seurat_clusters", raster =FALSE)

NucSeq_dblts <- RunAzimuth(NucSeq_dblts, reference = "data/reference/heart/")

color_scheme <- readRDS("results/application/GC151435/processed/color_scheme.rds")
DimPlot(NucSeq_dblts, reduction = "umap", group.by = "predicted.celltype.l1", raster =FALSE, cols = color_scheme)


#### 07 calculate cell composition with doublets ####
plt_cell_count_doublets<-NucSeq_dblts@meta.data[,c("sample", "predicted.celltype.l1")]
colnames(plt_cell_count_doublets)<-c("Sample", "Celltype")

ggplot(plt_cell_count_doublets,aes(x=Sample,fill=Celltype))+
  geom_bar(position="fill")+
  scale_fill_manual(values = color_scheme) +
  scale_x_discrete(labels = c("AF10" = "Donor 1", 
                              "AF25" = "Donor 2", 
                              "AF64" = "Donor 3", 
                              "AF65" = "Donor 4")) +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())

#### 08 removing the doublets, subset, normalizing and annotating #### 
NucSeq
Idents(NucSeq) <- "All_doublets"
NucSeq <- subset(x = NucSeq, idents = "Singlets")

NucSeq <- NormalizeData(NucSeq)
NucSeq <- FindVariableFeatures(NucSeq)
NucSeq <- ScaleData(NucSeq)
NucSeq <- RunPCA(NucSeq)
NucSeq <- FindNeighbors(NucSeq, dims = 1:30)
NucSeq <- FindClusters(NucSeq)
NucSeq <- RunUMAP(NucSeq, dims = 1:30)

DimPlot(NucSeq, reduction = "umap", group.by = "sample", raster =FALSE)
DimPlot(NucSeq, reduction = "umap", group.by = "seurat_clusters", raster =FALSE)

NucSeq <- RunAzimuth(NucSeq, reference = "data/reference/heart/")
DimPlot(NucSeq, reduction = "umap", group.by = "predicted.celltype.l1", raster =FALSE, cols = color_scheme)


#### 09 calculate cell composition ####
plt_cell_count<-NucSeq@meta.data[,c("sample", "predicted.celltype.l1")]
colnames(plt_cell_count)<-c("Sample", "Celltype")

ggplot(plt_cell_count,aes(x=Sample,fill=Celltype))+
  geom_bar(position="fill")+
  scale_fill_manual(values = color_scheme) +
  scale_x_discrete(labels = c("AF10" = "Donor 1", 
                              "AF25" = "Donor 2", 
                              "AF64" = "Donor 3", 
                              "AF65" = "Donor 4")) +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())
