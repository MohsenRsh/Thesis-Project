#### Visium HD Manuscript
### 

## METHOD: Single Cell Flex Analysis

## Load required Packages
library(Seurat)
# library(Azimuth)
# library(BPCells)

# Load Auxiliary Functions
source("AuxFunctions.R")

# Read Aggr'd dataset as a Seurat V5 object
FlexOutPath <- "crc_dataset/Chromium_sc_flex/outs/" # Path to cellranger aggr output folder
Flex.mat <- Read10X_h5(filename = paste0(FlexOutPath,"HumanColonCancer_Flex_Multiplex_count_filtered_feature_bc_matrix.h5")) 
# Flex.mat <- ConvertEnsembleToSymbol(mat = Flex.mat, species = "human") 

# Read MetaData.csv file to be used as MetaData ( Patient, etc)
MetaData<-read.csv(paste0(FlexOutPath,"SingleCell_MetaData.csv"))

rownames(MetaData)<-MetaData$Barcode

# Create Seurat Object
ColonCancer_Flex<-CreateSeuratObject(Flex.mat,meta.data = MetaData)

# Add % MT
ColonCancer_Flex[["MT.percent"]] <- PercentageFeatureSet(ColonCancer_Flex, pattern = "^MT-")

# Add UMAP projection
ColonCancer_Flex[["umap"]] <- CreateDimReducObject(embeddings = as.matrix(MetaData[,c("UMAP1","UMAP2")]), key = "UMAP_", assay = DefaultAssay(ColonCancer_Flex))


## NOTA BENE: i "remove" secondo i metadata degli autori sono 19103
##            i "remove" runnando il LORO SCRIPT con le LORO THRESHOLD sono 19070

# continuo usando i loro risultati
ColonCancer_Flex$QCFilter <- MetaData$QCFilter

# Remove bcs that failed QC
ColonCancer_Flex<-subset(ColonCancer_Flex,cells=colnames(ColonCancer_Flex)[ColonCancer_Flex$QCFilter=="Keep"])
# MetaData <- MetaData[MetaData$QCFilter=="Keep",]

# SCtrasform for normalization TROPPO GRANDE
# library(sctransform)
# ColonCancer_Flex <- SCTransform(ColonCancer_Flex, vars.to.regress = "MT.percent",
#                                 new.assay.name = "SCT")

# Seurat V5 Processing (We will sketch and then project) at ~ 15% of the full data
ColonCancer_Flex <- NormalizeData(ColonCancer_Flex)
ColonCancer_Flex <- FindVariableFeatures(ColonCancer_Flex)

ColonCancer_Flex <- SketchData(object = ColonCancer_Flex,ncells = 37000,
                               features = VariableFeatures(ColonCancer_Flex),
                               method = "LeverageScore",sketched.assay = "sketch")

DefaultAssay(ColonCancer_Flex) <- "sketch"
ColonCancer_Flex <- FindVariableFeatures(ColonCancer_Flex)
ColonCancer_Flex <- ScaleData(ColonCancer_Flex)
ColonCancer_Flex <- RunPCA(ColonCancer_Flex)
# ElbowPlot(ColonCancer_Flex,ndims=40)
# ColonCancer_Flex <- FindNeighbors(ColonCancer_Flex, dims = 1:25)
# ColonCancer_Flex <- FindClusters(ColonCancer_Flex, resolution = 0.6)

# Find Markers for level 2 clusters
ColonCancer_Flex<-SetIdent(ColonCancer_Flex,value = "Level2")
markers_level2<-FindAllMarkers(ColonCancer_Flex,min.diff.pct = 0.2,logfc.threshold = 0.2,only.pos = T)
markers_level2<-markers_level2[markers_level2$p_val_adj<0.05,]

# Find Markers for level 1 clusters
ColonCancer_Flex<-SetIdent(ColonCancer_Flex,value = "Level1")
markers_level1<-FindAllMarkers(ColonCancer_Flex,min.diff.pct = 0.2,logfc.threshold = 0.2,only.pos = T)
markers_level1<-markers_level1[markers_level1$p_val_adj<0.05,]



# Project to Full dataset

ColonCancer_Flex <- Seurat:::ProjectData(
  object = ColonCancer_Flex,
  assay = "RNA",
  full.reduction = "pca.full",
  sketched.assay = "sketch",
  sketched.reduction = "pca",
  umap.model = "umap",
  dims = 1:50,
  refdata = list(ProjectAll_L1 = "Level1",ProjectAll_L2 = "Level2")
)

