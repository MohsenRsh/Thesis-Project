


#Now, with Azimuth, I annotate the cells. I will use the Azimuth package to annotate the cells in the filtered Seurat 
#object. The Azimuth package provides a reference-based cell type annotation method for single-cell RNA-seq data.
library(Azimuth)


# If you haven't run SCTransform yet, do it first
obj1205_sct <- SCTransform(str1205ForAzimuth)



# The RunAzimuth function can take a Seurat object as input
azd801205 <- RunAzimuth(str1205ForAzimuth, reference = "humancortexref")

azd801205[["predicted.class"]]


# Plot the UMAP with Azimuth annotations
# The `DimPlot` function is used to visualize the UMAP with cell type annotations.
#all libraries are loaded
DimPlot(azd801205, group.by = "predicted.class", label = TRUE) +
  ggtitle("Azimuth Cell Type Annotations") +
  theme_minimal()

DimPlot(azd801205, group.by = "", label = TRUE) +
  ggtitle("Azimuth Cell Type Annotations") +
  theme_minimal()

azl1 <- DimPlot(azd801205, group.by = "predicted.cluster", label = TRUE, label.size = 3) 
azl1

sceAz <- as.SingleCellExperiment(azd801205)


azd801205[[]]

glimpse(azd801205[[]])

#Add a 'sampleID' column to the metadata but do it randomely and assign three different samples
azd801205$sampleID <- sample(c("1205-1", "1205-2", "1205-3"), size = ncol(azd801205), replace = TRUE)

#removing the last part of the rownames
azd801205 <- RenameCells(azd801205, new.names = gsub("-1", "", colnames(azd801205)))

# saving as matrix
azd801205 <- as.matrix(azd801205@assays$RNA@counts)


sceAz <- as.SingleCellExperiment(azd801205)

#save the object
saveRDS(azd801205, file = "D:/Thesis/Analysis/scData/1205_d80/azd801205.rds")



