##Cell type assignment

library(HGNChelper)

# load gene set preparation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# load cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

# DB file
db_ <- "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx"

tissue <- "Brain"


# prepare gene sets
gs_list <- gene_sets_prepare(db_, tissue)


# check Seurat object version (scRNA-seq matrix extracted differently in Seurat v4/v5)
seurat_package_v5 <- isFALSE('counts' %in% names(attributes(d801205_umap[["RNA"]])));
print(sprintf("Seurat object %s is used", ifelse(seurat_package_v5, "v5", "v4")))

# extract scaled scRNA-seq matrix
scRNAseqData_scaled <- if (seurat_package_v5) as.matrix(d801205_umap[["RNA"]]$scale.data) else as.matrix(d801205_umap[["RNA"]]@scale.data)


# run ScType
es.max <- sctype_score(scRNAseqData = scRNAseqData_scaled, scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)


# merge by cluster
cL_resutls <- do.call("rbind", lapply(unique(d801205_umap@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(d801205_umap@meta.data[d801205_umap@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(d801205_umap@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores <- cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores) 


d801205_umap@meta.data$sctype_classification = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  d801205_umap@meta.data$sctype_classification[d801205_umap@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

DimPlot(d801205_umap, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'sctype_classification')
