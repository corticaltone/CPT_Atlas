library(Seurat)
#library(SeuratData)
library(patchwork)
library(ggplot2)
library(cowplot)
library(dplyr)

#metadata setup
# Load the ChP datasets

experiment <- "Cars"

tumors <- c("C1", "C3", "C4", "C5", "T1", "T3", "T4", "T5", "T6", "T7", "T8", "T9", "T10", "T12", "T13")
#tumors <- c("T5", "T6", "T7", "T12", "T9")

ifnb.list <- list()

data.folder = "/omics/odcf/analysis/OE0519_projects/chptumor/NewAnalysis/Omega/MappingTest/Individual/"
#load seurat objects
for (tumor in tumors){
  file_loc <- paste0(data.folder,tumor,".Robj")
  data <- load(file_loc)

  ifnb.list <- append(ifnb.list, object)
  print(tumor)
}


# normalize and identify variable features for each dataset independently
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = ifnb.list)

#create integrated object
choroid.anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features)
choroid.combined <- IntegrateData(anchorset = choroid.anchors)


# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(choroid.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
choroid.combined <- ScaleData(choroid.combined, verbose = FALSE)
choroid.combined <- RunPCA(choroid.combined, npcs = 30, verbose = FALSE)
choroid.combined <- RunUMAP(choroid.combined, reduction = "pca", dims = 1:30)
choroid.combined <- RunTSNE(choroid.combined, reduction = "pca", dims = 1:30)
choroid.combined <- FindNeighbors(choroid.combined, reduction = "pca", dims = 1:30)
choroid.combined <- FindClusters(choroid.combined, resolution = 0.5)
save(choroid.combined, file = "choroid.combined.Robj")
# Visualization

p1 <- DimPlot(choroid.combined, reduction = "umap", group.by = "type")
p2 <- DimPlot(choroid.combined, reduction = "umap", label = TRUE, repel = TRUE)
png(filename = paste0(experiment, "DimplotPap.png"))
p1 + p2
dev.off()


VlnPlot(choroid.combined, features = c("OTX2", "CD163", "CD247", "CD79A", "MKI67", "VWF",  "PECAM1"), split.by = "type", group.by = "seurat_clusters", pt.size = 0, combine = FALSE)
dev.off()





cn <- length(table(Car.combined$seurat_clusters)) -1
for (cluster in 0:cn) {
  markers <- FindMarkers(Car.combined, ident.1 = cluster)
  write.table(markers, file = paste0("Cluster", cluster, "markersCar.csv"), sep = ",")
}

png(filename = paste0(experiment, "Vln.png"))
VlnPlot(choroid.combined, features = c(c("CD247", "PTPRC", "CD163", "GFAP", "RBFOX3", "MAP2", "PECAM1", "VWF", "PLVAP", "ADAM12", "COL3A1", "ANGPT1", "VCAM1", "AQP1", "OTX2", "KCNJ13", "ENPP2", "HTR2C", "FOXJ1", "GMNC", "ASPM", "MKI67")), pt.size = 0, stack = T)

dev.off()



Cars <- subset(Car.combined, subset = seurat_clusters %in% c(0,1,2,3,5,6, 7,8,9, 11))
Cars <- RenameIdents(Cars, '0' = "Epithelial-2", '1' = "Epithelial-2", '2' = "Epithelial-2", '3' = "Macrophage", '5' = "Endothelial", '6' = "Mesenchymal", '7' = "Proliferative", '8' = "Epithelial-1", '9' = "Epithelial-1", '11' = "Tcell")

DimPlot(Cars, reduction = "tsne", cols = c("brown1", "orange1", "darkgreen"))


