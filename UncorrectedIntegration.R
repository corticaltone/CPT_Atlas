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

ifnb.list <- list()
samples <- list()
data.folder = "/omics/odcf/analysis/OE0519_projects/chptumor/NewAnalysis/Omega/MappingTest/Individual/"
#load seurat objects
for (tumor in tumors){
  file_loc <- paste0(data.folder,tumor,".Robj")
  data <- load(file_loc)
  if(ncol(object) > 6000) {
    object <- subset(object, cells = sample(Cells(object), 6000)) #limit sample contribution to 6K cells
  }

  ifnb.list <- append(ifnb.list, object)
  print(tumor)
}


# normalize and identify variable features for each dataset independently
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

choroid.combined <- merge(ifnb.list[[1]], y = c(ifnb.list[[2]], ifnb.list[[3]], ifnb.list[[4]], ifnb.list[[5]],ifnb.list[[6]], ifnb.list[[7]], ifnb.list[[8]], ifnb.list[[9]], ifnb.list[[10]], ifnb.list[[11]], ifnb.list[[12]], ifnb.list[[13]],ifnb.list[[14]], ifnb.list[[15]]), add.cell.ids = c("a", "b", "c","d", "e", "f","g", "h", "i","j", "k", "l","m", "n", "o"), project = "PBMC15K")


# Run the standard workflow for visualization and clustering
choroid.combined <- ScaleData(choroid.combined, verbose = FALSE)
choroid.combined <- FindVariableFeatures(choroid.combined)
choroid.combined <- RunPCA(choroid.combined, npcs = 30, verbose = FALSE)
choroid.combined <- RunUMAP(choroid.combined, reduction = "pca", dims = 1:30)
save(choroid.combined, file = "CC.Robj")
choroid.combined <- RunTSNE(choroid.combined, reduction = "pca", dims = 1:30)
choroid.combined <- FindNeighbors(choroid.combined, reduction = "pca", dims = 1:30)
choroid.combined <- FindClusters(choroid.combined, resolution = 0.5)
choroidUC.combined <- choroid.combined
save(choroidUC.combined, file = "UC.Robj")
# Visualization

p1 <- DimPlot(choroid.combined, reduction = "umap", group.by = "type")
p2 <- DimPlot(choroid.combined, reduction = "umap", label = TRUE, repel = TRUE)
png(filename = paste0(experiment, "DimplotUC.png"))
p1 + p2
dev.off()

png("UCMarkers.png")
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

celltypes <- levels(Idents(choroid.combined))

for (cell in celltypes) {
  marks <- FindMarkers(choroid.combined, ident.1 = cell)
  write.table(marks, file = paste0("CelltypeMarkers-", cell, ".csv"), sep = )
}


markers <- c("OTX2", "LMX1A", "PRLR", "SLC4A10", "ENPP2")
DimPlot(choroid.combined, label = TRUE, group.by = choroid.combined$m)

table(choroid.combined$sample)
choroid.combined <- RenameIdents(choroid.combined, 'ChP1' = "ChP", 'ChP2' = "ChP", 'ChP3' = "ChP", 'ChP4' = "ChP", 'CP4' = "adult", 'CP1' = "adult", 'CP2' = "pedB", 'CP3' = "adult", 'CC1' = "pedB", 'CC2' = "pedB", 'CC3' = "pedB", 'CC4' = "pedB", 'CC5' = "pedB", 'aCP1' = "pedB", 'aCP2' = "pedB")
DimPlot(choroid.combined, label = FALSE, reduction = "tsne", cols = c("#33FF00", "#E69F00", "#56B4E9")) + theme(legend.position = "bottom")


