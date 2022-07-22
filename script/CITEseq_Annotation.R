setwd("~/Summer2022/Hao_PBMC")

library(Seurat)
pbmc <- readRDS("~/Summer2022/Hao_PBMC/metadata/filtered_PBMC_singlets.rds")

FeaturePlot(pbmc, 
            reduction = "umap", 
            features ="CD3D", 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = F)

DimPlot(pbmc,
        reduction = "umap",
        label = TRUE,
        label.size = 4)



###subset T cells-------------------------------------------------------------------------------------
DefaultAssay(pbmc)<-"SCT"
# Determine the clusters for various resolutions  
pbmc <- RunPCA(pbmc, npcs=40,verbose = F)
pbmc <- FindNeighbors(object = pbmc, 
                      dims = 1:40)
pbmc <- FindClusters(object = pbmc,
                     resolution = 1.5)
# Plot the UMAP
DimPlot(pbmc,
        reduction = "umap",
        label = TRUE,
        label.size = 4)

# visualize one or the other
DefaultAssay(pbmc) <- "ADT"
FeaturePlot(pbmc, "CD3-1", cols = c("lightgrey", "darkgreen")) 
DefaultAssay(pbmc) <- "SCT"
FeaturePlot(pbmc, "CD3G") 



pbmc_T<-subset(pbmc,idents = c(20,16,3,11,26,22,0,7,1,31,9,13))
# Normalize the counts
DefaultAssay(pbmc_T) <- "RNA"
pbmc_T <- NormalizeData(pbmc_T)
# run sctransform
pbmc_T <- SCTransform(pbmc_T, vars.to.regress = c("mitoRatio","S.Score", "G2M.Score"), verbose = F)
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
pbmc_T <- NormalizeData(pbmc_T,assay = "HTO",normalization.method = "CLR")
pbmc_T <- NormalizeData(pbmc_T, normalization.method = "CLR", margin = 2, assay = "ADT")

# These are now standard steps in the Seurat workflow for visualization and clustering
pbmc_T <- RunPCA(pbmc_T, npcs=30,verbose = F)
#ElbowPlot(pbmc_T)
pbmc_T <- FindNeighbors(pbmc_T, dims = 1:30, verbose = F)
pbmc_T <- FindClusters(pbmc_T, resolution = 2)
pbmc_T <- RunUMAP(pbmc_T, dims = 1:30, verbose = F)
DimPlot(pbmc_T,
        reduction = "umap",
        label = TRUE,
        label.size = 5)


#adjust the resolution
DefaultAssay(pbmc_T)<-"SCT"
pbmc_T <- RunPCA(pbmc_T, npcs=30,verbose = F)
pbmc_T <- FindNeighbors(pbmc_T, dims = 1:30, verbose = F)
pbmc_T <- FindClusters(pbmc_T, resolution = 3)
DimPlot(pbmc_T,
        reduction = "umap",
        label = TRUE,
        label.size = 5)


#subset CD4T cells----------------------------------------------------------------------------
pbmc_CD4T<-subset(pbmc_T,idents = c(0,1,2,4,5,6,8,9,11,13,24,28,30)) 
DimPlot(pbmc_CD4T,
        reduction = "umap",
        label = TRUE,
        label.size = 5)

# Normalize the counts
DefaultAssay(pbmc_CD4T) <- "RNA"
pbmc_CD4T <- NormalizeData(pbmc_CD4T)
# run sctransform
pbmc_CD4T <- SCTransform(pbmc_CD4T, vars.to.regress = c("mitoRatio","S.Score", "G2M.Score"), verbose = F)
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
pbmc_CD4T <- NormalizeData(pbmc_CD4T,assay = "HTO",normalization.method = "CLR")
pbmc_CD4T <- NormalizeData(pbmc_CD4T, normalization.method = "CLR", margin = 2, assay = "ADT")

# These are now standard steps in the Seurat workflow for visualization and clustering
pbmc_CD4T <- RunPCA(pbmc_CD4T, npcs=30,verbose = F)
ElbowPlot(pbmc_CD4T)
pbmc_CD4T <- FindNeighbors(pbmc_CD4T, dims = 1:30, verbose = F)
pbmc_CD4T <- FindClusters(pbmc_CD4T, resolution = 2)
pbmc_CD4T <- RunUMAP(pbmc_CD4T, dims = 1:30, verbose = F)
DimPlot(pbmc_CD4T,
        reduction = "umap",
        label = TRUE,
        label.size = 5)






# visualize one or the other
DefaultAssay(pbmc_CD4T) <- "ADT"
FeaturePlot(pbmc_CD4T, "CD69", cols = c("lightgrey", "darkgreen")) 
DefaultAssay(pbmc_CD4T) <- "SCT"
FeaturePlot(pbmc_CD4T, c("CLIC3")) 




#umap based on ADT & RNA(WNN) ------------------------------------------------------------------------------------------
DefaultAssay(pbmc_CD4T) <- 'ADT'
# we will use all ADT features for dimensional reduction
# we set a dimensional reduction name to avoid overwriting the 
library(dplyr)
VariableFeatures(pbmc_CD4T) <- rownames(pbmc_CD4T[["ADT"]])
pbmc_CD4T <- NormalizeData(pbmc_CD4T, normalization.method = 'CLR', margin = 2) %>% 
  SCTransform() %>% RunPCA(reduction.name = 'apca')

pbmc_CD4T <- FindMultiModalNeighbors(
  pbmc_CD4T, reduction.list = list("pca", "apca"), 
  dims.list = list(1:30, 1:18), modality.weight.name = "RNA.weight"
)

#WNN
pbmc_CD4T <- RunUMAP(pbmc_CD4T, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
pbmc_CD4T <- FindClusters(pbmc_CD4T, graph.name = "wsnn", algorithm = 3, resolution = 2, verbose = FALSE)
DimPlot(pbmc_CD4T, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 2.5) 

#RNA
pbmc_CD4T <-RunUMAP(pbmc_CD4T, reduction = 'pca', dims = 1:30, assay = 'RNA', 
              reduction.name = 'rna.umap', reduction.key = 'rnaUMAP_')
DimPlot(pbmc_CD4T, reduction = 'rna.umap', label = TRUE, 
        repel = TRUE, label.size = 2.5)

#ADT
pbmc_CD4T <-RunUMAP(pbmc_CD4T, reduction = 'apca', dims = 1:18, assay = 'ADT', 
              reduction.name = 'adt.umap', reduction.key = 'adtUMAP_')
DimPlot(pbmc_CD4T, reduction = 'adt.umap', label = TRUE, 
        repel = TRUE, label.size = 2.5)

#----------------------------------------------------------------------
FeaturePlot(pbmc_CD4T, features = c("adt_CD45RA","adt_CD16","adt_CD161"),
                  reduction = 'wnn.umap', max.cutoff = 2, 
                  cols = c("lightgrey","darkgreen"), ncol = 3)
FeaturePlot(pbmc_CD4T, features = c("rna_TRDC","rna_MPO","rna_AVP"), 
                  reduction = 'wnn.umap', max.cutoff = 3, ncol = 3)











Tmarkers<-c("CD3D","CD3E","CD3G","CD4","CD8A","CD8B","TRAC","TRBC2","TRBC1",
  "CCR7","FHIT","LEF1","SELL",
  "IL7R","TMSB10","ITGB1","LTB","AQP3","MAL","CCL5","GZMK","GZMH","FGFBP2","GNLY","NKG7","LINC02446","ANXA1","CXCR4","KLRB1","GZMA","TXN","TNFRSF4","TNFRSF1A","TNFRSF1B","TNFRSF6B","	TNFRSF8","LYAR",
  "NKG7","GNLY","CCL4","CCL5","GZMB","GZMH","IL32",
  "CD44","CD69",
  "RTKN2","FOXP3","AC133644.2","IL2RA","CCL5","KLRD1","GZMK","CST7","TRGC2","TIGIT","CTLA4",
  "KLRB1","NKG7","GZMK","SLC4A10","NCR3",
  "TRDC","TRGC1","KLRC1","NKG7")


