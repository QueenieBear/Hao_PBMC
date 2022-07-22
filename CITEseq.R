setwd("~/Summer2022/Hao_PBMC")
library(Seurat)
library(ggplot2)
library(sctransform)

#RNA data
# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "rawdata/RNA/")
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data,min.cells = 10, min.features = 200) 

# create a new assay to store ADT information
pbmc.adt <- Read10X(data.dir = "rawdata/ADT/")
adt_assay <- CreateAssayObject(counts = pbmc.adt)
# add this assay to the previously created Seurat object
pbmc[["ADT"]] <- adt_assay


# create a new assay to store HTO information
pbmc.hto <- Read10X(data.dir = "rawdata/HTO/")
hto_assay <- CreateAssayObject(counts = pbmc.hto)
# add this assay to the previously created Seurat object
pbmc[["HTO"]] <- hto_assay



# Validate that the object now contains multiple assays
Assays(pbmc)
# Extract a list of features measured in the ADT/HTO assay
rownames(pbmc[["ADT"]])
rownames(pbmc[["HTO"]])


#QC ------------------------------------------------------------------------------------------------------------

# Compute percent mito ratio
pbmc$mitoRatio <- PercentageFeatureSet(object = pbmc, pattern = "^MT-")
pbmc$mitoRatio<-pbmc@meta.data$mitoRatio/100



pdf("figure/Vln_nGene&nUMImitoRatio.pdf", width = 50, height = 20)
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA","nCount_ADT","nFeature_ADT", "mitoRatio"), 
        ncol = 2,pt.size =0)
dev.off()


#filter at cell level-------------------------------------------------------------------------------------
# Filter out low quality reads using selected thresholds 
filtered_pbmc <- subset(x = pbmc, 
                        subset= (nFeature_RNA <=6000 ) &
                          (nFeature_RNA >=500 ) &
                          (nCount_RNA >= 200) & 
                          (nCount_ADT <= 50000) & 
                          (mitoRatio <= 0.1))

#filter at gene level-------------------------------------------------------------------------------------
# default assay is RNA
DefaultAssay(pbmc) <- "RNA"
# Extract counts
counts <- GetAssayData(object = filtered_pbmc, slot = "counts")
# Output a logical vector for every gene on whether the more than zero counts per cell
nonzero <- counts > 0
# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
keep_genes <- Matrix::rowSums(nonzero) >= 10
# Only keeping those genes expressed in more than 10 cells
filtered_counts <- counts[keep_genes, ]

# Reassign to filtered Seurat object
filtered_PBMC<- CreateSeuratObject(filtered_counts, meta.data = filtered_pbmc@meta.data)
filtered_PBMC[["ADT"]]<-filtered_pbmc[["ADT"]]
filtered_PBMC[["HTO"]]<-filtered_pbmc[["HTO"]]



##reassess QC matrix-------------------------------------------------------------------------------------
pdf("figure/Vln_afterQC.pdf", width = 50, height = 20)
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA","nCount_ADT","nFeature_ADT", "mitoRatio"), 
        ncol = 2,pt.size =0)
dev.off()

metadata_clean <- filtered_PBMC@meta.data
pdf("figure/nGene&nCell_afterQC.pdf", width = 40, height = 7)
# Visualize the number of cell counts per sample
p1<-metadata_clean %>% 
  ggplot(aes(x=orig.ident, fill=orig.ident)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("nCells")

# Visualize the distribution of genes detected per cell via boxplot
p2<-metadata_clean %>% 
  ggplot(aes(x=orig.ident, y=nFeature_RNA, fill=orig.ident)) + 
  geom_boxplot() + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("nGene per cell")

p1+p2
dev.off()


pdf("figure/density_afterQC.pdf", width = 40, height = 5)
# Visualize the number UMIs/transcripts per cell
p3<-metadata_clean %>% 
  ggplot(aes(color=orig.ident, x=nCount_RNA, fill= orig.ident)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") 

# Visualize the distribution of genes detected per cell via histogram
p4<-metadata_clean %>% 
  ggplot(aes(color=orig.ident, x=nFeature_RNA, fill= orig.ident)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = c(200,5000),linetype="dashed")

# Visualize the distribution of mitochondrial gene expression detected per cell
p5<-metadata_clean %>% 
  ggplot(aes(color=orig.ident, x=mitoRatio, fill=orig.ident)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.1, linetype="dashed")
p3+p4+p5
dev.off()




#integration -------------------------------------------------------------------------------
#seperate samples
#reading HTO data
yhto <- as.data.frame(filtered_PBMC@assays$HTO@data)
rownames(yhto)#24samples
table(colSums(yhto)== 0)
pbmc.htos <- yhto[colSums(yhto)>0]
dim(pbmc.htos)

#take intersection
pbmc.umis <- filtered_PBMC@assays$RNA@data
pbmc.adts <- filtered_PBMC@assays$ADT@data

joint.bcs <- intersect(colnames(pbmc.umis), colnames(pbmc.htos))
length(joint.bcs)# 153,956 cells
joint.bcs <- intersect(joint.bcs, colnames(pbmc.adts))
length(joint.bcs)# 153,956 cells
length(colnames(pbmc.umis))#153956
length(colnames(pbmc.htos))#153956
length(colnames(pbmc.adts))#153956
all.equal(colnames(pbmc.umis), colnames(pbmc.htos))
all.equal(colnames(pbmc.umis), colnames(pbmc.adts))
#Subset RNA and HTO counts by joint cell barcodes
#pbmc.umis <- pbmc.umis[, joint.bcs]

saveRDS(filtered_PBMC,"metadata/filtered_PBMC.rds")

#normalization and scaling data
DefaultAssay(filtered_PBMC) <- "RNA"
# Normalize the counts
filtered_PBMC1 <- NormalizeData(filtered_PBMC)
# Score cells for cell cycle
filtered_PBMC1 <- CellCycleScoring(filtered_PBMC1, 
                                   g2m.features = cc.genes$g2m.genes, 
                                   s.features = cc.genes$s.genes)
filtered_PBMC1 <- FindVariableFeatures(filtered_PBMC1, selection.method = "mean.var.plot")
filtered_PBMC1 <- ScaleData(filtered_PBMC1, features = VariableFeatures(filtered_PBMC1))

# Normalize HTO data, here we use centered log-ratio (CLR) transformation
filtered_PBMC1 <- NormalizeData(filtered_PBMC1,assay = "HTO",normalization.method = "CLR")
filtered_PBMC1 <- NormalizeData(filtered_PBMC1, normalization.method = "CLR", margin = 2, assay = "ADT")

filtered_PBMC1 <- RunPCA(filtered_PBMC1, npcs=30,verbose = F)
ElbowPlot(filtered_PBMC1)
filtered_PBMC1 <- FindNeighbors(filtered_PBMC1, dims = 1:30, verbose = F)
filtered_PBMC1 <- FindClusters(filtered_PBMC1, verbose = F)
filtered_PBMC1 <- RunUMAP(filtered_PBMC1, dims = 1:30, verbose = F)
DimPlot(filtered_PBMC1,
        reduction = "umap",
        label = TRUE,
        label.size = 5)
saveRDS(filtered_PBMC1,"metadata/filtered_PBMC1.rds")


# Now, we will visualize CD3 levels for RNA and protein By setting the default assay, we can
# visualize one or the other
DefaultAssay(filtered_PBMC1) <- "ADT"
p1 <- FeaturePlot(filtered_PBMC1, "TRBC1", cols = c("lightgrey", "darkgreen")) + ggtitle("TRBC1 protein")
DefaultAssay(filtered_PBMC1) <- "RNA"
p2 <- FeaturePlot(filtered_PBMC1, "TRBC1") + ggtitle("TRBC1 RNA")
# place plots side-by-side
p1 | p2

###HTODemux--------------------------------------------------------------------------------------------
filtered_PBMC1_batch1<-subset(filtered_PBMC1,orig.ident %in% c("L1","L2","L3","L4","L5"))
filtered_PBMC1_batch2<-subset(filtered_PBMC1,orig.ident %in% c("E2L1","E2L2","E2L3","E2L4","E2L5","E2L6","E2L7","E2L8"))
filtered_PBMC1_batch1@assays$HTO@data<-filtered_PBMC1_batch1@assays$HTO@data[1:12,]
filtered_PBMC1_batch2@assays$HTO@data<-filtered_PBMC1_batch2@assays$HTO@data[13:24,]
filtered_PBMC1_batch1 <- HTODemux(filtered_PBMC1_batch1,
                         assay = "HTO", 
                         positive.quantile = 0.99)
filtered_PBMC1_batch2 <- HTODemux(filtered_PBMC1_batch2,
                                  assay = "HTO", 
                                  positive.quantile = 0.99)

filtered_PBMC1_batch2
table(filtered_PBMC1_batch1$HTO_classification.global)
filtered_PBMC1_batch2@assays
table(filtered_PBMC1_batch1$hash.ID)

# Group cells based on the max HTO signal
Idents(filtered_PBMC1_batch2) <- "HTO_maxID"
pdf("figure/RidgePlot_HTO.pdf",width = 15,height = 15)
RidgePlot(filtered_PBMC1_batch2, assay = "HTO", features = rownames(filtered_PBMC1_batch2[["HTO"]])[1:12], ncol = 3)
dev.off()

#analyze HTODemux error without splitting data by batches---------------------------------------------------------
df<-as.data.frame(filtered_PBMC1@assays$HTO@counts)
df1<-df[ , apply(df, 2, function(x) any(x==0))]
#batch1
df2<-df[1:12,]
df2.1<-df2[ , apply(df2, 2, function(x) any(x==0))]
L1.0<-df2[,startsWith(colnames(df2),"L1")]
L1<-df2.1[,startsWith(colnames(df2.1),"L1")]
L1[apply(L1, c(1,2), function(x) any(x==0))]<-"a"
L1[apply(L1, c(1,2), function(x) any(x!=0&x!="a"))]<-"b"
L1.1<-data.matrix(L1)
L1.1[apply(L1.1, c(1,2), function(x) any(x==2))]<-0
L1.2<-rowSums(L1.1)
L1.2<-data.matrix(L1.2)
colnames(L1.2)<-"L1"

Batch1<-cbind(L1.2,L2.2,L3.2,L4.2,L5.2,E2L1.2,E2L2.2,E2L3.2,E2L4.2,E2L5.2,E2L6.2,E2L7.2,E2L8.2)
Batch1.1<-rbind(Batch1,c(ncol(L1.0),ncol(L2.0),ncol(L3.0),ncol(L4.0),ncol(L5.0),ncol(E2L1.0),ncol(E2L2.0),ncol(E2L3.0),ncol(E2L4.0),ncol(E2L5.0),ncol(E2L6.0),ncol(E2L7.0),ncol(E2L8.0)))
rownames(Batch1.1)[13]<-"TotalCellNumber"
Batch1_prop<-t(apply(Batch1.1, 1, function(x) x /Batch1.1[13,]))
Batch1_prop<-Batch1_prop[1:12,]
#round(Batch1_prop,2)





#batch2
df3<-df[13:24,]
df3.1<-df3[ , apply(df3, 2, function(x) any(x==0))]
L1.0<-df3[,startsWith(colnames(df3),"L1")]
L1<-df3.1[,startsWith(colnames(df3.1),"L1")]
L1[apply(L1, c(1,2), function(x) any(x==0))]<-"a"
L1[apply(L1, c(1,2), function(x) any(x!=0&x!="a"))]<-"b"
L1.1<-data.matrix(L1)
L1.1[apply(L1.1, c(1,2), function(x) any(x==2))]<-0
L1.2<-rowSums(L1.1)
L1.2<-data.matrix(L1.2)
colnames(L1.2)<-"L1"

Batch2<-cbind(L1.2,L2.2,L3.2,L4.2,L5.2,E2L1.2,E2L2.2,E2L3.2,E2L4.2,E2L5.2,E2L6.2,E2L7.2,E2L8.2)
Batch2.1<-rbind(Batch2,c(ncol(L1.0),ncol(L2.0),ncol(L3.0),ncol(L4.0),ncol(L5.0),ncol(E2L1.0),ncol(E2L2.0),ncol(E2L3.0),ncol(E2L4.0),ncol(E2L5.0),ncol(E2L6.0),ncol(E2L7.0),ncol(E2L8.0)))
rownames(Batch2.1)[13]<-"TotalCellNumber"
Batch2_prop<-t(apply(Batch2.1, 1, function(x) x /Batch2.1[13,]))
Batch2_prop<-Batch2_prop[1:12,]
round(Batch2_prop,2)

Batch1and2_prop<-rbind(Batch1_prop,Batch2_prop)
Batch1and2_prop<-round(Batch1and2_prop,3)
library(pheatmap)
pdf("figure/Cellprop_0count.pdf",width = 6.5,height = 8)
pheatmap(Batch1and2_prop, 
         scale = "none", 
         cluster_cols = F, 
         cluster_rows = F,
         angle_col = 45,
         fontsize_row = 8,
         fontsize_col = 8,
         main="Proportion of Cells with 0 count",
         border_color = NA,
         display_numbers = T,
         number_format = "%.3f",
         number_color="white",
         legend=F,
         #legend_breaks=c(0,10,20,50,100,200,400,1000,2000),
         #color = c(brewer.pal(9, "Blues")[9:5]),
         #breaks =c(0,10,20,50,100,200,400,1000,2000),
         cellwidth=30,
         cellheight=20)
dev.off()

Batch1and2_num<-rbind(Batch1.1,Batch2.1)
Batch1and2_num<-Batch1and2_num[-13,]
pdf("figure/Cellnum_0count.pdf",width = 6.5,height = 8)
pheatmap(Batch1and2_num, 
         scale = "none", 
         cluster_cols = F, 
         cluster_rows = F,
         angle_col = 45,
         fontsize_row = 8,
         fontsize_col = 8,
         main="Number of Cells with 0 count",
         border_color = NA,
         display_numbers = T,
         number_format = "%.0f",
         number_color="white",
         legend=F,
         #legend_breaks=c(0,10,20,50,100,200,400,1000,2000),
         #color = c(brewer.pal(9, "Blues")[9:5]),
         #breaks =c(0,50,100,200,300,400,500,600,700,800,10000),
         cellwidth=30,
         cellheight=20)
dev.off()




pdf("figure/UMAP.pdf", width = 90, height = 5)
DimPlot(filtered_PBMC1,
        reduction = "umap",
       # split.by = "orig.ident",
        label = F,
        label.size = 0
        )
dev.off()



#subset data------------------------------------------------------------------------------------------------------------------------------
filtered_PBMC2<-subset(filtered_PBMC1_batch2,
                       subset = 
                       (orig.ident %in% c("E2L1","E2L2","E2L3")&
                       (HTO_classification.global=="Singlet")))


# Normalize the counts
DefaultAssay(filtered_PBMC2) <- "RNA"
filtered_PBMC2 <- NormalizeData(filtered_PBMC2)
# run sctransform
filtered_PBMC2 <- SCTransform(filtered_PBMC2, vars.to.regress = c("mitoRatio","S.Score", "G2M.Score"), verbose = F)
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
filtered_PBMC2 <- NormalizeData(filtered_PBMC2,assay = "HTO",normalization.method = "CLR")
filtered_PBMC2 <- NormalizeData(filtered_PBMC2, normalization.method = "CLR", margin = 2, assay = "ADT")

# These are now standard steps in the Seurat workflow for visualization and clustering
filtered_PBMC2 <- RunPCA(filtered_PBMC2, npcs=30,verbose = F)
ElbowPlot(filtered_PBMC2)
filtered_PBMC2 <- FindNeighbors(filtered_PBMC2, dims = 1:30, verbose = F)
filtered_PBMC2 <- FindClusters(filtered_PBMC2, verbose = F)
filtered_PBMC2 <- RunUMAP(filtered_PBMC2, dims = 1:30, verbose = F)
DimPlot(pbmc,
        reduction = "umap",
        label = TRUE,
        label.size = 5)









#doublet deletion---------------------------------------------------------------------------------------
library(DoubletFinder)

# pK Identification (no ground-truth) 
sweep.res.list_PBMC <- paramSweep_v3(filtered_PBMC2, PCs = 1:30, sct = T)
sweep.stats_PBMC <- summarizeSweep(sweep.res.list_PBMC, GT = FALSE)
bcmvn_PBMC <- find.pK(sweep.stats_PBMC)

ggplot(bcmvn_PBMC,aes(pK,BCmetric,group=1))+
  geom_point()+
  geom_line()+
  theme(axis.text.x=element_text(angle = 45))

library(dplyr)
pK<-bcmvn_PBMC %>%
  dplyr::filter(BCmetric == max(BCmetric)) %>%
  select(pK)

pK<-as.numeric(as.character(pK[[1]]))

# Homotypic Doublet Proportion Estimate 
annotations <- filtered_PBMC2@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)           
nExp_poi <- round(0.076*nrow(filtered_PBMC2@meta.data))  #7.6% is from https://uofuhealth.utah.edu/huntsman/shared-resources/gba/htg/single-cell/genomics-10x.php
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder with varying classification stringencies 
filtered_PBMC2 <- doubletFinder_v3(filtered_PBMC2, PCs = 1:30, pN = 0.25, pK = pK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = T)


#visualize doublets
names(filtered_PBMC2@meta.data)
DimPlot(filtered_PBMC2,reduction = "umap",group.by = "DF.classifications_0.25_0.12_1812")
VlnPlot(filtered_PBMC2, features = "nFeature_RNA", group.by = "DF.classifications_0.25_0.12_1812", pt.size = 0)

#number of singlet and doublet---------------------------------------------------------------------------------------
table(filtered_PBMC2@meta.data$DF.classifications_0.25_0.12_1812)

metadata <- filtered_PBMC2@meta.data
names(metadata)
colnames(metadata)[24] <- "doublet_finder"
filtered_PBMC2@meta.data <- metadata 

# subset and save
filtered_PBMC.singlets <- subset(filtered_PBMC2, doublet_finder == "Singlet")##whether do I need to rescale de subset data?
filtered_PBMC2
filtered_PBMC.singlets
DimPlot(filtered_PBMC.singlets,reduction = "umap")

saveRDS(filtered_PBMC.singlets, "metadata/filtered_PBMC_singlets.rds")

# featureplot the second density peak
metadata_clean <- filtered_PBMC.singlets@meta.data
metadata_clean %>% 
  ggplot(aes(color=orig.ident, x=nFeature_RNA, fill= orig.ident)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = c(500,6000),linetype="dashed")

FeaturePlot(filtered_PBMC.singlets, 
            reduction = "umap", 
            features ="G2M.Score", 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = F)

VlnPlot(filtered_PBMC.singlets, features = c("mitoRatio"),pt.size = 0)



