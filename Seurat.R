library(Seurat)
library(cowplot)
library(Matrix)
library(dplyr)
library(ggplot2)
library(reshape2)

# data ----
E1<-Read10X(data.dir = "./data/E1_count/outs/filtered_feature_bc_matrix")
E2<-Read10X(data.dir = "./data/E2_count/outs/filtered_feature_bc_matrix")
E3<-Read10X(data.dir = "./data/E3_count/outs/filtered_feature_bc_matrix")
E4<-Read10X(data.dir = "./data/E4_count/outs/filtered_feature_bc_matrix")
N1<-Read10X(data.dir = "./data/N1_count/outs/filtered_feature_bc_matrix")
N2<-Read10X(data.dir = "./data/N2_count/outs/filtered_feature_bc_matrix")
N3<-Read10X(data.dir = "./data/N3_count/outs/filtered_feature_bc_matrix")
N4<-Read10X(data.dir = "./data/N4_count/outs/filtered_feature_bc_matrix")
N5<-Read10X(data.dir = "./data/N5_count/outs/filtered_feature_bc_matrix")

E1<-CreateSeuratObject(counts = E1, project = "E1",min.cells = 3, min.features = 200)
E2<-CreateSeuratObject(counts = E2, project = "E2",min.cells = 3, min.features = 200)
E3<-CreateSeuratObject(counts = E3, project = "E3",min.cells = 3, min.features = 200)
E4<-CreateSeuratObject(counts = E4, project = "E4",min.cells = 3, min.features = 200)
N1<-CreateSeuratObject(counts = N1, project = "N1",min.cells = 3, min.features = 200)
N2<-CreateSeuratObject(counts = N2, project = "N2",min.cells = 3, min.features = 200)
N3<-CreateSeuratObject(counts = N3, project = "N3",min.cells = 3, min.features = 200)
N4<-CreateSeuratObject(counts = N4, project = "N4",min.cells = 3, min.features = 200)
N5<-CreateSeuratObject(counts = N5, project = "N5",min.cells = 3, min.features = 200)

E1$group<-"POP"
E2$group<-"POP"
E3$group<-"POP"
E4$group<-"POP"
N1$group<-"CTRL"
N2$group<-"CTRL"
N3$group<-"CTRL"
N4$group<-"CTRL"
N5$group<-"CTRL"
N1$sample <- "N1"
N2$sample <- "N2"
N3$sample <- "N3"
N4$sample <- "N4"
N5$sample <- "N5"
E1$sample <- "E1"
E2$sample <- "E2"
E3$sample <- "E3"
E4$sample <- "E4"

E1[["percent.mt"]]<-PercentageFeatureSet(E1,pattern = '^ACOX1$|^ND1$|^ND2$|^COX1$|^COX2$|^ATP8$|^ATP6$|^COX3$|^ND3$|^ND4L$|^ND4$|^ND5$|^ND6$|^CYTB$')
P.E1<-VlnPlot(E1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
E2[["percent.mt"]]<-PercentageFeatureSet(E2,pattern = "^ACOX1$|^ND1$|^ND2$|^COX1$|^COX2$|^ATP8$|^ATP6$|^COX3$|^ND3$|^ND4L$|^ND4$|^ND5$|^ND6$|^CYTB$")
P.E2<-VlnPlot(E2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
E3[["percent.mt"]]<-PercentageFeatureSet(E3,pattern = "^ACOX1$|^ND1$|^ND2$|^COX1$|^COX2$|^ATP8$|^ATP6$|^COX3$|^ND3$|^ND4L$|^ND4$|^ND5$|^ND6$|^CYTB$")
P.E3<-VlnPlot(E3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
E4[["percent.mt"]]<-PercentageFeatureSet(E4,pattern = "^ACOX1$|^ND1$|^ND2$|^COX1$|^COX2$|^ATP8$|^ATP6$|^COX3$|^ND3$|^ND4L$|^ND4$|^ND5$|^ND6$|^CYTB$")
P.E4<-VlnPlot(E4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
N1[["percent.mt"]]<-PercentageFeatureSet(N1,pattern = "^ACOX1$|^ND1$|^ND2$|^COX1$|^COX2$|^ATP8$|^ATP6$|^COX3$|^ND3$|^ND4L$|^ND4$|^ND5$|^ND6$|^CYTB$")
P.N1<-VlnPlot(N1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
N2[["percent.mt"]]<-PercentageFeatureSet(N2,pattern = "^ACOX1$|^ND1$|^ND2$|^COX1$|^COX2$|^ATP8$|^ATP6$|^COX3$|^ND3$|^ND4L$|^ND4$|^ND5$|^ND6$|^CYTB$")
P.N2<-VlnPlot(N2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
N3[["percent.mt"]]<-PercentageFeatureSet(N3,pattern = "^ACOX1$|^ND1$|^ND2$|^COX1$|^COX2$|^ATP8$|^ATP6$|^COX3$|^ND3$|^ND4L$|^ND4$|^ND5$|^ND6$|^CYTB$")
P.N3<-VlnPlot(N3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
N4[["percent.mt"]]<-PercentageFeatureSet(N4,pattern = "^ACOX1$|^ND1$|^ND2$|^COX1$|^COX2$|^ATP8$|^ATP6$|^COX3$|^ND3$|^ND4L$|^ND4$|^ND5$|^ND6$|^CYTB$")
P.N4<-VlnPlot(N4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
N5[["percent.mt"]]<-PercentageFeatureSet(N5,pattern = "^ACOX1$|^ND1$|^ND2$|^COX1$|^COX2$|^ATP8$|^ATP6$|^COX3$|^ND3$|^ND4L$|^ND4$|^ND5$|^ND6$|^CYTB$")
P.N5<-VlnPlot(N5, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#QC
E1 <- subset(E1, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 10 & nCount_RNA >1000)
E2 <- subset(E2, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 10 & nCount_RNA >1000)
E3 <- subset(E3, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 10 & nCount_RNA >1000)
E4 <- subset(E4, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 10 & nCount_RNA >1000)
N1 <- subset(N1, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 10 & nCount_RNA >1000)
N2 <- subset(N2, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 10 & nCount_RNA >1000)
N3 <- subset(N3, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 10 & nCount_RNA >1000)
N4 <- subset(N4, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 10 & nCount_RNA >1000)
N5 <- subset(N5, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 10 & nCount_RNA >1000)

#Normalize
E1 <- NormalizeData(E1)
E2 <- NormalizeData(E2)
E3 <- NormalizeData(E3)
E4 <- NormalizeData(E4)
N1 <- NormalizeData(N1)
N2 <- NormalizeData(N2)
N3 <- NormalizeData(N3)
N4 <- NormalizeData(N4)
N5 <- NormalizeData(N5)

#Variable genes
E1 <- FindVariableFeatures(E1, selection.method = "vst", nfeatures = 2000)
E2 <- FindVariableFeatures(E2, selection.method = "vst", nfeatures = 2000)
E3 <- FindVariableFeatures(E3, selection.method = "vst", nfeatures = 2000)
E4 <- FindVariableFeatures(E4, selection.method = "vst", nfeatures = 2000)
N1 <- FindVariableFeatures(N1, selection.method = "vst", nfeatures = 2000)
N2 <- FindVariableFeatures(N2, selection.method = "vst", nfeatures = 2000)
N3 <- FindVariableFeatures(N3, selection.method = "vst", nfeatures = 2000)
N4 <- FindVariableFeatures(N4, selection.method = "vst", nfeatures = 2000)
N5 <- FindVariableFeatures(N5, selection.method = "vst", nfeatures = 2000)


## DoubletFinder ----
library(DoubletFinder)
# E1 0.0831----
E1 <- ScaleData(object = E1, verbose = FALSE)
E1 <- RunPCA(object = E1, npcs = 10, verbose = FALSE)
E1 <- FindNeighbors(object = E1, reduction = "pca", dims = 1:10)
E1 <- FindClusters(E1, resolution = 0.5)
E1 <- RunTSNE(object = E1, reduction = "pca", dims = 1:10)
sub.1 <- paramSweep_v3(E1, PCs = 1:10,sct = FALSE)
sub.2 <- summarizeSweep(sub.1, GT = FALSE)
sub.3 <- find.pK(sub.2)
mpK<-as.numeric(as.vector(sub.3$pK[which.max(sub.3$BCmetric)]))
pdf("./Doublet/pK_E1.pdf")
plot(sub.3$pK, sub.3$BCmetric, type='b')
dev.off()
annotations <- E1@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           
nExp_poi <- round(0.0831*nrow(E1@meta.data)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
seu_kidney <- doubletFinder_v3(E1, PCs = 1:10, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
seu_kidney <- doubletFinder_v3(seu_kidney, PCs = 1:10, pN = 0.25, pK = mpK, nExp = nExp_poi.adj, reuse.pANN = colnames(seu_kidney@meta.data)[9], sct = FALSE)
E1.Doubl <- seu_kidney

# E2 0.0426----
E2 <- ScaleData(object = E2, verbose = FALSE)
E2 <- RunPCA(object = E2, npcs = 10, verbose = FALSE)
E2 <- FindNeighbors(object = E2, reduction = "pca", dims = 1:10)
E2 <- FindClusters(E2, resolution = 0.5)
E2 <- RunTSNE(object = E2, reduction = "pca", dims = 1:10)
sub.1 <- paramSweep_v3(E2, PCs = 1:10,sct = FALSE)
sub.2 <- summarizeSweep(sub.1, GT = FALSE)
sub.3 <- find.pK(sub.2)
mpK<-as.numeric(as.vector(sub.3$pK[which.max(sub.3$BCmetric)]))
pdf("./Doublet/pK_E2.pdf")
plot(sub.3$pK, sub.3$BCmetric, type='b')
dev.off()
annotations <- E2@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           
nExp_poi <- round(0.0426*nrow(E2@meta.data)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
seu_kidney <- doubletFinder_v3(E2, PCs = 1:10, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
seu_kidney <- doubletFinder_v3(seu_kidney, PCs = 1:10, pN = 0.25, pK = mpK, nExp = nExp_poi.adj, reuse.pANN = colnames(seu_kidney@meta.data)[9], sct = FALSE)
E2.Doubl <- seu_kidney

# E3 0.0759----
E3 <- ScaleData(object = E3, verbose = FALSE)
E3 <- RunPCA(object = E3, npcs = 10, verbose = FALSE)
E3 <- FindNeighbors(object = E3, reduction = "pca", dims = 1:10)
E3 <- FindClusters(E3, resolution = 0.5)
E3 <- RunTSNE(object = E3, reduction = "pca", dims = 1:10)
sub.1 <- paramSweep_v3(E3, PCs = 1:10,sct = FALSE)
sub.2 <- summarizeSweep(sub.1, GT = FALSE)
sub.3 <- find.pK(sub.2)
mpK<-as.numeric(as.vector(sub.3$pK[which.max(sub.3$BCmetric)]))
pdf("./Doublet/pK_E3.pdf")
plot(sub.3$pK, sub.3$BCmetric, type='b')
dev.off()
annotations <- E3@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           
nExp_poi <- round(0.0759*nrow(E3@meta.data)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
seu_kidney <- doubletFinder_v3(E3, PCs = 1:10, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
seu_kidney <- doubletFinder_v3(seu_kidney, PCs = 1:10, pN = 0.25, pK = mpK, nExp = nExp_poi.adj, reuse.pANN = colnames(seu_kidney@meta.data)[9], sct = FALSE)
E3.Doubl <- seu_kidney

# E4 0.0797----
E4 <- ScaleData(object = E4, verbose = FALSE)
E4 <- RunPCA(object = E4, npcs = 10, verbose = FALSE)
E4 <- FindNeighbors(object = E4, reduction = "pca", dims = 1:10)
E4 <- FindClusters(E4, resolution = 0.5)
E4 <- RunTSNE(object = E4, reduction = "pca", dims = 1:10)
sub.1 <- paramSweep_v3(E4, PCs = 1:10,sct = FALSE)
sub.2 <- summarizeSweep(sub.1, GT = FALSE)
sub.3 <- find.pK(sub.2)
mpK<-as.numeric(as.vector(sub.3$pK[which.max(sub.3$BCmetric)]))
pdf("./Doublet/pK_E4.pdf")
plot(sub.3$pK, sub.3$BCmetric, type='b')
dev.off()
annotations <- E4@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           
nExp_poi <- round(0.0797*nrow(E4@meta.data)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
seu_kidney <- doubletFinder_v3(E4, PCs = 1:10, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
seu_kidney <- doubletFinder_v3(seu_kidney, PCs = 1:10, pN = 0.25, pK = mpK, nExp = nExp_poi.adj, reuse.pANN = colnames(seu_kidney@meta.data)[9], sct = FALSE)
E4.Doubl <- seu_kidney

# N1 0.065----
N1 <- ScaleData(object = N1, verbose = FALSE)
N1 <- RunPCA(object = N1, npcs = 10, verbose = FALSE)
N1 <- FindNeighbors(object = N1, reduction = "pca", dims = 1:10)
N1 <- FindClusters(N1, resolution = 0.5)
N1 <- RunTSNE(object = N1, reduction = "pca", dims = 1:10)
sub.1 <- paramSweep_v3(N1, PCs = 1:10,sct = FALSE)
sub.2 <- summarizeSweep(sub.1, GT = FALSE)
sub.3 <- find.pK(sub.2)
mpK<-as.numeric(as.vector(sub.3$pK[which.max(sub.3$BCmetric)]))
pdf("./Doublet/pK_N1.pdf")
plot(sub.3$pK, sub.3$BCmetric, type='b')
dev.off()
annotations <- N1@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           
nExp_poi <- round(0.065*nrow(N1@meta.data)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
seu_kidney <- doubletFinder_v3(N1, PCs = 1:10, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
seu_kidney <- doubletFinder_v3(seu_kidney, PCs = 1:10, pN = 0.25, pK = mpK, nExp = nExp_poi.adj, reuse.pANN = colnames(seu_kidney@meta.data)[9], sct = FALSE)
N1.Doubl <- seu_kidney

# N2 0.0455----
N2 <- ScaleData(object = N2, verbose = FALSE)
N2 <- RunPCA(object = N2, npcs = 10, verbose = FALSE)
N2 <- FindNeighbors(object = N2, reduction = "pca", dims = 1:10)
N2 <- FindClusters(N2, resolution = 0.5)
N2 <- RunTSNE(object = N2, reduction = "pca", dims = 1:10)
sub.1 <- paramSweep_v3(N2, PCs = 1:10,sct = FALSE)
sub.2 <- summarizeSweep(sub.1, GT = FALSE)
sub.3 <- find.pK(sub.2)
mpK<-as.numeric(as.vector(sub.3$pK[which.max(sub.3$BCmetric)]))
pdf("./Doublet/pK_N2.pdf")
plot(sub.3$pK, sub.3$BCmetric, type='b')
dev.off()
annotations <- N2@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           
nExp_poi <- round(0.0455*nrow(N2@meta.data)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
seu_kidney <- doubletFinder_v3(N2, PCs = 1:10, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
seu_kidney <- doubletFinder_v3(seu_kidney, PCs = 1:10, pN = 0.25, pK = mpK, nExp = nExp_poi.adj, reuse.pANN = colnames(seu_kidney@meta.data)[9], sct = FALSE)
N2.Doubl <- seu_kidney

# N3 0.0106----
N3 <- ScaleData(object = N3, verbose = FALSE)
N3 <- RunPCA(object = N3, npcs = 10, verbose = FALSE)
N3 <- FindNeighbors(object = N3, reduction = "pca", dims = 1:10)
N3 <- FindClusters(N3, resolution = 0.5)
N3 <- RunTSNE(object = N3, reduction = "pca", dims = 1:10)
sub.1 <- paramSweep_v3(N3, PCs = 1:10,sct = FALSE)
sub.2 <- summarizeSweep(sub.1, GT = FALSE)
sub.3 <- find.pK(sub.2)
mpK<-as.numeric(as.vector(sub.3$pK[which.max(sub.3$BCmetric)]))
pdf("./Doublet/pK_N3.pdf")
plot(sub.3$pK, sub.3$BCmetric, type='b')
dev.off()
annotations <- N3@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           
nExp_poi <- round(0.0106*nrow(N3@meta.data)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
seu_kidney <- doubletFinder_v3(N3, PCs = 1:10, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
seu_kidney <- doubletFinder_v3(seu_kidney, PCs = 1:10, pN = 0.25, pK = mpK, nExp = nExp_poi.adj, reuse.pANN = colnames(seu_kidney@meta.data)[9], sct = FALSE)
N3.Doubl <- seu_kidney

# N4 0.0685----
N4 <- ScaleData(object = N4, verbose = FALSE)
N4 <- RunPCA(object = N4, npcs = 10, verbose = FALSE)
N4 <- FindNeighbors(object = N4, reduction = "pca", dims = 1:10)
N4 <- FindClusters(N4, resolution = 0.5)
N4 <- RunTSNE(object = N4, reduction = "pca", dims = 1:10)
sub.1 <- paramSweep_v3(N4, PCs = 1:10,sct = FALSE)
sub.2 <- summarizeSweep(sub.1, GT = FALSE)
sub.3 <- find.pK(sub.2)
mpK<-as.numeric(as.vector(sub.3$pK[which.max(sub.3$BCmetric)]))
pdf("./Doublet/pK_N4.pdf")
plot(sub.3$pK, sub.3$BCmetric, type='b')
dev.off()
annotations <- N4@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           
nExp_poi <- round(0.0685*nrow(N4@meta.data)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
seu_kidney <- doubletFinder_v3(N4, PCs = 1:10, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
seu_kidney <- doubletFinder_v3(seu_kidney, PCs = 1:10, pN = 0.25, pK = mpK, nExp = nExp_poi.adj, reuse.pANN = colnames(seu_kidney@meta.data)[9], sct = FALSE)
N4.Doubl <- seu_kidney

# N5 0.0791----
N5 <- ScaleData(object = N5, verbose = FALSE)
N5 <- RunPCA(object = N5, npcs = 10, verbose = FALSE)
N5 <- FindNeighbors(object = N5, reduction = "pca", dims = 1:10)
N5 <- FindClusters(N5, resolution = 0.5)
N5 <- RunTSNE(object = N5, reduction = "pca", dims = 1:10)
sub.1 <- paramSweep_v3(N5, PCs = 1:10,sct = FALSE)
sub.2 <- summarizeSweep(sub.1, GT = FALSE)
sub.3 <- find.pK(sub.2)
mpK<-as.numeric(as.vector(sub.3$pK[which.max(sub.3$BCmetric)]))
pdf("./Doublet/pK_N5.pdf")
plot(sub.3$pK, sub.3$BCmetric, type='b')
dev.off()
annotations <- N5@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           
nExp_poi <- round(0.0791*nrow(N5@meta.data)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
seu_kidney <- doubletFinder_v3(N5, PCs = 1:10, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
seu_kidney <- doubletFinder_v3(seu_kidney, PCs = 1:10, pN = 0.25, pK = mpK, nExp = nExp_poi.adj, reuse.pANN = colnames(seu_kidney@meta.data)[9], sct = FALSE)
N5.Doubl <- seu_kidney

colnames(E1.Doubl@meta.data)[10] <- 'DF'
colnames(E2.Doubl@meta.data)[10] <- 'DF'
colnames(E3.Doubl@meta.data)[10] <- 'DF'
colnames(E4.Doubl@meta.data)[10] <- 'DF'
colnames(N1.Doubl@meta.data)[10] <- 'DF'
colnames(N2.Doubl@meta.data)[10] <- 'DF'
colnames(N3.Doubl@meta.data)[10] <- 'DF'
colnames(N4.Doubl@meta.data)[10] <- 'DF'
colnames(N5.Doubl@meta.data)[10] <- 'DF'
E1.final <- subset(E1.Doubl,subset = DF=='Singlet')
E2.final <- subset(E2.Doubl,subset = DF=='Singlet')
E3.final <- subset(E3.Doubl,subset = DF=='Singlet')
E4.final <- subset(E4.Doubl,subset = DF=='Singlet')
N1.final <- subset(N1.Doubl,subset = DF=='Singlet')
N2.final <- subset(N2.Doubl,subset = DF=='Singlet')
N3.final <- subset(N3.Doubl,subset = DF=='Singlet')
N4.final <- subset(N4.Doubl,subset = DF=='Singlet')
N5.final <- subset(N5.Doubl,subset = DF=='Singlet')

# integration ----
anchors=FindIntegrationAnchors(object.list = list(N1.final,N2.final,N3.final,N4.final,N5.final,
                                                  E1.final,E2.final,E3.final,E4.final),
                               dims = 1:50)
combined=IntegrateData(anchorset = anchors, dims = 1:50)
# DefaultAssay(object = combined) <- "RNA"
DefaultAssay(object = combined) <- "integrated"

combined <- ScaleData(combined,verbose = FALSE)
combined <- RunPCA(combined, verbose = FALSE,npcs = 50)
combined <- FindNeighbors(object = combined, dims = 1:50)
combined <- FindClusters(object = combined, resolution = 0.2)

# RunUMAP ----
sample.umap <- RunUMAP(object = combined, dims = 1:50)
DimPlot(object = sample.umap, reduction = 'umap',label = TRUE)

# FindAllMarkers ----
all.markers <- FindAllMarkers(sample.umap, logfc.threshold = 0.25)
# ENSMMUG00000040413 KRT13
# ENSMMUG00000047772 CD79A
# ENSMMUG00000063725 C1QA
