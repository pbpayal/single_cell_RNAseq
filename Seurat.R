
install.packages('tidyverse')
library(tidyverse)  
install.packages('metap')
library('metap')
library(patchwork)
install.packages('Seurat')
library(Seurat)
install.packages("dplyr")
install.packages("magrittr") # package installations are only needed the first time you use it
install.packages("dplyr")    # alternative installation of the %>%
library(magrittr) # needs to be run every time you start R and want to use %>%
library(dplyr) 
library(ggplot2)
library(RColorBrewer)

#########################
######## Control ########
#########################
sample_300 <- Read10X(data.dir = "/Documents/2021/RUN2/300/filtered_feature_bc_matrix")
sample_300_sobj <- CreateSeuratObject(counts = sample_300, project = "sample_300", min.cells = 3, min.features = 200)
dim(sample_300_sobj)
# [1] 19410  5976
sample_300_sobj$stim <- "CTRL"
sample_300_sobj[["percent.mt"]] <- PercentageFeatureSet(sample_300_sobj, pattern = "^MT-")
sample_300_sub <- subset(sample_300_sobj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 )
sample_300_sub <- NormalizeData(sample_300_sub, verbose = FALSE)
sample_300_sub <- FindVariableFeatures(sample_300_sub, selection.method = "vst", nfeatures = 2000)
dim(sample_300_sobj)
# store mitochondrial percentage in object meta data
#sample_300_sub <- PercentageFeatureSet(sample_300_sobj, pattern = "^MT-", col.name = "percent.mt")
# run sctransform
#sample_300_sub <- SCTransform(sample_300_sub, vars.to.regress = "percent.mt", verbose = FALSE)
# ### Mitochondrial genes
# grep("^MT-",rownames(sample_300_sobj@assays$RNA@counts),value = TRUE)
# PercentageFeatureSet(sample_300_sobj,pattern="^MT-") -> sample_300_sobj$percent.MT
# head(sample_300_sobj$percent.MT)
# ##### Ribosomal Genes ####
# grep("^RP[LS]",rownames(sample_300_sobj@assays$RNA@counts),value = TRUE)
# PercentageFeatureSet(sample_300_sobj,pattern="^RP[LS]") -> sample_300_sobj$percent.Ribosomal
# head(sample_300_sobj$percent.Ribosomal)
# DotPlot(sample_300_sobj, features=c("nCount_RNA","percent.MT", "percent.Ribosomal"))
# [1] 19410  5976
sample_301 <- Read10X(data.dir = "/Documents/2021/RUN2/301/filtered_feature_bc_matrix")
sample_301_sobj <- CreateSeuratObject(counts = sample_301, project = "sample_301", min.cells = 3, min.features = 200)
dim(sample_301_sobj)
# [1] 17775  5967
sample_301_sobj$stim <- "CTRL"
sample_301_sobj[["percent.mt"]] <- PercentageFeatureSet(sample_301_sobj, pattern = "^MT-")
sample_301_sub <- subset(sample_301_sobj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5 )
sample_301_sub <- NormalizeData(sample_301_sub, verbose = FALSE)
sample_301_sub <- FindVariableFeatures(sample_301_sub, selection.method = "vst", nfeatures = 2000)
# # store mitochondrial percentage in object meta data
# sample_301_sub <- PercentageFeatureSet(sample_301_sobj, pattern = "^MT-", col.name = "percent.mt")
# # run sctransform
# # sample_301_sub <- SCTransform(sample_301_sub, vars.to.regress = "percent.mt", verbose = FALSE)
# 
sample_165 <- Read10X(data.dir = "/Documents/2021/RUN4/165/filtered_feature_bc_matrix")
sample_165_sobj <- CreateSeuratObject(counts = sample_165, project = "sample_165", min.cells = 3, min.features = 200)
dim(sample_165_sobj)
# # [1] 20737 24694
sample_165_sobj$stim <- "CTRL"
sample_165_sobj[["percent.mt"]] <- PercentageFeatureSet(sample_165_sobj, pattern = "^MT-")
sample_165_sub <- subset(sample_165_sobj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 )
sample_165_sub <- NormalizeData(sample_165_sub, verbose = FALSE)
sample_165_sub <- FindVariableFeatures(sample_165_sub, selection.method = "vst", nfeatures = 2000)
# # store mitochondrial percentage in object meta data
# sample_165_sub <- PercentageFeatureSet(sample_165_sobj, pattern = "^MT-", col.name = "percent.mt")
# # run sctransform
# sample_165_sub <- SCTransform(sample_165_sub, vars.to.regress = "percent.mt", verbose = FALSE)

sample_295 <- Read10X(data.dir = "/Documents/2021/RUN4/295/filtered_feature_bc_matrix")
sample_295_sobj <- CreateSeuratObject(counts = sample_295, project = "sample_295", min.cells = 3, min.features = 200)
dim(sample_295_sobj)
# [1] 21600 25303
sample_295_sobj$stim <- "CTRL"
sample_295_sobj[["percent.mt"]] <- PercentageFeatureSet(sample_295_sobj, pattern = "^MT-")
sample_295_sub <- subset(sample_295_sobj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)
sample_295_sub <- NormalizeData(sample_295_sub, verbose = FALSE)
sample_295_sub <- FindVariableFeatures(sample_295_sub, selection.method = "vst", nfeatures = 2000)
# # store mitochondrial percentage in object meta data
# sample_295_sub <- PercentageFeatureSet(sample_295_sobj, pattern = "^MT-", col.name = "percent.mt")
# # # run sctransform
# # sample_295_sub <- SCTransform(sample_295_sub, vars.to.regress = "percent.mt", verbose = FALSE)


sample_259 <- Read10X(data.dir = "/Documents/2021/259/filtered_feature_bc_matrix/")
sample_259_sobj <- CreateSeuratObject(counts = sample_259, project = "sample_259", min.cells = 3, min.features = 200)
dim(sample_259_sobj)
# [1] 21600 25303
sample_259_sobj$stim <- "CTRL"
sample_259_sobj[["percent.mt"]] <- PercentageFeatureSet(sample_259_sobj, pattern = "^MT-")
sample_259_sobj <- subset(sample_259_sobj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)
sample_259_sobj <- NormalizeData(sample_259_sobj, verbose = FALSE)
sample_259_sobj <- FindVariableFeatures(sample_259_sobj, selection.method = "vst", nfeatures = 2000)

sample_260 <- Read10X(data.dir = "/Documents/2021/260/filtered_feature_bc_matrix/")
sample_260_sobj <- CreateSeuratObject(counts = sample_260, project = "sample_260", min.cells = 3, min.features = 200)
dim(sample_260_sobj)
# [1] 21600 25303
sample_260_sobj$stim <- "CTRL"
sample_260_sobj[["percent.mt"]] <- PercentageFeatureSet(sample_260_sobj, pattern = "^MT-")
sample_260_sobj <- subset(sample_260_sobj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)
sample_260_sobj <- NormalizeData(sample_260_sobj, verbose = FALSE)
sample_260_sobj <- FindVariableFeatures(sample_260_sobj, selection.method = "vst", nfeatures = 2000)
# store mitochondrial percentage in object meta data
sample_260_sobj <- PercentageFeatureSet(sample_260_sobj, pattern = "^MT-", col.name = "percent.mt")


#########################
########## MD ###########
#########################
sample_326 <- Read10X(data.dir = "/Documents/2021/RUN1/326/filtered_feature_bc_matrix")
sample_326_sobj <- CreateSeuratObject(counts = sample_326, project = "sample_326", min.cells = 3, min.features = 200)
dim(sample_326_sobj)
# [1] 18392  4782
sample_326_sobj$stim <- "MD"
sample_326_sobj[["percent.mt"]] <- PercentageFeatureSet(sample_326_sobj, pattern = "^MT-")
sample_326_sub <- subset(sample_326_sobj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 )
sample_326_sub <- NormalizeData(sample_326_sub, verbose = FALSE)
sample_326_sub <- FindVariableFeatures(sample_326_sub, selection.method = "vst", nfeatures = 2000)
# store mitochondrial percentage in object meta data
# sample_326_sub <- PercentageFeatureSet(sample_326_sobj, pattern = "^MT-", col.name = "percent.mt")
# # run sctransform
# sample_326_sub <- SCTransform(sample_326_sub, vars.to.regress = "percent.mt", verbose = FALSE)

sample_319 <- Read10X(data.dir = "/Documents/2021/RUN1/319/filtered_feature_bc_matrix")
sample_319_sobj <- CreateSeuratObject(counts = sample_319, project = "sample_319", min.cells = 3, min.features = 200)
dim(sample_319_sobj)
# [1] 17644  3375
sample_319_sobj$stim <- "MD"
sample_319_sobj[["percent.mt"]] <- PercentageFeatureSet(sample_319_sobj, pattern = "^MT-")
sample_319_sub <- subset(sample_319_sobj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 )
sample_319_sub <- NormalizeData(sample_319_sub, verbose = FALSE)
sample_319_sub <- FindVariableFeatures(sample_319_sub, selection.method = "vst", nfeatures = 2000)
# store mitochondrial percentage in object meta data
# sample_319_sub <- PercentageFeatureSet(sample_319_sobj, pattern = "^MT-", col.name = "percent.mt")
# # run sctransform
# sample_319_sub <- SCTransform(sample_319_sub, vars.to.regress = "percent.mt", verbose = FALSE)

sample_155 <- Read10X(data.dir = "/Documents/2021/RUN4/155/filtered_feature_bc_matrix")
sample_155_sobj <- CreateSeuratObject(counts = sample_155, project = "sample_155", min.cells = 3, min.features = 200)
dim(sample_155_sobj)
# [1] 19114  9793
sample_155_sobj$stim <- "MD"
sample_155_sobj[["percent.mt"]] <- PercentageFeatureSet(sample_155_sobj, pattern = "^MT-")
sample_155_sub <- subset(sample_155_sobj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5 )
sample_155_sub <- NormalizeData(sample_155_sub, verbose = FALSE)
sample_155_sub <- FindVariableFeatures(sample_155_sub, selection.method = "vst", nfeatures = 2000)
# # store mitochondrial percentage in object meta data
# sample_155_sub <- PercentageFeatureSet(sample_155_sobj, pattern = "^MT-", col.name = "percent.mt")
# # run sctransform
# sample_155_sub <- SCTransform(sample_155_sub, vars.to.regress = "percent.mt", verbose = FALSE)

sample_132 <- Read10X(data.dir = "/Documents/2021/RUN4/132/filtered_feature_bc_matrix")
sample_132_sobj <- CreateSeuratObject(counts = sample_132, project = "sample_132", min.cells = 3, min.features = 200)
dim(sample_132_sobj)
# [1] 21405 12925
sample_132_sobj$stim <- "MD"
sample_132_sobj[["percent.mt"]] <- PercentageFeatureSet(sample_132_sobj, pattern = "^MT-")
sample_132_sub <- subset(sample_132_sobj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5 )
sample_132_sub <- NormalizeData(sample_132_sub, verbose = FALSE)
sample_132_sub <- FindVariableFeatures(sample_132_sub, selection.method = "vst", nfeatures = 2000)
# # store mitochondrial percentage in object meta data
# sample_132_sub <- PercentageFeatureSet(sample_132_sobj, pattern = "^MT-", col.name = "percent.mt")
# # run sctransform
# sample_132_sub <- SCTransform(sample_132_sub, vars.to.regress = "percent.mt", verbose = FALSE)


sample_292 <- Read10X(data.dir = "/Documents/2021/RUN2/292/filtered_feature_bc_matrix/")
sample_292_sobj <- CreateSeuratObject(counts = sample_292, project = "sample_292", min.cells = 3, min.features = 200)
dim(sample_132_sobj)
# [1] 21405 12925
sample_132_sobj$stim <- "MD"
sample_132_sobj[["percent.mt"]] <- PercentageFeatureSet(sample_132_sobj, pattern = "^MT-")
sample_132_sub <- subset(sample_132_sobj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5 )
sample_132_sub <- NormalizeData(sample_132_sub, verbose = FALSE)
sample_132_sub <- FindVariableFeatures(sample_132_sub, selection.method = "vst", nfeatures = 2000)


# Run Single cell steps
immune.anchors.2 <- FindIntegrationAnchors(object.list = list(sample_300_sub, sample_301_sub,sample_165_sub, sample_295_sub,
                                                            sample_326_sub,sample_319_sub,sample_155_sub,sample_132_sub), dims = 1:20)
immune.combined.integrated <- IntegrateData(anchorset = immune.anchors.2, dims = 1:20)
dim(immune.combined.integrated)
saveRDS(immune.combined.integrated, file = "/Documents/2021/Part2/immune.combined_integrated_part2.rds")
immune.combined.integrated_1 <- readRDS("/Documents/2021/Part2/immune.combined_integrated_part2.rds")
# [1]  2000 31584
table(immune.combined.integrated_1$orig.ident)
# sample_132 sample_155 sample_165 sample_295 sample_300 sample_301 sample_319 sample_326 
# 3004       3576       7042       9932       1572       3841        827       1790 
table(immune.combined.integrated_1@meta.data$stim)
# CTRL    MD 
# 22387  9197
## Table of cells per sample


no_of_cells_per_sample <- read.csv("/Documents/2021/Part2/no_of_cells_per_sample.csv")


# Simple Bar Plot
coul <- brewer.pal(8,"Set1", c(4)) 
colors <- rep(c("purple","green"),each=4)
#counts_per_sample <- table(immune.combined.integrated$orig.ident)
barplot(as.numeric(no_of_cells_per_sample$Log10_no_of_cells),main="Number of cells per sample",decreasing=TRUE, 
        horiz = FALSE,
        xlab = "Samples",
        ylab="log10(Cell number)", 
        names.arg=c("Control-1", "Control-2", "Control-3", "Control-4", 
                                               "Leigh-1","Leigh-2", "Leigh-3", "Leigh-4"), 
        las=0, cex.names=0.6, cex.axis=0.7,
        col=colors)
        

       
print(head(immune.combined.integrated@meta.data))
print(tail(immune.combined.integrated@meta.data))

### Ribosomal genes
DefaultAssay(immune.combined.integrated) <- "RNA"
grep("^RP[LS]",rownames(immune.combined.integrated@assays$RNA@counts),value = TRUE)
PercentageFeatureSet(immune.combined.integrated,pattern="^RP[LS]") -> immune.combined.integrated$percent.Ribosomal
head(immune.combined.integrated$percent.Ribosomal)
### Mitochondrial genes
grep("^MT-",rownames(immune.combined.integrated@assays$RNA@counts),value = TRUE)
PercentageFeatureSet(immune.combined.integrated,pattern="^MT-") -> immune.combined.integrated$percent.MT
head(immune.combined.integrated$percent.MT)

DefaultAssay(immune.combined.integrated_1) <- "integrated"

# Run the standard workflow for visualization and clustering
immune.combined.integrated_1 <- ScaleData(immune.combined.integrated_1, verbose = FALSE)
immune.combined.integrated_1 <- RunPCA(immune.combined.integrated_1, npcs = 20, verbose = FALSE)

immune.combined.integrated_jack <- JackStraw(immune.combined.integrated_1, num.replicate = 100)
immune.combined.integrated_jack <- ScoreJackStraw(immune.combined.integrated_jack, dims = 1:20)
JackStrawPlot(immune.combined.integrated_jack, dims = 1:15)
ElbowPlot(immune.combined.integrated_jack)
#save_plot(paste0(output_directory,outputPrefix, "_Elbowplot.png"),ElbowPlot(md_sub_pca_jack), base_height=8,base_width=7)

# t-SNE and Clustering
immune.combined.integrated_1 <- RunUMAP(immune.combined.integrated_1, reduction = "pca", dims = 1:20)
immune.combined.integrated_1 <- FindNeighbors(immune.combined.integrated_1, reduction = "pca", dims = 1:20)
# Resolution
immune.combined.integrated_1 <- FindClusters(immune.combined.integrated_1, resolution = 0.2)
saveRDS(immune.combined.integrated_1, file = "/Documents/2021/Part2/immune.combined_integrated_part2_res0.2_unannotated.rds")
immune.combined.integrated_1 <- readRDS("/Documents/2021/Part2/immune.combined_integrated_part2_res0.2_unannotated.rds")
DefaultAssay(immune.combined.integrated_1)
# Visualization
p1 <- DimPlot(immune.combined.integrated_1, reduction = "umap", group.by = "stim")
p1
p2 <- DimPlot(immune.combined.integrated_1, reduction = "umap", label = TRUE)
p2
p1 + p2


DimPlot(immune.combined.integrated_1, reduction = "umap", split.by = "stim")


DefaultAssay(immune.combined.integrated_1) <- "RNA"
comb_cluster_markers_2 <- FindAllMarkers(immune.combined.integrated_1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
output_directory <- "/Documents/2021/Part2/"
outputPrefix <- "rerun"
write.csv(comb_cluster_markers_2,file = paste0(output_directory,outputPrefix, "part2_res0.4_markers.csv"))

comb_cluster_markers_2 %>%
  group_by(cluster) %>% 
  slice_max(n=5, order_by = avg_log2FC) -> top5_cluster_markers_2
top5_cluster_markers_2
outputPrefix <- "Part2_"
outputPrefix <- "rerun_res0.4"
output_directory <- "/Documents/2021/Part2/"
write.csv(top5_cluster_markers_2,file = paste0(output_directory,outputPrefix, "top5_res0.2_markers.csv"))

comb_cluster_markers_2 %>%
  group_by(cluster) %>% 
  slice_max(n=10, order_by = avg_log2FC) -> top10_cluster_markers_2
top10_cluster_markers_2
write.csv(top10_cluster_markers_2,file = paste0(output_directory,outputPrefix, "_top10_res0.4_markers.csv"))

# Annotate clusters
new.cluster.ids <- c("Erythroid cells", "Memory Tcells/Tcells", "DC/Monocyte",
                     "NK" ,"Bcells", "CD8 Tcells","MAIT T cell",
                    "Plasma","Plasmacytoid DC", "RBC","Platelets")
names(new.cluster.ids) <- levels(immune.combined.integrated_1)
immune.combined.integrated_1 <- RenameIdents(immune.combined.integrated_1, new.cluster.ids)
# Rename the seurat clusters in metadata
immune.combined.integrated_1[["seurat_clusters"]] <- Idents(object = immune.combined.integrated_1)
# # Rename classes.
# stem.combined <- RenameIdents(object = stem.combined, `0` = "your cell type", `1` = "your other cell type", `2` = "your last cell type")

DimPlot(immune.combined.integrated_1, reduction = "umap", label = TRUE, 
        pt.size = 0.6, label.size = 2.5) 
saveRDS(immune.combined.integrated_1, file = "/Documents/2021/Part2/immune.combined_integrated_part2_res0.2_annotated.rds")
immune.combined.integrated_1 <- readRDS("/Documents/2021/Part2/immune.combined_integrated_part2_res0.2_annotated.rds")

immune.combined.integrated_1 <- readRDS("/Documents/Payal_scripts/2021_Singlecell_Peter/immune.combined_integrated_part2_res0.2_annotated.rds")

# Visualization
p1 <- DimPlot(immune.combined.integrated_1, reduction = "umap", group.by = "stim")
p1
p2 <- DimPlot(immune.combined.integrated_1, reduction = "umap", label = TRUE)
p2
p1 + p2

setwd("/Documents/2021/2021_Peter_SC/")
png("clusters_res0.2_highres.png", width = 30, height = 14, units = 'in', res = 300)
print(DimPlot(immune.combined.integrated_1, reduction = "umap", label = TRUE, label.size = 10, pt.size = 0.6) + theme(text = element_text(face = "bold")))
# Make plot
p2
#LabelClusters(p2, id = "ident",  fontface = "bold", color = "black")
dev.off()

features = c("MS4A1","IL7R","LYZ","GNLY",
             "CD8B","KLRB1","MZB1","GZMB", 
             "HBB","PPBP" )


DimPlot(immune.combined.integrated_1, reduction = "umap", split.by = "stim")

FeaturePlot(immune.combined.integrated_1, features = c("MS4A1", "IGKC", "HLA-DRA"))
FeaturePlot(immune.combined.integrated_1, features = c("GZMB", "GNLY"))

Idents(immune.combined.integrated_1)
DefaultAssay(immune.combined.integrated_1) <- "RNA"
#DefaultAssay(immune.combined.integrated_1) <- "integrated"
features1 = c("MT-ND1","MT-ND2","MT-ND3","MT-ND4","MT-ND4L","MT-ND5","MT-ND","MT-CYB",
             "MT-CO1","MT-CO2","MT-CO3","MT-ATP6","MT-ATP8","MT-RNR1","MT-TA", "MT-TC","MT-TD","MT-TE")
#The following requested variables were not found: MT-ND, MT-RNR1, MT-TA, MT-TC, MT-TD, MT-TE
features2 =c("MT-TF","MT-TG","MT-TH","MT-TI","MT-TK","MT-TL1","MT-TL2",
             "MT-TM","MT-TN","MT-TP","MT-TQ","MT-TR","MT-TS1","MT-TS2",
             "MT-TT","MT-TV","MT-TW","MT-TY")
#The following requested variables were not found: MT-TF, MT-TG, MT-TH, MT-TI, 
# MT-TK, MT-TL1, MT-TL2, MT-TM, MT-TN, MT-TP,MT-TQ, MT-TR, MT-TS1, MT-TS2, MT-TT, MT-TV, MT-TW, MT-TY
DotPlot(immune.combined.integrated_1, features=features1)
DotPlot(immune.combined.integrated_1, features=features2)
FeaturePlot(immune.combined.integrated_1, features=features1)
#Nuc
# AD
features3 =c("SDHC","SDHAF2","CYCS","RTN4IP1","DNA2","CHCHD10","DLST","PPOX","IDH2")
# CYCS was duplicate
DotPlot(immune.combined.integrated_1, features=features3)
FeaturePlot(immune.combined.integrated_1, features=features3)
#AR
features4 = c("NDUFA2","NDUFA6","NDUFA9","NDUFA10","NDUFA11","NDUFA12",
              "NDUFA13","NDUFB3","NDUFB8","NDUFB9","NDUFB10","NDUFS1","NDUFS2",
              "NDUFS3","NDUFS4","NDUFS6","NDUFS7","NDUFS8","NDUFV1","NDUFV2")

DotPlot(immune.combined.integrated_1, features=features4)
FeaturePlot(immune.combined.integrated_1, features=features4)
features5 = c("ACAD9","FOXRED1","NDUFAF1","NDUFAF2","NDUFAF3","NDUFAF4","NDUFAF5","NDUFAF6","NDUFAF8",
              "NUBPL","TIMMDC1","TMEM126B","SDHAF1","CYC1","UQCRB","UQCRC2","UQCRFS1","UQCRQ","BCS1L")
DotPlot(immune.combined.integrated_1, features=features5)
FeaturePlot(immune.combined.integrated_1, features=features5)
features6 =c("LYRM7","TTC19","UQCC2","UQCC3","COQ2","COQ4",
             "COQ5","COQ6","COQ7","COQ8A","COQ8B","COQ9","PDSS1",
             "PDSS2","COX4I1","COX4I2","COX5A","COX6A1","COX6B1",
             "COX8A","NDUFA4","CEP89")
DotPlot(immune.combined.integrated_1, features=features6)
FeaturePlot(immune.combined.integrated_1, features=features6)
features7 =c("COA3","COA5","COA6","COA7","COX10","COX14","COX15","COX20",
             "FASTKD2","OXA1L","PET100","PET117","SCO1","SCO2","SURF1","TACO1")
DotPlot(immune.combined.integrated_1, features=features7)
FeaturePlot(immune.combined.integrated_1, features=features7)

features8 =c("ATP5F1A","ATP5F1D","ATP5F1E","ATPAF2","ATP5MD","TMEM70",
             "FBXL4","MGME1","MPV17","RNASEH1","TMEM65","ABAT","DGUOK","SAMHD1",
             "SUCLA2","SUCLG1","TK2","TYMP","ELAC2","ERAL1")
DotPlot(immune.combined.integrated_1, features=features8)
FeaturePlot(immune.combined.integrated_1, features=features8)

features9 =c("GTPBP3","LRPPRC","MRM2","MTO1","MTPAP","NSUN3","PNPT1","PUS1",
             "TRIT1","TRMT5","TRMT10C","TRMU","TRNT1","C12orf65","C1QBP",
             "GFM1","GFM2","GUF1","RMND1","TSFM","TUFM")
DotPlot(immune.combined.integrated_1, features=features9)
FeaturePlot(immune.combined.integrated_1, features=features9)


features10 =c("MRPL3","MRPL12","MRPL44","MRPS2","MRPS7","MRPS14","MRPS16",
              "MRPS22","MRPS23","MRPS25","MRPS28","MRPS34",
              "PTCD3","AARS2","CARS2","DARS2","EARS2","FARS2",
              "GARS","GATB","GATC","HARS2","IARS2")
DotPlot(immune.combined.integrated_1, features=features10)
FeaturePlot(immune.combined.integrated_1, features=features10)

features11 =c("KARS","LARS2","MARS2","MTFMT","NARS2","PARS2","QRSL1","RARS2",
              "SARS2","TARS2","VARS2","WARS2","YARS2","NME3","MFF","STAT2","SLC25A46",
              "C19orf70","MIEF2","YME1L1","APOPT1","MICU1","MICU2","CLPB","CLPP",
              "HTRA2","LONP1","MIPEP")
DotPlot(immune.combined.integrated_1, features=features11)
FeaturePlot(immune.combined.integrated_1, features=features11)

features12 =c("PITRM1","PMPCB","SACS","AGK","DNAJC19","GFER","PAM16","PISD","PMPCA",
              "SERAC1","TIMM22","TIMM50","TOMM70","XPNPEP3","ETFA","ETFB","ETFDH","HADHA",
              "HADHB","SLC22A5","SLC25A20","ACAT1")
DotPlot(immune.combined.integrated_1, features=features12)
FeaturePlot(immune.combined.integrated_1, features=features12)

features13 =c("HMGCL","HMGCS2","OXCT1","CA5A","PC","ACO2","FH","IDH3A","IDH3B","MDH2",
              "GOT2","MDH1","SLC25A13","DLAT","DLD","PDHB","PDHX","PDP1","MPC1","SLC25A3",
              "SLC25A10","SLC25A12","SLC25A26","BOLA3","FDX1L","FDXR","FXN")
DotPlot(immune.combined.integrated_1, features=features13)
#The following requested variables were not found: HMGCS2, CA5A, FDX1L
FeaturePlot(immune.combined.integrated_1, features=features13)

features14 =c("GLRX5","IBA57","ISCA1","ISCA2","ISCU","LYRM4",
              "NFS1","NFU1","BTD","HLCS","SLC19A3","SLC25A19",
              "TPK1","LIAS","LIPT1","LIPT2","MCAT")
DotPlot(immune.combined.integrated_1, features=features14)
#The following requested variables were not found: HMGCS2, CA5A, FDX1L
FeaturePlot(immune.combined.integrated_1, features=features14)





# Make the dataet ready by adding custom column to accomodate celltype and experimnet information
immune.combined.integrated_1$celltype.exp <- paste(Idents(immune.combined.integrated_1), immune.combined.integrated_1$stim, sep = "_")
immune.combined.integrated_1$celltype <- Idents(immune.combined.integrated_1)
Idents(immune.combined.integrated_1) <- "celltype.exp"
head(Idents(immune.combined.integrated_1))
tail(Idents(immune.combined.integrated_1))

# DEG
output_directory <- "/Documents/2021/Part2/"
outputPrefix <- "rerun_int"
DefaultAssay(immune.combined.integrated_1) <- "RNA"
DefaultAssay(immune.combined.integrated_1) <- "integrated"
# ALL DEGs
immune.deg <- FindMarkers(immune.combined.integrated_1, ident.1 = "CTRL", ident.2 = "MD", verbose = TRUE, group.by="stim", logfc.threshold = log(2))
head(immune.deg)
write.csv(immune.deg,file = paste0(output_directory,outputPrefix, "all_DEGs.csv"))
# DEGs per cluster
immune.deg.cluster0 <- FindMarkers(immune.combined.integrated_1, ident.1 = "CTRL", ident.2 = "MD", verbose = TRUE, group.by="stim", subset.ident = "Erythroid cells", logfc.threshold = log(2))
write.csv(immune.deg.cluster0,file = paste0(output_directory,outputPrefix, "Erythroid_cells_cluster0_DEGs.csv"))
immune.deg.cluster1 <- FindMarkers(immune.combined.integrated_1, ident.1 = "CTRL", ident.2 = "MD", verbose = TRUE, group.by="stim", subset.ident = "Memory Tcells/Tcells", logfc.threshold = log(2))
write.csv(immune.deg.cluster1,file = paste0(output_directory,outputPrefix, "Memory_Tcells_and_Tcells_cluster1_DEGs.csv"))
immune.deg.cluster2 <- FindMarkers(immune.combined.integrated_1, ident.1 = "CTRL", ident.2 = "MD", verbose = TRUE, group.by="stim", subset.ident = "DC/Monocyte", logfc.threshold = log(2))
write.csv(immune.deg.cluster2,file = paste0(output_directory,outputPrefix, "DC_Monocyte_cluster2_DEGs.csv"))
immune.deg.cluster3 <- FindMarkers(immune.combined.integrated_1, ident.1 = "CTRL", ident.2 = "MD", verbose = TRUE, group.by="stim", subset.ident = "NK", logfc.threshold = log(2))
write.csv(immune.deg.cluster3,file = paste0(output_directory,outputPrefix, "NK_cluster3_DEGs.csv"))
immune.deg.cluster4 <- FindMarkers(immune.combined.integrated_1, ident.1 = "CTRL", ident.2 = "MD", verbose = TRUE, group.by="stim", subset.ident = "Bcells", logfc.threshold = log(2))
write.csv(immune.deg.cluster4,file = paste0(output_directory,outputPrefix, "Bcells_cluster4_DEGs.csv"))
immune.deg.cluster5 <- FindMarkers(immune.combined.integrated_1, ident.1 = "CTRL", ident.2 = "MD", verbose = TRUE, group.by="stim", subset.ident = "CD8 Tcells", logfc.threshold = log(2))
write.csv(immune.deg.cluster5,file = paste0(output_directory,outputPrefix, "CD8_Tcells_cluster5_DEGs.csv"))
immune.deg.cluster6 <- FindMarkers(immune.combined.integrated_1, ident.1 = "CTRL", ident.2 = "MD", verbose = TRUE, group.by="stim", subset.ident = "MAIT T cell", logfc.threshold = log(2))
write.csv(immune.deg.cluster6,file = paste0(output_directory,outputPrefix, "MAIT_T_cell_cluster6_DEGs.csv"))
immune.deg.cluster7 <- FindMarkers(immune.combined.integrated_1, ident.1 = "CTRL", ident.2 = "MD", verbose = TRUE, group.by="stim", subset.ident = "Plasma", logfc.threshold = log(2))
write.csv(immune.deg.cluster7,file = paste0(output_directory,outputPrefix, "Plasma_cluster7_DEGs.csv"))
immune.deg.cluster8 <- FindMarkers(immune.combined.integrated_1, ident.1 = "CTRL", ident.2 = "MD", verbose = TRUE, group.by="stim", subset.ident = "Plasmacytoid DC", logfc.threshold = log(2))
write.csv(immune.deg.cluster8,file = paste0(output_directory,outputPrefix, "Plasmacytoid DC_cluster8_DEGs.csv"))
immune.deg.cluster9 <- FindMarkers(immune.combined.integrated_1, ident.1 = "CTRL", ident.2 = "MD", verbose = TRUE, group.by="stim", subset.ident = "RBC", logfc.threshold = log(2))
write.csv(immune.deg.cluster9,file = paste0(output_directory,outputPrefix, "RBC_cluster9_DEGs.csv"))
immune.deg.cluster10 <- FindMarkers(immune.combined.integrated_1, ident.1 = "CTRL", ident.2 = "MD", verbose = TRUE, group.by="stim", subset.ident = "Platelets", logfc.threshold = log(2))
write.csv(immune.deg.cluster10,file = paste0(output_directory,outputPrefix, "Platelets_cluster10_DEGs.csv"))


FeaturePlot(bcells, features = c("MTRNR2L8","PLCG2", "HLA-DQA2", "MS4A1", "IGKC"), split.by = "stim")


# Subset on the expression level of a gene/feature
bcells <- subset(immune.combined.integrated_1, idents = "Bcells")
DimPlot(bcells, reduction = "umap", label = TRUE)
immune.deg.cluster4 <- FindMarkers(bcells, ident.1 = "CTRL", ident.2 = "MD", verbose = TRUE, group.by="stim", logfc.threshold = log(2))
immune.deg.cluster4
avg.bcells <- log1p(AverageExpression(bcells, verbose = FALSE)$RNA)
avg.bcells
avg.bcells$gene <- rownames(avg.bcells)
p2 <- ggplot(avg.bcells, aes(CTRL, MD)) + geom_point() + ggtitle("Bcells")

# Marker Dot plot
features = c("MS4A1","IL7R","LYZ","GNLY",
             "CD8B","KLRB1","MZB1","GZMB", 
             "HBB","PPBP" )
DotPlot(immune.combined.integrated, features = features) + RotatedAxis()
DotPlot(immune.combined.integrated, features = "MZB1") + RotatedAxis()
# Feature Plot
FeaturePlot(immune.combined.integrated,features = features )
FeaturePlot(immune.combined.integrated,features = "MZB1" )
# How many cells per group per condition diagram
x <- table(Idents(immune.combined.integrated), immune.combined.integrated$stim)
# CTRL   MD
# 0  8287 1582
# 1  3044 2804
# 2  3441  776
# 3  2814 1101
# 4  2060  647
# 5  1110 1222
# 6   735  734
# 7   721   55
# 8    82  144
# 9    28   89
# 10   65   43
# CTRL   MD
# Erythroid cells      8287 1582
# Memory Tcells/Tcells 3044 2804
# DC/Monocyte          3441  776
# NK                   2814 1101
# Bcells               2060  647
# CD8 Tcells           1110 1222
# MAIT T cell           735  734
# Plasma                721   55
# Plasmacytoid DC        82  144
# RBC                    28   89
# Platelets              65   43
x <- as.data.frame(x)
x$Var1 <- as.character(x$Var1)

ggplot(x, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  xlab("Sample") +
  ylab("Proportion") +
  scale_fill_manual(values = brewer.pal(12, "Paired")) +
  theme(legend.title = element_blank())

# How many cells are in each replicate?
table(immune.combined.integrated$orig.ident)
# sample_132 sample_155 sample_165 sample_295 sample_300 sample_301 sample_319 sample_326 
# 3004       3576       7042       9932       1572       3841        827       1790 

# What proportion of cells are in each cluster?
prop.table(table(Idents(immune.combined.integrated_1)))

# How does cluster membership vary by replicate?
table(Idents(immune.combined.integrated_1), immune.combined.integrated_1$orig.ident)
# sample_132 sample_155 sample_165 sample_295 sample_300 sample_301 sample_319 sample_326
# 0         375        641       1791       5526        154        816        275        291
# 1         588       1403        874       1130        670        370        263        550
# 2         445        130        994       1476         67        904         11        190
# 3         604        262       1285        546        182        801         29        206
# 4         212        261       1358        329        105        268         17        157
# 5         360        550        272        478        185        175        158        154
# 6         191        291        233        123        203        176         65        187
# 7          29          2        181        267          0        273          7         17
# 8         112         12         22         25          2         33          0         20
# 9          65         10          0         15          4          9          2         12
# 10         23         14         32         17          0         16          0          6

# Cell Type Frequency per sample 
y <- table(Idents(immune.combined.integrated_1), immune.combined.integrated_1$orig.ident)
y <- as.data.frame(y)
write_csv(y, "/Documents/2021/Part2/celltype_frequency_per_sample.csv")
# make V1 an ordered factor
y$Var1 <- as.character(y$Var1)
#Turn your 'treatment' column into a character vector
y$Var2 <- as.character(y$Var2)
#Then turn it back into a factor with the levels in the correct order
y$Var2 <- factor(y$Var2, levels=c("sample_300", "sample_301", "sample_165","sample_295",
                                 "sample_326","sample_319","sample_155","sample_132"))



ggplot(y, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  xlab("Sample") +
  ylab("Proportion") +
  scale_fill_manual(values = brewer.pal(12, "Paired")) +
  theme(legend.title = element_blank())


# ## Find out Bcells in controls and MD group
# # Controls
# sum(y[which(y[,1]=="Bcells" & y[,2]==c("sample_300","sample_301","sample_165","sample_295")),3])
# # MD
# sum(y[which(y[,1]=="Bcells" & y[,2]==c("sample_326","sample_319","sample_155","sample_132")),3])

control_bcells <- "2060"
md_bcells <- "647"
df_bcells <- data.frame(celltype=c("control_bcells","md_bcells"),
                        number=c(2060,647))

# Basic barplot
p<-ggplot(data=df_bcells, aes(x=celltype, y=number)) +
  geom_bar(stat="identity", fill="steelblue")+
  geom_text(aes(label=number),vjust=1.6, color="white", size=3.5)+
  theme_minimal()
p

# Change barplot fill colors by groups
p<-ggplot(df_bcells, aes(x=celltype, y=number, fill=celltype)) +
  geom_bar(stat="identity")+
  geom_text(aes(label=number),vjust=1.6, color="white", size=3.5)+
  theme_minimal()
p

################################
########### MTRNR2L8 ###########
################################
# Percentage of cells expressing MTRNR2L8 per condition
#https:/www.biostars.org/p/9492208/#9492351
DefaultAssay(immune.combined.integrated_1) <- "RNA"
FeaturePlot(immune.combined.integrated_1, features = c("MTRNR2L8"), split.by = "stim",
            max.cutoff = 3 ,slot = "data",
            cols = c("grey", "red")) & theme(legend.position = "right")
features2 <- c("MTRNR2L8")
group_column <- "stim"
perc_exp <- DotPlot(immune.combined.integrated, features=features2, group.by=group_column)$data[, c("features.plot", "id", "pct.exp")]
perc_exp
library("ggplot2")
ggplot(perc_exp, aes(x=id, y=pct.exp, fill=id)) +
  geom_col() +
  facet_wrap(~features.plot)

# Percentage of cells expressing MTRNR2L8 per cluster
a <- DotPlot(object = immune.combined.integrated, features = features2)
a$data
ggplot(a$data, aes(x=id, y=pct.exp, fill=id)) +
  geom_col() +
  facet_wrap(~features.plot)

# Percentage of cells expressing MTRNR2L8 per cluster per condition
a <- DotPlot(object = immune.combined.integrated, features = features2, split.by = "stim")
a$data
write.csv(as.data.frame(a$data), "/Documents/2021/Part2/MTRNR2L8_pct.exp_all_clusters.csv")
ggplot(a$data, aes(x=id, y=pct.exp, fill=id)) +
  geom_col() + theme(axis.title.x=element_blank(),
                     axis.text.x=element_blank(),
                     axis.ticks.x=element_blank()) +
  facet_wrap(~features.plot) 

################################
########### IGHG1 #############
################################
DefaultAssay(immune.combined.integrated) <- "RNA"
features3 <- c("IGHG1")
FeaturePlot(immune.combined.integrated, features = features3, split.by = "stim", max.cutoff = 3, 
            cols = c("grey", "red"))
group_column <- "stim"
perc_exp <- DotPlot(immune.combined.integrated, features=features3, group.by=group_column)$data[, c("features.plot", "id", "pct.exp")]
perc_exp
library("ggplot2")
ggplot(perc_exp, aes(x=id, y=pct.exp, fill=id)) +
  geom_col() +
  facet_wrap(~features.plot)

# Percentage of cells expressing IGHG1 per cluster
a <- DotPlot(object = immune.combined.integrated, features = features3)
a$data
ggplot(a$data, aes(x=id, y=pct.exp, fill=id)) +
  geom_col() +
  facet_wrap(~features.plot)

# Percentage of cells expressing IGHG1 per cluster per condition
a <- DotPlot(object = immune.combined.integrated, features = features3, split.by = "stim")
a$data
#write.csv(as.data.frame(a$data), "/Documents/2021/New/Merge2/Merge2_clusters/2_8/MTRNR2L8_pct.exp_all_clusters.csv")
ggplot(a$data, aes(x=id, y=pct.exp, fill=id)) +
  geom_col() + theme(axis.title.x=element_blank(),
                     axis.text.x=element_blank(),
                     axis.ticks.x=element_blank()) +
  facet_wrap(~features.plot) 



# feature_scatter_integrated_test_plot <- FeatureScatter(immune.combined.integrated, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# feature_scatter_integrated_test_plot
# 
# DotPlot(immune.combined.integrated, features=c("nCount_RNA","percent.MT", "percent.Ribosomal"))


#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
# Subset B cells and Plasmacytoid dendritic cells/Bcells
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################

DefaultAssay(immune.combined.integrated_1) <- "integrated"
bells_part2 <- subset(x = immune.combined.integrated, idents = "4")
bells_part2 <- subset(x = immune.combined.integrated_1, idents = "Bcells")
bells_part2
dim(bells_part2)
#[1] 23477  2707
#23477 genes with 2707 Bcells
#[1] 2000 2707
unique(bells_part2@meta.data$seurat_clusters)
Idents(bells_part2)
# MTRNR2L8 expression in Bcell cluster
DefaultAssay(immune.combined.integrated_1) <- "RNA"
# Bar Plot
b_1_part2 <- DotPlot(object = immune.combined.integrated, features = "MTRNR2L8")
b_1_part2$data
ggplot(b_1_part2$data, aes(x=id, y=pct.exp, fill=id)) +
  geom_col() +
  facet_wrap(~features.plot)
# Split by condition
b_part2 <- DotPlot(object = immune.combined.integrated, features = "MTRNR2L8", split.by = "stim")
b_part2$data
ggplot(b_part2$data, aes(x=id, y=pct.exp, fill=id)) +
  geom_col() +
  facet_wrap(~features.plot)
# Feature Plot
FeaturePlot(bells_part2, features = "MTRNR2L8", split.by = "stim", max.cutoff = 3, 
            cols = c("grey", "red"))

DefaultAssay(bells_part2) <- "RNA"
bcell_sub <- FindVariableFeatures(bells_part2, selection.method = "vst", nfeatures = 2000)
top10_bcell_sub <- head(VariableFeatures(bcell_sub), 10)
top10_bcell_sub
bcell_sub_var1 <- VariableFeaturePlot(bcell_sub)

bcell_sub_var2 <- LabelPoints(plot = bcell_sub_var1, points = top10_bcell_sub, repel = TRUE)
bcell_sub_var2
bcell_sub_all_genes <- rownames(bcell_sub)
bcell_sub_scaled <- ScaleData(bcell_sub, features = bcell_sub_all_genes)
bcell_sub_PCA <- RunPCA(bcell_sub_scaled, features = VariableFeatures(object = bcell_sub_scaled))
bcell_sub_jack <- JackStraw(bcell_sub_PCA, num.replicate = 100)
bcell_sub_jack <- ScoreJackStraw(bcell_sub_jack, dims = 1:20)
ElbowPlot(bcell_sub_jack)
# t-SNE and Clustering
bcell_sub_UMAP <- RunUMAP(bcell_sub_PCA, reduction = "pca", dims = 1:30)
bcell_sub_nbr <- FindNeighbors(bcell_sub_UMAP, reduction = "pca", dims = 1:30)
# Resolution
bcell_sub_clusters <- FindClusters(bcell_sub_nbr, resolution = 0.3)
saveRDS(bcell_sub_clusters, file = "/Documents/2021/Part2/Bcell/Bcell_clusters_0.2.rds")
bcell_sub_clusters <- readRDS(file = "/Documents/2021/Part2/Bcell/Bcell_clusters_0.2.rds")

# Visualization
bcell_sub_p1 <- DimPlot(bcell_sub_clusters, reduction = "umap", group.by = "stim")
bcell_sub_p2 <- DimPlot(bcell_sub_clusters, reduction = "umap", label = TRUE)
bcell_sub_p2
bcell_sub_p1 + bcell_sub_p2


DimPlot(bcell_sub_clusters, reduction = "umap", split.by = "stim")
# All Markers
DefaultAssay(bcell_sub_clusters) <- "RNA"
bcell_sub_cluster_markers <- FindAllMarkers(bcell_sub_clusters, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(bcell_sub_cluster_markers,file = "/Documents/2021/Part3_paper/bcell_sub_clusters_0.2_all_markers.csv")
write.csv(bcell_sub_cluster_markers,file = "/Documents/2021/Part3_paper/bcell_sub_clusters_0.3_all_markers.csv")
write.csv(bcell_sub_cluster_markers,file = "/Documents/2021/Part3_paper/bcell_sub_clusters_0.4_all_markers.csv")
write.csv(bcell_sub_cluster_markers,file = "/Documents/2021/Part3_paper/bcell_sub_clusters_0.5_all_markers.csv")

bcell_sub_cluster_markers %>%
  group_by(cluster) %>% 
  slice_max(n=10, order_by = avg_log2FC) -> top10_bcell_sub_cluster_markers
top10_bcell_sub_cluster_markers
write.csv(top10_bcell_sub_cluster_markers,file = "/Documents/2021/Part3_paper/bcell_sub_clusters_0.2_top10_markers.csv")
write.csv(top10_bcell_sub_cluster_markers,file = "/Documents/2021/Part3_paper/bcell_sub_clusters_0.3_top10_markers.csv")
write.csv(top10_bcell_sub_cluster_markers,file = "/Documents/2021/Part3_paper/bcell_sub_clusters_0.4_top10_markers.csv")
write.csv(top10_bcell_sub_cluster_markers,file = "/Documents/2021/Part3_paper/bcell_sub_clusters_0.5_top10_markers.csv")

#write.csv(top10_bcell_sub_cluster_markers,file = "/Documents/2021/Part2/Bcell/bcell_sub_clusters_0.9_top10_2_markers.csv")
# DEG
output_directory <- "/Documents/2021/Part2/Bcell/"
DefaultAssay(bcell_sub_clusters) <- "RNA"
bcell_sub_cluster0 <- FindMarkers(bcell_sub_clusters, ident.1 = "CTRL", ident.2 = "MD", verbose = TRUE, group.by="stim", subset.ident = "0", logfc.threshold = log(2))
write.csv(bcell_sub_cluster0,file = paste0(output_directory,outputPrefix, "_bcell_sub_cluster_0_DEGs.csv"))
bcell_sub_cluster1 <- FindMarkers(bcell_sub_clusters, ident.1 = "CTRL", ident.2 = "MD", verbose = TRUE, group.by="stim", subset.ident = "1", logfc.threshold = log(2))
write.csv(bcell_sub_cluster1,file = paste0(output_directory,outputPrefix, "_bcell_sub_cluster_1_DEGs.csv"))
bcell_sub_cluster2 <- FindMarkers(bcell_sub_clusters, ident.1 = "CTRL", ident.2 = "MD", verbose = TRUE, group.by="stim", subset.ident = "2", logfc.threshold = log(2))
write.csv(bcell_sub_cluster2,file = paste0(output_directory,outputPrefix, "_bcell_sub_cluster_2_DEGs.csv"))
bcell_sub_cluster3 <- FindMarkers(bcell_sub_clusters, ident.1 = "CTRL", ident.2 = "MD", verbose = TRUE, group.by="stim", subset.ident = "3", logfc.threshold = log(2))
write.csv(bcell_sub_cluster3,file = paste0(output_directory,outputPrefix, "_bcell_sub_cluster_3_DEGs.csv"))
bcell_sub_cluster4 <- FindMarkers(bcell_sub_clusters, ident.1 = "CTRL", ident.2 = "MD", verbose = TRUE, group.by="stim", subset.ident = "4", logfc.threshold = log(2))
write.csv(bcell_sub_cluster4,file = paste0(output_directory,outputPrefix, "_bcell_sub_cluster_4_DEGs.csv"))
bcell_sub_cluster5 <- FindMarkers(bcell_sub_clusters, ident.1 = "CTRL", ident.2 = "MD", verbose = TRUE, group.by="stim", subset.ident = "5", logfc.threshold = log(2))
write.csv(bcell_sub_cluster5,file = paste0(output_directory,outputPrefix, "_bcell_sub_cluster_5_DEGs.csv"))
bcell_sub_cluster6 <- FindMarkers(bcell_sub_clusters, ident.1 = "CTRL", ident.2 = "MD", verbose = TRUE, group.by="stim", subset.ident = "6", logfc.threshold = log(2))
write.csv(bcell_sub_cluster6,file = paste0(output_directory,outputPrefix, "_bcell_sub_cluster_6_DEGs.csv"))
bcell_sub_cluster7 <- FindMarkers(bcell_sub_clusters, ident.1 = "CTRL", ident.2 = "MD", verbose = TRUE, group.by="stim", subset.ident = "7", logfc.threshold = log(2))
write.csv(bcell_sub_cluster7,file = paste0(output_directory,outputPrefix, "_bcell_sub_cluster_7_DEGs.csv"))


## BCells and Plasma
b_plasma <- subset(x = immune.combined.integrated_1, idents = c("Bcells","Plasma"))
b_plasma
dim(b_plasma)
#[1] 23477  3483
Idents(b_plasma)

DefaultAssay(b_plasma) <- "RNA"

b_plasma <- FindVariableFeatures(b_plasma, selection.method = "vst", nfeatures = 2000)
# top10_b_plasma_sub <- head(VariableFeatures(b_plasma), 10)
# top10_b_plasma_sub
# b_plasma_sub_var1 <- VariableFeaturePlot(b_plasma)
# 
# bcell_sub_var2 <- LabelPoints(plot = bcell_sub_var1, points = top10_bcell_sub, repel = TRUE)
# bcell_sub_var2
# bcell_sub_all_genes <- rownames(bcell_sub)
b_plasma_sub_all_genes <- rownames(b_plasma)
b_plasma_sub_scaled <- ScaleData(b_plasma, features = b_plasma_sub_all_genes)
b_plasma_sub_PCA <- RunPCA(b_plasma_sub_scaled, features = VariableFeatures(object = b_plasma_sub_scaled))
# t-SNE and Clustering
b_plasma_sub_UMAP <- RunUMAP(b_plasma_sub_PCA, reduction = "pca", dims = 1:20)
b_plasma_sub_nbr <- FindNeighbors(b_plasma_sub_UMAP, reduction = "pca", dims = 1:20)
# Resolution
b_plasma_sub_clusters <- FindClusters(b_plasma_sub_nbr, resolution = 0.5)
saveRDS(b_plasma_sub_clusters, file = "/Documents/2021/Part2/Bcell/Bcell_plasma_clusters_0.5.rds")
b_plasma_sub_clusters <- readRDS(file = "/Documents/2021/Part2/Bcell/Bcell_plamsa_clusters_0.5.rds")
# Visualization
b_plasma_sub_p1 <- DimPlot(b_plasma_sub_clusters, reduction = "umap", group.by = "stim")
b_plasma_sub_p2 <- DimPlot(b_plasma_sub_clusters, reduction = "umap", label = TRUE)
b_plasma_sub_p1 + b_plasma_sub_p2
b_plasma_sub_p2

DimPlot(b_plasma_sub_clusters, reduction = "umap", split.by = "stim")
# All Markers
DefaultAssay(b_plasma_sub_clusters) <- "RNA"
b_plasma_sub_cluster_markers <- FindAllMarkers(b_plasma_sub_clusters, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
b_plasma_sub_cluster_markers %>%
  group_by(cluster) %>% 
  slice_max(n=10, order_by = avg_log2FC) -> top10_bcell_plasma_sub_cluster_markers
top10_bcell_plasma_sub_cluster_markers
write.csv(top10_bcell_plasma_sub_cluster_markers,file = "/Documents/2021/Part2/Bcell_Plasma/bcell_plamsa_sub_clusters_2_0.5_top10_2_markers.csv")

######################################################################################
######################################################################################
######################################################################################
######################################################################################

# Subset on the expression level of a gene/feature
DefaultAssay(immune.combined.integrated_1) <- "RNA"
#DefaultAssay(immune.combined.integrated_1) <- "integrated"
bplasma_part3 <- subset(x = immune.combined.integrated_1, subset = CD79A >3 | 
                          MS4A1 > 3 | MZB1 >3  ) 
dim(bplasma_part3)
#[1] 23477  2707

bplasma_part3 <- subset(x = immune.combined.integrated_1, idents = c("Bcells","Plasma") ) 
bplasma_part3
dim(bplasma_part3)
# > dim(bplasma_part3)
# [1] 23453  3483
bplasma_part3 <- subset(x = immune.combined.integrated_1, idents = c("Bcells") ) 
bplasma_part3
dim(bplasma_part3)
# [1] 23477  2707

unique(bplasma_part3@meta.data$seurat_clusters)
Idents(bplasma_part3)
# Filter MALAT1
bplasma_part3 <- bplasma_part3[!grepl("MALAT1", rownames(bplasma_part3)), ]
dim(bplasma_part3)
# 23476  1536
# Filter Mitocondrial
bplasma_part3 <- bplasma_part3[!grepl("^MT-", rownames(bplasma_part3)), ]
dim(bplasma_part3)
# 23463  1536
# # Filter Ribossomal gene (optional if that is a problem on your data) data.filt
# data.filt  <- bplasma_part3[!grepl("^RPL", rownames(bplasma_part3)), ]
# dim(bplasma_part3)
# # Filter Ribossomal gene (optional if that is a problem on your data) data.filt
# data.filt  <- bplasma_part3[!grepl("^RPS", rownames(bplasma_part3)), ]
# dim(bplasma_part3)
# # [1] 23463  1536
# Filter Hemoglobin gene (optional if that is a problem on your data)
bplasma_part3 <- bplasma_part3[!grepl("^HB[^(P)]", rownames(bplasma_part3)), ]
dim(bplasma_part3)
# [1] 23453  1536
# # Filter Erythroid precursor cells (optional if that is a problem on your data)
# bplasma_part3 <- bplasma_part3[!grepl("ADAM10", rownames(bplasma_part3)), ]
# bplasma_part3 <- bplasma_part3[!grepl("SLC25A37", rownames(bplasma_part3)), ]
# bplasma_part3 <- bplasma_part3[!grepl("SEMA3E", rownames(bplasma_part3)), ]
# bplasma_part3 <- bplasma_part3[!grepl("XIST", rownames(bplasma_part3)), ]
# bplasma_part3 <- bplasma_part3[!grepl("NABP1", rownames(bplasma_part3)), ]
# bplasma_part3 <- bplasma_part3[!grepl("GPM6A", rownames(bplasma_part3)), ]
# bplasma_part3 <- bplasma_part3[!grepl("KCNJ3", rownames(bplasma_part3)), ]
# bplasma_part3 <- bplasma_part3[!grepl("MT-ND6", rownames(bplasma_part3)), ]
# bplasma_part3 <- bplasma_part3[!grepl("STAG2", rownames(bplasma_part3)), ]
# dim(bplasma_part3)
DefaultAssay(bplasma_part3) <- "RNA"

bplasma_part3 <- FindVariableFeatures(bplasma_part3, selection.method = "vst", nfeatures = 2000)
# top10_b_plasma_sub <- head(VariableFeatures(b_plasma), 10)
# top10_b_plasma_sub
# b_plasma_sub_var1 <- VariableFeaturePlot(b_plasma)
# 
# bcell_sub_var2 <- LabelPoints(plot = bcell_sub_var1, points = top10_bcell_sub, repel = TRUE)
# bcell_sub_var2
# bcell_sub_all_genes <- rownames(bcell_sub)
b_plasma_sub_all_genes <- rownames(bplasma_part3)
b_plasma_sub_scaled <- ScaleData(bplasma_part3, features = b_plasma_sub_all_genes)
b_plasma_sub_PCA <- RunPCA(b_plasma_sub_scaled, features = VariableFeatures(object = b_plasma_sub_scaled))
# t-SNE and Clustering
b_plasma_sub_UMAP <- RunUMAP(b_plasma_sub_PCA, reduction = "pca", dims = 1:20)
b_plasma_sub_nbr <- FindNeighbors(b_plasma_sub_UMAP, reduction = "pca", dims = 1:20)
# Resolution
b_plasma_sub_clusters <- FindClusters(b_plasma_sub_nbr, resolution = 0.2)
# saveRDS(b_plasma_sub_clusters, file = "/Documents/2021/Part2/Bcell/Bcell_plasma_filtered_clusters_0.5.rds")
# b_plasma_sub_clusters <- readRDS(file = "/Documents/2021/Part2/Bcell/Bcell_plamsa__filtered_clusters_0.5.rds")
saveRDS(b_plasma_sub_clusters, file = "/Documents/2021/Part2/Bcell_Plasma/Subset_CD79A_MS4A1_MZB1/Bcell_plasma_filtered_clusters_0.5.rds")
b_plasma_sub_clusters <- readRDS(file = "/Documents/2021/Part2/Bcell_Plasma/Subset_CD79A_MS4A1_MZB1/Bcell_plamsa__filtered_clusters_0.5.rds")

# Visualization
b_plasma_sub_p1 <- DimPlot(b_plasma_sub_clusters, reduction = "umap", group.by = "stim")
b_plasma_sub_p2 <- DimPlot(b_plasma_sub_clusters, reduction = "umap", label = TRUE)
b_plasma_sub_p1 + b_plasma_sub_p2
b_plasma_sub_p2

DimPlot(b_plasma_sub_clusters, reduction = "umap", split.by = "stim")
# All Markers
DefaultAssay(b_plasma_sub_clusters) <- "RNA"
b_plasma_sub_cluster_markers <- FindAllMarkers(b_plasma_sub_clusters, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(b_plasma_sub_cluster_markers,file = "/Documents/2021/Part2/Bcell_Plasma/Part2/bcell_plamsa_filtered_sub_clusters_0.5_markers.csv")
write.csv(b_plasma_sub_cluster_markers,file = "/Documents/2021/Part2/Bcell_Plasma/Subset_CD79A_MS4A1_MZB1/bcell_plamsa_filtered_sub_clusters_0.5_markers.csv")
write.csv(b_plasma_sub_cluster_markers,file = "/Documents/2021/Part2/Bcell_Plasma/Part2/Filtered_integrated/bcell_plamsa_filtered_sub_clusters_0.5_markers.csv")
write.csv(b_plasma_sub_cluster_markers,file = "/Documents/2021/Part2/Bcell_Plasma/Part2/Filtered_integrated/bcell_plamsa_filtered_sub_clusters_0.1_markers.csv")
write.csv(b_plasma_sub_cluster_markers,file = "/Documents/2021/Part2/Bcell_Plasma/Part2/Resolution_0.2/bcell_plamsa_filtered_sub_clusters_0.2_markers.csv")
write.csv(b_plasma_sub_cluster_markers,file = "/Documents/2021/Part2/Bcell_Plasma/Part2/Resolution_0.4/bcell_plamsa_filtered_sub_clusters_0.4_markers.csv")
write.csv(b_plasma_sub_cluster_markers,file = "/Documents/2021/Part2/Bcell/Resolution_0.2/bcell_plamsa_filtered_sub_clusters_0.2_markers.csv")
write.csv(b_plasma_sub_cluster_markers,file = "/Documents/2021/Part2/Bcell/Resolution_0.3/bcell_plamsa_filtered_sub_clusters_0.3_markers.csv")
write.csv(b_plasma_sub_cluster_markers,file = "/Documents/2021/Part2/Bcell/Resolution_0.4/bcell_plamsa_filtered_sub_clusters_0.4_markers.csv")
write.csv(b_plasma_sub_cluster_markers,file = "/Documents/2021/Part2/Bcell/Resolution_0.7/bcell_plamsa_filtered_sub_clusters_0.7_markers.csv")

b_plasma_sub_cluster_markers %>%
  group_by(cluster) %>% 
  slice_max(n=10, order_by = avg_log2FC) -> top10_bcell_plasma_sub_cluster_markers
top10_bcell_plasma_sub_cluster_markers
write.csv(top10_bcell_plasma_sub_cluster_markers,file = "/Documents/2021/Part2/Bcell_Plasma/Part2/bcell_plamsa_filtered_sub_clusters_0.5_top10_markers.csv")
write.csv(top10_bcell_plasma_sub_cluster_markers,file = "/Documents/2021/Part2/Bcell_Plasma/Subset_CD79A_MS4A1_MZB1/bcell_plamsa_filtered_sub_clusters_0.5_top10_markers.csv")
write.csv(top10_bcell_plasma_sub_cluster_markers,file = "/Documents/2021/Part2/Bcell_Plasma/Part2/Filtered_integrated/bcell_plamsa_filtered_sub_clusters_0.5_top10_markers.csv")
write.csv(top10_bcell_plasma_sub_cluster_markers,file = "/Documents/2021/Part2/Bcell_Plasma/Part2/Resolution_0.2/bcell_plamsa_filtered_sub_clusters_0.2_top10_markers.csv")


FeaturePlot(b_plasma_sub_clusters, features = c('SLC25A37'))
FeaturePlot(b_plasma_sub_clusters, features = c('CD79A'))
FeaturePlot(b_plasma_sub_clusters, features = c('CD27'))
FeaturePlot(b_plasma_sub_clusters, features = c('TCL1A'))
FeaturePlot(b_plasma_sub_clusters, features = c('FCER2'))

features = c('CD24','TCL1A', 'MS4A1', 'MTRNR2L8',
             'CD27', 'MZB1', 'ITGAX','IGHA2',
             'RHOB','IL2RA', 'IGHM','CD79A', 
             'JCHAIN', 'CD93',"HLA-DRA", 
             "FCER2", "CD69", "IGHD", "CD19")
features2 = c('TCL1A', 'MS4A1','IGHM',"HLA-DRA",
             'CD27', 'MZB1','IGHA2',"IGLV2-23", "IGHV4-39",
             'IL2RA','CD79A', 
             'JCHAIN', "FCER2", "IGHD")

DotPlot(bcell_sub_clusters, features = features2) + RotatedAxis()
FeaturePlot(bcell_sub_clusters, features = features)
# CD24 -> naïve Bcells
# CD27 -> classical_memoryBcells
# IGHD ->classical_memoryBcells
# CD19 -> classical_memoryBcells
# TCL1A -> transitional Bcells
# MZB1 -> Plasma cells
# IGHA2 -> DN1
# JCHAIN  -> DN1
# ITGAX -> DN2
# RHOB -> DN2_3
# FCER2 -> DN4 
# IL2RA -> matureBcells
# Ifitm2 , Ifitm3 -> proBcells
# IGHM -> naive_memoryBcell
# CD79A -> all_Bcells
# HLA-DRA -> plasmablast
# CD53 -> BCell develoment regulator

features = c("CD24","CD27","TCL1A","MZB1")
# Ridge plots - from ggridges. Visualize single cell expression distributions in each cluster
RidgePlot(b_plasma_sub_clusters, features = features, ncol = 2)

# bplasma_part4 <- subset(x = b_plasma_sub_clusters, idents = c(1,2,3,4)) 
# bplasma_part4 <- FindVariableFeatures(bplasma_part4, selection.method = "vst", nfeatures = 2000)
# # top10_b_plasma_sub <- head(VariableFeatures(b_plasma), 10)
# # top10_b_plasma_sub
# # b_plasma_sub_var1 <- VariableFeaturePlot(b_plasma)
# # 
# # bcell_sub_var2 <- LabelPoints(plot = bcell_sub_var1, points = top10_bcell_sub, repel = TRUE)
# # bcell_sub_var2
# # bcell_sub_all_genes <- rownames(bcell_sub)
# b_plasma_sub_all_genes <- rownames(bplasma_part4)
# b_plasma_sub_scaled <- ScaleData(bplasma_part3, features = b_plasma_sub_all_genes)
# b_plasma_sub_PCA <- RunPCA(b_plasma_sub_scaled, features = VariableFeatures(object = b_plasma_sub_scaled))
# # t-SNE and Clustering
# b_plasma_sub_UMAP <- RunUMAP(b_plasma_sub_PCA, reduction = "pca", dims = 1:20)
# b_plasma_sub_nbr <- FindNeighbors(b_plasma_sub_UMAP, reduction = "pca", dims = 1:20)
# # Resolution
# b_plasma_sub_clusters <- FindClusters(b_plasma_sub_nbr, resolution = 0.5)
# #saveRDS(b_plasma_sub_clusters, file = "/Documents/2021/Part2/Bcell/Bcell_plasma_filtered_clusters_0.5.rds")
# #b_plasma_sub_clusters <- readRDS(file = "/Documents/2021/Part2/Bcell/Bcell_plamsa__filtered_clusters_0.5.rds")
# # Visualization
# b_plasma_sub_p1 <- DimPlot(b_plasma_sub_clusters, reduction = "umap", group.by = "stim")
# b_plasma_sub_p2 <- DimPlot(b_plasma_sub_clusters, reduction = "umap", label = TRUE)
# b_plasma_sub_p1 + b_plasma_sub_p2
# b_plasma_sub_p2
# 
# DimPlot(b_plasma_sub_clusters, reduction = "umap", split.by = "stim")
# # All Markers
# DefaultAssay(b_plasma_sub_clusters) <- "RNA"
# b_plasma_sub_cluster_markers <- FindAllMarkers(b_plasma_sub_clusters, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# #write.csv(b_plasma_sub_cluster_markers,file = "/Documents/2021/Part2/Bcell_Plasma/Part2/bcell_plamsa_filtered_sub_clusters_0.5_markers.csv")
# b_plasma_sub_cluster_markers %>%
#   group_by(cluster) %>% 
#   slice_max(n=10, order_by = avg_log2FC) -> top10_bcell_plasma_sub_cluster_markers
# top10_bcell_plasma_sub_cluster_markers
# #write.csv(top10_bcell_plasma_sub_cluster_markers,file = "/Documents/2021/Part2/Bcell_Plasma/Part2/bcell_plamsa_filtered_sub_clusters_0.5_top10_markers.csv")
