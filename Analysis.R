setwd("~/Documents/Projects/")
library(tidyverse)
##install.packages("SoupX")
library(SoupX)
library(Seurat)
library(parallel)
library(ggplot2)
library(cowplot)
library(grid)
library(gridExtra)

##install.packages("hdf5r")
##library(hdf5r)
library(RColorBrewer)
`%nin%` <- Negate(`%in%`)
getPalette <- colorRampPalette(brewer.pal(9, "Set1"))

dataDir <- "./FL/"

samples_names <- list.files("./FL")
samples_dir <- list.files("./FL", full.names = T)
# sc = load10X(samples_dir[5])
# 
filt.matrix  = Seurat::Read10X(file.path(samples_dir[5], "filtered_feature_bc_matrix"))
raw.matrix = Seurat::Read10X(file.path(samples_dir[5], "raw_feature_bc_matrix"))
filt_genes <- rownames(filt.matrix)
# Subset raw.matrix to keep only the genes in filt.matrix
raw.matrix_subset <- raw.matrix[rownames(raw.matrix) %in% filt_genes, ]
sc = SoupChannel(raw.matrix_subset, filt.matrix)
seurat_obj <- Read10X(paste0(samples_dir[5],"/filtered_feature_bc_matrix/"))

seurat_obj <- CreateSeuratObject(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj), verbose = FALSE)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30, verbose = FALSE)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5, verbose = FALSE)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:30, verbose = FALSE)

# addin clusters and DR to soupx object
clusters <- data.frame("clusters" = seurat_obj$seurat_clusters, row.names = colnames(seurat_obj))
umap <- as.data.frame(seurat_obj@reductions$umap@cell.embeddings)
sc = setClusters(sc, setNames(clusters$clusters, rownames(clusters)))
sc = setDR(sc, umap)

if (missing(rho)){
  
  acinar_markers <- c('CTRB1','KLK1','RBPJL','PTF1A','CELA3A','PRSS1','SPINK1','ZG16','CEL','CELA2A','CPB1','CELA1','RNASE1','AMY2B','CPA2','CPA1','CELA3B','PNLIP','CTRB2','PLA2G1B','PRSS2','CLPS','REG1A','SYCN','PNLIPRP1','CTRC','REG3A','SERPINA3','PRSS3','REG1B','CFB','GDF15','MUC1','C15ORF48','DUOXA2','AKR1C3','OLFM4','GSTA1','LGALS2','PDZK1IP1','RARRES2','CXCL17','UBD','GSTA2','ANPEP','LYZ','ANGPTL4','ALDOB')
  acinar_markers <- acinar_markers[acinar_markers %in% row.names(sc$toc)]
  rbc_markers <- c('HEXA','AQP1','CCNA2','PKLR','LOX','ARG1','HBB','HBA1','HBA2','HBG1','HBG2','CEACAM1','CD36','GYPA','THBS1','ITGA4','GYPB','AHSP','HBM','FECH','EPOR','KIT','BPGM','DCAF12','HBD','HEMGN','HMBS','NCOA4','RBM38','UCP2','UBB','TMCC2','SLC4A1','SLC25A39','SLC25A37','FOSB','HTATSF1','EPB41','RSAD2','COX6B2','APOA1','ASNS','ALAD','S100A9','AHSG','S100A8','ERMAP','KCNN4','PAN3','PFAS','VAMP5','PNPO')
  rbc_markers <- rbc_markers[rbc_markers %in% row.names(sc$toc)]
  beta_markers <- c('GCGR','JPH3','CD40','HAMP','EZH1','NTRK1','PDX1','SLC2A2','NKX6-2','FXYD2','NPY','INS','RIMS1','MAFA','EFNA5','LMX1A','NKX2-2','NKX6-1','PAX4','IAPP','PCSK2','G6PC2','SLC30A8','PCSK1','GJD2','SCGN','IGF2','SYT13','FFAR2','NPTX2','PFKFB2','EDARADD','HOPX','SH3GL2','ADCYAP1','SCGB2A1','CASR','MAFB','PAX6','NEUROD1','ISL1','TGFBR3','SMAD9','SIX3','SIX2','BMP5','PIR','STXBP5','DLK1','MEG3','RGS16')
  beta_markers <- beta_markers[beta_markers %in% row.names(sc$toc)]
  
  nonExpressedGeneList = list(acinar_markers  = acinar_markers,
                              rbc_markers = rbc_markers,
                              beta_markers = beta_markers)
  
  useToEst = estimateNonExpressingCells(sc, nonExpressedGeneList = nonExpressedGeneList)
  
  ##sc = calculateContaminationFraction(sc, nonExpressedGeneList = nonExpressedGeneList, useToEst = useToEst)
  print(paste0("The estimated contamination fraction is: ", round(sc$metaData$rho[1],3)))
  
  if (missing(extra)){
    
    note = "no extra value was added to the contamination fraction as it wasn't provided"
    if (cap){
      new_rho = min(sc$metaData$rho[1],0.2)
      sc = setContaminationFraction(sc, new_rho)
      out = adjustCounts(sc, roundToInt = T, verbose = F)
    } else {out = adjustCounts(sc, roundToInt = T, verbose = F)}
    
  } else {
    
    new_rho = sc$metaData$rho[1]+extra
    note = "added extra 10% to the contamination fraction as provided"
    if (cap){
      new_rho = min(new_rho,0.2)
      sc = setContaminationFraction(sc, new_rho)
      out = adjustCounts(sc, roundToInt = T, verbose = F)
    } else {
      sc = setContaminationFraction(sc, new_rho, forceAccept = T)
      out = adjustCounts(sc, roundToInt = T, verbose = F)
    }
  }
  
  # if (sc$metaData$rho[1] > 0.2) {
  #   note = "no extra value was added to the contamination fraction, Contamination fraction is already high (> 20%)"
  #   new_rho = min(sc$metaData$rho[1],0.2)
  #   sc = setContaminationFraction(sc, new_rho)
  #   out = adjustCounts(sc, roundToInt = T, verbose = F)
  
  # } else {
  # print(paste0("Adding ",extra," to the estimated fraction"))
  # new_rho = sc$metaData$rho[1]+extra
  # sc = setContaminationFraction(sc, new_rho)
  # note = "added extra 10% to the contamination fraction as provided"
  # out = adjustCounts(sc, roundToInt = T, verbose = F)
  # }
  
  
} else {
  
  sc = setContaminationFraction(sc, rho)
  out = adjustCounts(sc, roundToInt = T, verbose = F)
  note = "Contamination Fraction was provided by the user"
  
}

TopContGenes <- data.frame("Genes" = row.names(sc$soupProfile),
                           "counts" = sc$soupProfile$counts,
                           "fraction" = sc$soupProfile$est)

Misc(seurat_obj, "TopContGenes") <- TopContGenes[order(TopContGenes$counts, decreasing = T),]
Misc(seurat_obj, "ContaminationFraction") <- round(sc$metaData$rho[1],3)
Misc(seurat_obj, "notes") <- note

seurat_obj[['CorrectedCounts']] <- CreateAssayObject(counts = out)
DefaultAssay(seurat_obj) <- 'CorrectedCounts'

#################################################################
##           Ambient RNA Correction (w +10% and Cap)           ##
#################################################################

clust <- makeCluster(detectCores())
clusterExport(clust,
              c("correct_reads"),
              envir=environment())
corrected_samples <-
  parLapply(clust,
            samples_dir,
            fun = function(x){
              print(x)
              correct_reads(data_dir = x, rho = 0.2, extra = 0.1, cap = TRUE, h5 = FALSE)}
  )

##correct_reads(data_dir = samples_dir[5], rho = 0.2, extra = 0.1, cap = TRUE, h5 = FALSE)
stopCluster(clust)

dir.create("objects", showWarnings = F)
names(corrected_samples) <- samples_names
saveRDS(object = corrected_samples, file = "./results/Fl_corrected_samples_w_extra_10_cap.rds")
remove(corrected_samples)

##sc = load10X(dataDir, verbose = F)

#################################################################
##                 Loading and Processing Data                 ##
#################################################################

corrected_samples <- readRDS("./results/Fl_corrected_samples_w_extra_10_cap.rds")

## Adding Metadata
# for (sample in names(corrected_samples)){
#   corrected_samples[[sample]]$sample_name <- samples_info[sample, "sample_name"]
#   corrected_samples[[sample]]$patient_id <- samples_info[sample, "patient_id"]
#   corrected_samples[[sample]]$DiseaseState <- samples_info[sample, "DiseaseState"]
#   corrected_samples[[sample]]$runId <- samples_info[sample, "runId"]
#   corrected_samples[[sample]]$location <- samples_info[sample, "location"]
#   corrected_samples[[sample]]$lesion <- samples_info[sample, "lesion"]
# }

names(corrected_samples)
corrected_samples[[1]]$sample_name = "D"
corrected_samples[[2]]$sample_name = "A"
corrected_samples[[3]]$sample_name = "B"
corrected_samples[[4]]$sample_name = "C"
corrected_samples[[5]]$sample_name = "E"
corrected_samples[[6]]$sample_name = "F"

corrected_samples[[3]] = NULL
corrected_samples$EIPA_PD1
## Merging all samples into one object
merged_samples <- corrected_samples[[1]]
corrected_samples[[1]] = NULL
corrected_samples$
merged_samples <- merge(merged_samples, corrected_samples[[2]])
corrected_samples[[2]] = NULL
corrected_samples$
merged_samples <- merge(merged_samples, corrected_samples[[3]])
merged_samples <- merge(merged_samples, corrected_samples[[4]])
merged_samples$sample_name[[2]]

rm(merged_samples)

merged_samples <- corrected_samples[[1]]
for (i in 2:length(corrected_samples)){
  merged_samples <- merge(merged_samples, corrected_samples[[i]])
}
## Saving merged object
merged_samples = readRDS(file = "./New Folder With Items/Gift-of-Life-Public-Repository-master/results/Fl_merged_samples_corrected_w_extra_10.rds")

#################################################################
##                           Data QC                           ##
#################################################################

merged_samples[['percent.mt_RNA']] <- PercentageFeatureSet(merged_samples, pattern = "^MT-", assay = "RNA")
merged_samples[['percent.mt_CorrectedCounts']] <- PercentageFeatureSet(merged_samples, pattern = "^MT-", assay = "CorrectedCounts")
Idents(merged_samples) <- "all_samples"

VlnPlot(merged_samples, features = c("nFeature_RNA","nFeature_CorrectedCounts"), pt.size = 0)
VlnPlot(merged_samples, features = c("nCount_RNA","nCount_CorrectedCounts"), pt.size = 0)
VlnPlot(merged_samples, features = c("percent.mt_RNA","percent.mt_CorrectedCounts"), pt.size = 0)
VlnPlot(merged_samples, group.by = "sample_name", features = "percent.mt_CorrectedCounts", pt.size = 0)

merged_samples <- subset(merged_samples, subset = nFeature_CorrectedCounts > 200 & percent.mt_CorrectedCounts < 15)

## Saving merged object
DefaultAssay(merged_samples) <- 'CorrectedCounts'

##################################################################
##                       rPCA Integration                       ##
##################################################################
unique(merged_samples@meta.data$sample_name)
samples.list <- SplitObject(merged_samples, split.by = "sample_name")
samples.list <- lapply(X = samples.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
features <- SelectIntegrationFeatures(object.list = samples.list)
samples.list <- lapply(X = samples.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

anchors <- FindIntegrationAnchors(object.list = samples.list, anchor.features = features, reduction = "rpca")
rm(merged_samples)
rm(samples.list)
rm(corrected_samples)
rpca_integrated_by_samples <- IntegrateData(anchorset = anchors)
DefaultAssay(rpca_integrated_by_samples) <- "integrated"
rpca_integrated_by_samples <- ScaleData(rpca_integrated_by_samples, verbose = FALSE)
rpca_integrated_by_samples <- RunPCA(rpca_integrated_by_samples, verbose = FALSE)
rpca_integrated_by_samples <- RunUMAP(rpca_integrated_by_samples, reduction = "pca", dims = 1:30)
rpca_integrated_by_samples <- FindNeighbors(rpca_integrated_by_samples, reduction = "pca", dims = 1:30)
rpca_integrated_by_samples <- FindClusters(rpca_integrated_by_samples, resolution = 0.5)
integrated_data = rpca_integrated_by_samples
DefaultAssay(integrated_data) <- "CorrectedCounts"
rm(rpca_integrated_by_samples)
markers <- FindAllMarkers(integrated_data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(markers, "Fl_markers.csv")

##-----------------------------------------
##  FeaturePlots for known markers
##-----------------------------------------

p = DimPlot(integrated_data, label = T, cols = getPalette(length(levels(integrated_data))),
        label.box = T, repel = T, raster = T) + ggtitle("Seurat Clusters")

ggsave(plot=p, filename=file.path(dataDir, "Fl_integrated_clusters_annot.pdf"), units = "in", width=10, height= 7, dpi=600)
##DimPlot(integrated_data, label = F, group.by = "DiseaseState", raster = T, order = T) + ggtitle("DiseaseState")
p1 = DimPlot(integrated_data, label = F, group.by = "sample_name", repel = T, raster = T, order = T) + ggtitle("Sample Name")
ggsave(plot=p1, filename=file.path(dataDir, "Flex_integrated_clusters_samplewise.pdf"), units = "in", width=10, height= 7, dpi=600)
DefaultAssay(integrated_data) <- "CorrectedCounts"
# features <- c("SPINK1", "PRSS3", "CTRB1", "PRSS1", "AMY2A",
#               "EPCAM","KRT19", "SOX9", "MUC1", "KRT18",
#               "CDH5", "VWF",
#               "PDGFRB", "PDGFRA", "PDPN", "ACTA2","TAGLN", "DPT", "RGS5",
#               "PTPRC", "CD14", "ITGAM", "MARCO", "APOE", "C1QA","CD68",
#               "CD3E", "CD4","CD8A",
#               "CD19", "CD79A",'MS4A1',
#               "HBA2", "HBA1",
#               "MKI67",
#               "INS",
#               'NCAM1', 'NEGR1', 'NRN1')

# Myh11", "Notch3", "Dcn", "Col1a1", "Col3a1",
#               "Col1a2","Cd3d", "Cd8a", "Nkg7", "Gzma",
#               "Siglecg", "H2-Ob",
#               "Itgae", "Cxcr2", "Mmp9", "Csf3r","Itgam", "Lyz2", "C1qb",
#               "Apoe", "Plvap", "Cldn5", "Pecam1", "Ube2c", "Mki67","Top2a",
#               "Ifi44", "Ifit3b","Ifit1",
#               "Ifit3", "Krt7",'Krt19',
#               "Epcam", "Muc1", "Chac1", "Trib3","Luc7l2", "Luc7l3", "Vegfa",
#               "Trf",
#               "Spp1",
#               "Aldh2",

features <- unique(c( 
  "ductal cell 1", "Ambp", "Cftr", "Mmp7",
              "ductal cell 2", "Krt19", "Krt7", "Tspan8", "Slpi",
              "acinar", "Prss1", "Ctrb1", "Ctrb2", "Reg1b",
              "endocrine", "Chgb", "Chga", "Ins", "Iapp",
              "Stellate cell", "Rgs5", "Acta2", "Pdgfrb", "Adirf",
              "Fibroblast", "Lum", "Dcn", "Col1a1","Cdh11","Pdgfra","Pdgfrb","Acta2", 'Col3a1', 
              'Rgs5', 'Igfbp7', 'Pdpn','Mcam','Il6','Apoe','Gli1','Gli2','Gli3','Pdgfa',
              "Endothelial","Cdh5", "Plvap", "Vwf", "Cldn5",
  "Macrophage", "Aif1", "Cd64", "Cd14", "Cd68",
  "Tcell", "Cd3d", "Cd3e", "Cd4", "Cd8",
  "Bcell", "Ms4a1", "Cd79a", "Cd79b", "Cd52",
  "Epi_Markers","Krt7","Krt8","Krt18","Krt19","Epcam","Cdh1",'Prss11', 'Ctrb2','Reg1a','Clu','Mki67','Krt8','Spink1','Krt19','Krt18','Tff1','Muc1',
  "T_Cell_Markers","Cd3e","Cd3g","Cd3d","Cd4","Il7r","Cd8a","Lef1",
  "Myeloid_Markers","Cd14","Itgam","Mnda","Mpeg1","Itgax",'Fcgr3a','Fcgr3b','Apoe',"C1qa","Marco","Lyz","Hla-Dra",
  "B_Cell_Markers","Cd79a","Ms4a1","Cd19",
  "RBC_Markers","Hba1","Hbb","Hba2",
  "NK_Markers","Ncr3","Fcgr3a","Ncam1","Klrf1","Klrc1","Cd38","Klrc1","Nkg7", "Gzmb","Gzma",
  "MAST CELLS",'Tpsab1','Cpa3',
  "Dendritic cells",'Itgae','Lyz','Clec9a','Batf3','Irf8','Ido1','Cd207','Cd1a'
  ,'Cd1C', 'Hla-Dra','Ccl22','Lamp3','Il22ra2','Cd101',
  "Trf", "Spp1","Aldh2",
  "Luc7l2", "Luc7l3", "Vegfa",
  "Chac1", "Trib3","Aldh112",
  "Ube2c", "Mki67","Top2a"
  
  
  
  ))

p2 = DotPlot(integrated_data, features = features, assay = 'CorrectedCounts') + coord_flip()
ggsave(plot=p2, filename=file.path(dataDir, "F1_integrated_dot_plot_markers.pdf"), units = "in", width=10, height= 14, dpi=600)
saveRDS(object = integrated_data, file = "./results/Fl_integrated_data_w_extra_10_cap.rds")

#################################################################
##                      Renaming Clusters                      ##
#################################################################

cluster_annotation <- readxl::read_xlsx("./cluster_annotation.xlsx", sheet = 3)

cell_type <- cluster_annotation$cell_type
names(cell_type) <- cluster_annotation$cluster
integrated_data <- RenameIdents(integrated_data, cell_type)
integrated_data$cell_types <- Idents(integrated_data)




##-----------------------------
##  Plotting renamed clusters
##-----------------------------

DimPlot(integrated_data, group.by = "seurat_clusters",order = T, raster = T, label= T)
p3 = DimPlot(integrated_data, label = T, group.by = "cell_types", cols = getPalette(length(levels(integrated_data))), order = T, raster = T,
            label.box = T, repel = T)
ggsave(plot=p3, filename=file.path(dataDir, "Fl_integrated_clusters_annotated.pdf"), units = "in", width=10, height= 7, dpi=600)

##DimPlot(integrated_data, group.by = "DiseaseState", order = T, raster = T)
DimPlot(integrated_data, group.by = "sample_name", order = T, raster = T)

# integrated_data$cell_types <- factor(integrated_data$cell_types,
#                                      levels = c("Acinar","Acinar_2","Ductal",
#                                                 "Endothelial","Fibroblast",
#                                                 "T-cells","B-cells","Myeloid","Granulocytes",
#                                                 "Mast","Cycling","Endocrine"))
# DotPlot(integrated_data,
#         assay = "CorrectedCounts",
#         features = c("SPINK1", "CTRB1", "PRSS1", "AMY2A",
#                      "KRT19", "SOX9", "MUC1", "KRT18","MUC2",
#                      "CDH5", "VWF",
#                      "PDGFRB", "PDGFRA", "PDPN", "ACTA2","TAGLN",
#                      "PTPRC", "CD14", "APOE", "C1QA", "CD68", "HLA-DRA",
#                      "CD3E", "CD4", "FOXP3", "CD8A", "GZMB",
#                      "CD19", "CD79A", 'MS4A1',
#                      "TPSAB1", "CPA3",
#                      "S100A8", "CSF3R","G0S2",
#                      "MKI67","INS"),
#         group.by = "cell_types") +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7, face = "bold"), axis.text.y = element_text(size = 8, face = "bold")) +
#   coord_flip() + ggtitle("Corrected Counts")


#################################################################
##         Plotting cell type abundance in each sample         ##
#################################################################

#splitting samples
samples.list <- unique(integrated_data$sample_name)
clusters <- lapply(samples.list, function(x){
  subset <- subset(integrated_data, subset = sample_name == x)
  dist <- data.frame(table(subset$cell_types))
  
  return(dist)
})

names(clusters) <- samples.list

#calculate relative freq (fractions) of each cell type
clusters_percent <- lapply(clusters, FUN = function(x){
  summ <- sum(x$Freq)
  x$Freq <- (x$Freq/summ)
  return(x)
})

#making things ggplot-friendly!
clusters_dist <- reshape2::melt(clusters, id.var = "Var1")
colnames(clusters_dist) <- c("cell_type","variable","value","sample")
clusters_percent_dist <- reshape2::melt(clusters_percent, id.var = "Var1")
colnames(clusters_percent_dist) <- c("cell_type","variable","value","sample")

#calculating Shannon Entropy as index for diversity
# entropy <- reshape::melt(lapply(clusters, FUN = function(x){
#   Entropy(x$Freq)
# }))
# entropy$value <- scale(entropy$value)

#plotting
p4 = ggplot(clusters_dist, aes(fill=cell_type, y = value, x = sample)) +
  geom_bar(position = "stack", stat = "identity") +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 45, size = 8, face = "bold", hjust = 1))+
  ylab("count") + xlab("sample") + ggtitle("Cell Types Abundance") +
  scale_fill_manual(values = getPalette(length(levels(integrated_data))))
ggsave(plot=p4, filename=file.path(dataDir, "Fl_integrated_clusters_annotated_abundances.pdf"), units = "in", width=10, height= 7, dpi=600)


p5 =ggplot(clusters_percent_dist, aes(fill=cell_type, y = value, x = sample)) +
  geom_bar(position = "stack", stat = "identity") +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 45, size = 8, face = "bold", hjust = 1))+
  ylab("count") + xlab("sample") + ggtitle("Relative Cell Types Abundance") +
  scale_fill_manual(values = getPalette(length(levels(integrated_data))))
ggsave(plot=p5, filename=file.path(dataDir, "Fl_integrated_clusters_annotated_relative_abundances.pdf"), units = "in", width=10, height= 7, dpi=600)

