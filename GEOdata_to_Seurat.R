## Use GEOquery to download data from GEO 
options(stringsAsFactors = F)
library(GEOquery)
library(Biobase)
library(dplyr)
library(tidyverse)
library(Seurat)
library(patchwork)
library(ggplot2)
library(cowplot)
Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 1000)
options(timeout = max(3000, getOption("timeout")))

## The processed scRNA-seq data is provided via the Gene Expression Omnibus (GEO), accession number GEO: GSE115978. 
## The raw scRNA-seq data is being deposited in dbGAP.
data <- getGEO('GSE115978', destdir="./Resources/", GSEMatrix = T) # 0 rows (featureds), 7186 columns (samples). Not accessible.
getGEOSuppFiles('GSE115978') # Download cell.annotations.csv.gz, counts.csv.gz, tpm.csv.gz files
d = read.csv("Data/GSE115978/GSE115978_tpm.csv") # GEO expression matrix
rownames(d) = d$X
d$X <- NULL

## Annotation files (Cluster)
mal = read.delim("Resources/tumors.mal_tsne_anno.txt", sep = "\t", header = TRUE)
mal = mal[2:nrow(mal),]
rownames(mal) = mal[,1]
mal[,1] <- NULL
nonmal = read.delim("Resources/tumors.nonmal_tsne_anno.txt", sep = "\t", header = TRUE)
nonmal = nonmal[2:nrow(nonmal),] 
rownames(nonmal) = nonmal[,1]
nonmal[,1] <- NULL

d.mal = d[,rownames(mal)]
d.nonmal = d[,rownames(nonmal)]

######

d.mal.angpt <- d.mal["ANGPT2",]
d.mal.angpt <- t(d.mal.angpt)
mal = cbind.data.frame(d.mal.angpt, mal)

d.nonmal.angpt <- d.nonmal["ANGPT2",]
d.nonmal.angpt <- t(d.nonmal.angpt)
nonmal = cbind.data.frame(d.nonmal.angpt, nonmal)

ggplot(data = mal, aes(x = tumor, y = ANGPT2)) + geom_boxplot()
ggplot(data = nonmal, aes(x = cell.type, y = ANGPT2)) + geom_boxplot()


################# malignants ##################

## Create Seurat object from GEO expression matrix
mal_sobj = CreateSeuratObject(counts = d.mal, project = "Malignant") # columns : cells/samples, rows : features/genes
non_sobj <- NormalizeData(object = non_sobj)
non_sobj <- FindVariableFeatures(object = non_sobj)
non_sobj <- ScaleData(object = non_sobj)
non_sobj <- RunPCA(object = non_sobj)
non_sobj <- FindNeighbors(object = non_sobj)
non_sobj <- FindClusters(object = non_sobj)
non_sobj <- RunTSNE(object = non_sobj)

non_sobj.dr <- CreateDimReducObject(embeddings = as.matrix(nonmal_embeds), assay = "RNA")
non_sobj@reductions$tsne <- non_sobj.dr
non_sobj@meta.data$cell.type <- nonmal$cell.type

DimPlot(non_sobj, reduction = "tsne", group.by = "cell.type", label = TRUE, pt.size = .1)

VlnPlot(non_sobj, features = angpt.group)
