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

d = read.csv("Resources/GSE115978/GSE115978_counts.csv.gz") # GEO expression matrix
rownames(d) = d$X
d$X <- NULL

## Create Seurat object from GEO expression matrix
sobj = CreateSeuratObject(counts = d, project = "GSE115978", 
                          min.cells = 10, min.features = 200) # columns : cells/samples, rows : features/genes
sobj[["percent.mt"]] <- PercentageFeatureSet(sobj, pattern = "^mt-")
sobj$log10GenesPerUMI <- log10(sobj$nFeature_RNA) / log10(sobj$nCount_RNA)

## Merge all data to visualize violinplot and compare
VlnPlot(sobj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
save_plot("./vlnplot_alldata_0402.pdf", p, base_height = 6, base_width = 12)

## Pre-processing (QC) and Normalizing Data
helpless <- subset(helpless, subset = nFeature_RNA < 5000 & percent.mt < 5 & log10GenesPerUMI >= 0.8)
helpless <- NormalizeData(helpless, normalization.method = "LogNormalize", scale.factor = 10000)
coghelp <- subset(coghelp, subset = nFeature_RNA < 5000 & percent.mt < 5 & log10GenesPerUMI >= 0.8)
coghelp <- NormalizeData(coghelp, normalization.method = "LogNormalize", scale.factor = 10000)
sephelp <- subset(sephelp, subset = nFeature_RNA < 5000 & percent.mt < 5 & log10GenesPerUMI >= 0.8)
sephelp <- NormalizeData(sephelp, normalization.method = "LogNormalize", scale.factor = 10000)
naive <- subset(naive, subset = nFeature_RNA < 5000 & percent.mt < 5 & log10GenesPerUMI >= 0.8) 
naive <- NormalizeData(naive, normalization.method = "LogNormalize", scale.factor = 10000)





