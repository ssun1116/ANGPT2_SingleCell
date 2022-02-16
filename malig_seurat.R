options(stringsAsFactors = F)
library(dplyr)
library(tidyverse)
library(Seurat)
library(patchwork)
library(ggplot2)
library(cowplot)

d = read.csv("Data/GSE115978/GSE115978_tpm.csv") # GEO expression matrix
rownames(d) = d$X
d$X <- NULL

## Annotation files (Cluster)
mal = read.delim("Resources/tumors.mal_tsne_anno.txt", sep = "\t", header = TRUE)
mal = mal[2:nrow(mal),]
rownames(mal) = mal[,1]
mal[,1] <- NULL

d.mal = d[,rownames(mal)]
d.mal.angpt <- d.mal["ANGPT2",]
d.mal.angpt <- t(d.mal.angpt)
mal = cbind.data.frame(d.mal.angpt, mal)
ggplot(data = mal, aes(x = tumor, y = ANGPT2)) + geom_boxplot()

nonmal = read.delim("Resources/tumors.nonmal_tsne_anno.txt", sep = "\t", header = TRUE)
nonmal = nonmal[2:nrow(nonmal),] 
rownames(nonmal) = nonmal[,1]
nonmal[,1] <- NULL

d.nonmal = d[,rownames(nonmal)]
d.nonmal.angpt <- d.nonmal["ANGPT2",]
d.nonmal.angpt <- t(d.nonmal.angpt)
nonmal = cbind.data.frame(d.nonmal.angpt, nonmal)
ggplot(data = nonmal, aes(x = cell.type, y = ANGPT2)) + geom_boxplot()


################# malignants ##################
## Create Seurat object from GEO expression matrix
mal_sobj = CreateSeuratObject(counts = d.mal, project = "Malignant") # columns : cells/samples, rows : features/genes

VlnPlot(mal_sobj, features = ("ANGPT2"), ncol = 3)

#mal_sobj <- NormalizeData(object = mal_sobj)
mal_sobj <- FindVariableFeatures(object = mal_sobj)
mal_sobj <- ScaleData(object = mal_sobj)
mal_sobj <- RunPCA(object = mal_sobj)
mal_sobj <- FindNeighbors(object = mal_sobj)
mal_sobj <- FindClusters(object = mal_sobj)

# Save current cluster identites in object@meta.data under 'clusterID'
# Run only if Seurat::FindClusters() was executed
# With Seurat 3.X, stashing identity classes can be accomplished with the following:
# mal_sobj[["clusterID"]] <- Idents(object = mal_sobj)

mal_sobj@meta.data$seurat_clusters <- mal$tumor
# Set cell identity to sample identity
mal_sobj <- Seurat::SetIdent(object = mal_sobj, value = "seurat_clusters")
# Find all sample specific marker genes 
markers <- Seurat::FindAllMarkers(object = mal_sobj)
#write_xlsx(markers, "Tables/FindAllMarkers_malignant.scobject_1104.xlsx")

mal_sobj <- RunTSNE(object = mal_sobj)
VlnPlot(mal_sobj, features = "ANGPT2")
RidgePlot(test_sobj, features = "ANGPT2")
DimPlot(test_sobj, reduction = "tsne", group.by = "seurat_clusters", label = TRUE, pt.size = .2)
FeaturePlot(test_sobj, features = "ANGPT2")

################################################################################################

df = mal_sobj@assays$RNA@counts
df2 = as.data.frame(df["ANGPT2",])
colnames(df2) = "ANGPT2"

quantile(df2$ANGPT2, probs = c(0, 0.5, 0.9, 0.99, 1))
# 0%        50%        90%        99%       100% 
# 0.00000000 0.03703073 0.66083737 2.84677069 5.35945089 

d.mal.angpt = as.data.frame(d.mal.angpt)
d.mal.angpt$level <- ifelse(d.mal.angpt$ANGPT2 <= 0.03703073, "low",
                            ifelse(d.mal.angpt$ANGPT2 >= 0.66083737, ifelse(d.mal.angpt$ANGPT2 >= 2.84677069, "very high", "high"),"moderate"))
d.mal.angpt = cbind.data.frame(d.mal.angpt, mal)

################################################################################################

ggplot(d.mal.angpt, aes(x = factor(level, levels = c("low", "moderate", "high", "very high")), y = ANGPT2)) + geom_boxplot()
ggplot(d.mal.angpt, aes(x = factor(level, levels = c("low", "moderate", "high", "very high")), y = as.numeric(OE.Immune_resistance))) + geom_boxplot()

test_sobj = mal_sobj
test_sobj@meta.data$seurat_clusters <- factor(d.mal.angpt$level, levels = c("low", "moderate", "high", "very high"))
# Set cell identity to sample identity
test_sobj <- Seurat::SetIdent(object = test_sobj, value = "seurat_clusters")
# Find all sample specific marker genes 
markers <- Seurat::FindAllMarkers(object = test_sobj)
#write_xlsx(markers, "Tables/FindAllMarkers_malignant.scobject_1104.xlsx")

test_sobj <- RunTSNE(object = test_sobj)
VlnPlot(test_sobj, features = "ANGPT2") + geom_boxplot(fill = "white", width = 0.3)
RidgePlot(test_sobj, features = "ANGPT2")
DimPlot(test_sobj, reduction = "tsne", group.by = "seurat_clusters", label = TRUE, pt.size = .2)
FeaturePlot(test_sobj, features = "ANGPT2")

VlnPlot(test_sobj, features = "ANGPT2") + geom_boxplot(fill = "white", width = 0.3)



################# nonmalignants ##################
## Create Seurat object from GEO expression matrix
nonmal.epi = nonmal %>% filter(cell.type == "Endothelial.cell")
d.nonmal = d[,rownames(nonmal.epi)]
d.nonmal.angpt <- d.nonmal["ANGPT2",]
d.nonmal.angpt <- t(d.nonmal.angpt)
nonmal.epi = cbind.data.frame(d.nonmal.angpt, nonmal.epi)

nonmal.epi_sobj = CreateSeuratObject(counts = d.nonmal, project = "nonmalignant") # columns : cells/samples, rows : features/genes
nonmal.epi_sobj <- NormalizeData(object = nonmal.epi_sobj)
nonmal.epi_sobj <- FindVariableFeatures(object = nonmal.epi_sobj)
nonmal.epi_sobj <- ScaleData(object = nonmal.epi_sobj)
nonmal.epi_sobj <- RunPCA(object = nonmal.epi_sobj)
nonmal.epi_sobj <- FindNeighbors(object = nonmal.epi_sobj)
nonmal.epi_sobj <- FindClusters(object = nonmal.epi_sobj)

# Save current cluster identites in object@meta.data under 'clusterID'
# Run only if Seurat::FindClusters() was executed
# With Seurat 3.X, stashing identity classes can be accomplished with the following:
# nonmal_sobj[["clusterID"]] <- Idents(object = nonmal_sobj)
test_non.epi.sobj = nonmal.epi_sobj
quantile(nonmal.epi$ANGPT2, c(0, 0.5, 0.9, 1))
nonmal.epi$level <- ifelse(nonmal.epi$ANGPT2 <= 0.2630264, "low",
                          ifelse(nonmal.epi$ANGPT2 >= 4.5930595, "high","moderate"))

test_non.epi.sobj@meta.data$seurat_clusters <- nonmal.epi$level
# Set cell identity to sample identity
test_non.epi.sobj <- Seurat::SetIdent(object = test_non.epi.sobj, value = "seurat_clusters")
# Find all sample specific marker genes 
markers <- Seurat::FindAllMarkers(object = test_non.epi.sobj)
write_xlsx(markers, "Tables/FindAllMarkers_epithelial.scobject_1110.xlsx")

test_non.sobj <- RunTSNE(object = test_non.sobj)
VlnPlot(test_non.epi.sobj, features = "ANGPT2") + geom_boxplot(fill = "white", width = 0.3)

DimPlot(test_non.sobj, reduction = "tsne", group.by = "seurat_clusters", label = TRUE, pt.size = .1)
RidgePlot(test_non.sobj, features = "ANGPT2")


marker_mal = read_excel("Tables/FindAllMarkers_malignant.scobject_1104.xlsx")
marker_nonmal = read_excel("Tables/FindAllMarkers_nonmalignant.scobject_1104.xlsx")
marker_mal = marker_mal[marker_mal$cluster=="Mel89",]
marker_nonmal = marker_nonmal[marker_nonmal$cluster=="Endothelial.cell",]




keep_genes = rownames(test_sobj)[rowSums(test_sobj) >= 10]  # sampling하고 나서 유전자 발현량이 0이 되는 유전자 찾기
mal_sobj = subset(test_sobj, features = keep_genes) #feature filtering
print (dim(mal_sobj)) # 14158 5000
mal_sobj = RenameCells(mal_sobj, add.cell.id = "Malignant")

keep_genes = rownames(test_non.sobj)[rowSums(test_non.sobj) >= 10]  # sampling하고 나서 유전자 발현량이 0이 되는 유전자 찾기
nonmal_sobj = subset(test_non.sobj, features = keep_genes) #feature filtering
print (dim(nonmal_sobj)) # 14158 5000
nonmal_sobj = RenameCells(nonmal_sobj, add.cell.id = "Non_Malignant")

d.sc = list()
d.sc[['Malignant']] <- mal_sobj
d.sc[['Non_Malignant']] <- nonmal_sobj



sc <- merge(d.sc[["Malignant"]], d.sc[['Non_Malignant']]) # 14740 15000. Malat1 removed.
sc@meta.data$orig.ident <- factor(sc@meta.data$orig.ident, levels = c("Malignant", "Non_Malignant"))
sc <- NormalizeData(sc, normalization.method = "LogNormalize", scale.factor = 10000)
sc <- FindVariableFeatures(sc, selection.method = "vst", nfeatures = 5000)
## Run ScaleData code and save the RDS data - remove after RunHarmony.
sc <- ScaleData(sc, features = rownames(sc), vars.to.regress = c("nFeature_RNA", "nCount_RNA")) ## nFeature,nCount가 다름.
sc <- RunPCA(sc, features = VariableFeatures(object = sc)) 


## Run Harmony
set.seed(42)
options(repr.plot.height = 3, repr.plot.width = 5)
sc1 <- RunHarmony(sc, "orig.ident", plot_convergence = T) # 4 iterations

## Perplexity : determine clustering strength (higher value -> more distinct clusters)
ndims = 10
sc1 <- RunTSNE(sc1, reduction = "harmony", dims = 1:ndims, min.dist = 0.1, perplexity = 100) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:ndims) %>% 
  FindClusters(resolution = 0.6) %>% 
  identity()

# find markers for every cluster compared to all remaining cells, report only the positive ones
markers <- FindAllMarkers(sc1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#write.table(markers, paste('Tables/table.FindAllMarkers_without_naive', 'Harmony', 'per100', '0406', 'txt', sep='.'), 
#            sep='\t', quote = F, row.names = F, col.names = T)
t = markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC)
#write.table(t, paste('Tables/table.Top15FindAllMarkers_without_naive', 'Harmony', 'per100', '0406', 'txt', sep='.'), 
#            sep='\t', quote = F, row.names = F, col.names = T)


df = sc1@assays$RNA@counts
df2 = as.data.frame(df["ANGPT2",])
colnames(df2) = "ANGPT2"

quantile(df2$ANGPT2, probs = c(0, 0.5, 0.9, 0.99, 1))
# 0%       50%       90%       99%      100% 
# 0.0000000 0.1505597 0.9900835 3.7908301 7.3632499 

df2$level <- ifelse(df2$ANGPT2 <= 0.1505597, "low",
                            ifelse(df2$ANGPT2 >= 0.9900835, "high","moderate"))
#d.mal.angpt = cbind.data.frame(d.mal.angpt, mal)


sc1@meta.data$seurat_clusters <- factor(df2$level, levels = c("low", "moderate", "high"))
# Set cell identity to sample identity
sc1 <- Seurat::SetIdent(object = sc1, value = "seurat_clusters")
# Find all sample specific marker genes 
markers <- Seurat::FindAllMarkers(object = sc1)
write_xlsx(markers, "Tables/FindAllMarkers_harmony_all.scobject_1110.xlsx")

test_sobj <- RunTSNE(object = test_sobj)
VlnPlot(sc1, features = "ANGPT2") + geom_boxplot(fill = "white", width = 0.3)
RidgePlot(sc1, features = "ANGPT2")
DimPlot(test_sobj, reduction = "tsne", group.by = "seurat_clusters", label = TRUE, pt.size = .2)
FeaturePlot(test_sobj, features = "ANGPT2")

VlnPlot(test_sobj, features = "ANGPT2") + geom_boxplot(fill = "white", width = 0.3)

########


d.nonmal = d[,rownames(nonmal)]
d.nonmal.angpt <- d.nonmal["ANGPT2",]
d.nonmal.angpt <- t(d.nonmal.angpt)
nonmal = cbind.data.frame(d.nonmal.angpt, nonmal)

# Save current cluster identites in object@meta.data under 'clusterID'
# Run only if Seurat::FindClusters() was executed
# With Seurat 3.X, stashing identity classes can be accomplished with the following:
# nonmal_sobj[["clusterID"]] <- Idents(object = nonmal_sobj)
test_non.sobj = nonmal_sobj
test_non.sobj@meta.data$seurat_clusters <- nonmal$cell.type
# Set cell identity to sample identity
test_non.sobj <- Seurat::SetIdent(object = test_non.sobj, value = "seurat_clusters")
# Find all sample specific marker genes 
markers <- Seurat::FindAllMarkers(object = test_non.sobj)
write_xlsx(markers, "Tables/FindAllMarkers_nonmalignant.scobject_1104.xlsx")

test_non.sobj <- RunTSNE(object = test_non.sobj)
DimPlot(test_non.sobj, reduction = "tsne", group.by = "seurat_clusters", label = TRUE, pt.size = .1)
RidgePlot(test_non.sobj, features = "ANGPT2")


marker_mal = read_excel("Tables/FindAllMarkers_malignant.scobject_1104.xlsx")
marker_nonmal = read_excel("Tables/FindAllMarkers_nonmalignant.scobject_1104.xlsx")
marker_mal = marker_mal[marker_mal$cluster=="Mel89",]
marker_nonmal = marker_nonmal[marker_nonmal$cluster=="Endothelial.cell",]








