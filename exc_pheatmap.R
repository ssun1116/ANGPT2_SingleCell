# load packages
library(ggplot2) # ggplot() for plotting
library(dplyr) # data reformatting
library(tidyr) # data reformatting
library(stringr) # string manipulation
library(pheatmap)
library(tidyverse)
library(WGCNA)
setwd("~/Dropbox/ANGPT2_SingleCell")

##

final.met = readRDS("Results/FINAL_Tcell_Exclusion_Met_1123.rds")
#final.pri = readRDS("Results/")
final.total = readRDS("Results/FINAL_Tcell_Exclusion_Total_1123.rds")
load("ImmuneResistance/Results/Resistance/Functional/functional.resistance.sig.RData")
results.paper = readRDS("ImmuneResistance/Results/Resistance/Exclusion/FINAL_Tcell_Exclusion.full.rds")

bulk = as.data.frame(read.table("Data/TCGA.SKCM.sampleMap_HiSeqV2", header = T)) # 475 obs.
rownames(bulk) = bulk$sample
bulk$sample <- NULL

##
bulk2 = scale(t(bulk), center = TRUE, scale = TRUE) %>% data.frame
quantile(bulk2$ANGPT2, c(0, 0.5, 0.9, 1))
bulk2$quant = ifelse(bulk2$ANGPT2 <= -0.07490557, "<50th", 
                     ifelse(bulk2$ANGPT2 >= 1.34707025, ">90th", "50-90th"))
quant_tot.50 = bulk2 %>% filter(bulk2$quant == "<50th") %>% rownames(.)
quant_tot.50.90 = bulk2 %>% filter(bulk2$quant == "50-90th") %>% rownames(.)
quant_tot.90 = bulk2 %>% filter(bulk2$quant == ">90th") %>% rownames(.)

bulk.met <- readRDS("Data/met_TCGA-SKCM_GDCdata.rds")
colnames(bulk.met) = gsub("-", ".", colnames(bulk.met)) 
colnames(bulk.met) = substr(colnames(bulk.met), 1, 15)
bulk.met = bulk[, colnames(bulk) %in% colnames(bulk.met)] # 366 obs.
bulk.met2 = scale(t(bulk.met), center = TRUE, scale = TRUE) %>% data.frame
quantile(bulk.met2$ANGPT2, c(0, 0.5, 0.9, 1))
bulk.met2$quant = ifelse(bulk.met2$ANGPT2 <= -0.08155082, "<50th", 
                         ifelse(bulk.met2$ANGPT2 >= 1.29747959, ">90th", "50-90th"))
quant_met.50 = bulk.met2 %>% filter(bulk.met2$quant == "<50th") %>% rownames(.)
quant_met.50.90 = bulk.met2 %>% filter(bulk.met2$quant == "50-90th") %>% rownames(.)
quant_met.90 = bulk.met2 %>% filter(bulk.met2$quant == ">90th") %>% rownames(.)


# Option 1. Use ANGPT2 in TCGA sample ordering 
bulk2 = bulk2 %>% data.frame %>% arrange(desc(ANGPT2)) 
ang.bulk = rownames(bulk2) # save sample order in string format
bulk.met2 = bulk.met2 %>% data.frame %>% arrange(desc(ANGPT2))
ang.bulk.met = rownames(bulk.met2) # save sample order in string format


# Option 2. Use ANGPT2 related genes in TCGA sample ordering
# related.genes = c("PDGFA","FILIP1","LAYN","CD36","LAMA4","LYVE1","TFPI","SEMA3A","MYCT1","DLL4","SDPR","ANGPT2")
# bulk$related.average <- rowMeans(bulk[,related.genes], na.rm=TRUE)
# bulk = bulk %>% arrange(desc(related.average)) 
# ave.bulk = rownames(bulk) # save sample order in string format



##

## Total data heatmap
exc.seed.up = final.total$sig.final$exc.seed.up
exc.seed.down = final.total$sig.final$exc.seed.down
exc.up = final.total$sig.final$exc.up
exc.down = final.total$sig.final$exc.down
paper.up = c("SERPINF1", "ISYNA1", "RPL6", "NOLC1", "RSL1D1", "RPL21", "HNRNPA1", "PFN1", "CDK4", "ILF2", "SOX4", "C17orf76-AS1", "PABPC1", "RUVBL2", "RPS24", "CACYBP", "DARS", "PAICS", "SMARCA4")
paper.down = c("HSPA1A", "C4A", "CTSL1", "LAMP2", "TAPBP", "HLA-A", "HLA-B", "HLA-H", "CD59", "B2M", "HLA-C", "FGFR1", "HLA-G", "CARD16", "CD47", "AHNAK", "CAV2", "CTSD", "TIMP1", "SLC5A3", "CST3", "CD151", "HLA-E", "MIA", "CD58", "CTSB", "S100A6", "EMP1", "HLA-F", "TSC22D3", "KCNN4", "MT2A")

bulk.up = bulk2[,colnames(bulk2) %in% exc.seed.up]
# bulk.up = cbind.data.frame(rownames(bulk.up), bulk.up)
#colnames(bulk.up)[1] = "sample"
# bulk.up = bulk.up %>% gather(colnames(bulk.up)[2:ncol(bulk.up)], key = 'gene', value = 'exp')

bulk.down = bulk2[,colnames(bulk2) %in% exc.seed.down]
# bulk.down = cbind.data.frame(rownames(bulk.down), bulk.down)
#colnames(bulk.down)[1] = "sample"
# bulk.down = bulk.down %>% gather(colnames(bulk.down)[2:ncol(bulk.down)], key = 'gene', value = 'exp')
 
## Paper version
bulk.up = bulk2[,colnames(bulk2) %in% paper.up]
bulk.up = cbind.data.frame(rownames(bulk.up), bulk.up)
colnames(bulk.up)[1] = "sample"
bulk.up = bulk.up %>% gather(colnames(bulk.up)[2:ncol(bulk.up)], key = 'gene', value = 'exp')

bulk.down = bulk2[,colnames(bulk2) %in% paper.down]
bulk.down = cbind.data.frame(rownames(bulk.down), bulk.down)
colnames(bulk.down)[1] = "sample"
bulk.down = bulk.down %>% gather(colnames(bulk.down)[2:ncol(bulk.down)], key = 'gene', value = 'exp')



## Combine bulk up and bulk down
bulk.tmp = rbind.data.frame(bulk.up, bulk.down)
#bulk.tmp = cbind.data.frame(bulk.up, bulk.down)
#bulk.tmp = bulk.tmp[,2:ncol(bulk.tmp)]
bulk.tmp2= bulk.tmp %>% as.tibble() %>% mutate_if(is.numeric, ~ifelse(. > 2, 2, .)) %>% mutate_if(is.numeric, ~ifelse(. < -2, -2, .)) %>% as.data.frame
rownames(bulk.tmp2) = rownames(bulk.tmp)


## Paper figure style
gene.order = c(exc.seed.up, exc.seed.down)
b <- c(-3, -1.5, 0, 1.5, 3)
textcol <- "grey40"
ggplot(bulk.tmp,aes(x= factor(gene, level = gene.order) , y=factor(sample, level = ang.bulk), fill = exp)) +
  geom_tile(show.legend = F) +
  guides(fill=guide_legend(title="Gene Expression\nLevels.")) +
  labs(x="",y="",title="")+
  scale_fill_gradientn(limits = c(-3,3),
                       colours=c("#0B0BF1", "#898BF5", "white", "#EE8D8B", "#EA3627"),breaks=b, labels=format(b)) +
  #  scale_fill_gradient2(low="#3182bd", mid = "white", high="#E64B35", midpoint = 0) +
  theme_bw(base_size=10)+
  theme(axis.text.x=element_text(size=10,
                                 color = "black",
                                 face = "italic",
                                 hjust=1,
                                 angle=90),
       axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        plot.margin=margin(0.7,0.4,0.1,0.2,"cm"),
        plot.title=element_text(colour=textcol,hjust=0,size=14,face="bold"),
        panel.border = element_rect(colour = "black",size = 0.6,fill = NA))

#ggsave("Figures/heatmap_final.total_exc.sig.pdf", last_plot(),width = 6.5, height = 4)

gene.order = c(paper.up, paper.down)
b <- c(-3, -1.5, 0, 1.5, 3)
textcol <- "grey40"
ggplot(bulk.tmp2,aes(x= factor(gene, level = gene.order) , y=factor(sample, level = ang.bulk), fill = exp)) +
  geom_tile(show.legend = F) +
  guides(fill=guide_legend(title="Gene Expression\nLevels.")) +
  labs(x="",y="",title="")+
  scale_fill_gradientn(limits = c(-3,3),
                       colours=c("#0B0BF1", "#898BF5", "white", "#EE8D8B", "#EA3627"),breaks=b, labels=format(b)) +
  #  scale_fill_gradient2(low="#3182bd", mid = "white", high="#E64B35", midpoint = 0) +
  theme_bw(base_size=10)+
  theme(axis.text.x=element_text(size=10,
                                 color = "black",
                                 face = "italic",
                                 hjust=1,
                                 angle=90),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        plot.margin=margin(0.7,0.4,0.1,0.2,"cm"),
        plot.title=element_text(colour=textcol,hjust=0,size=14,face="bold"),
        panel.border = element_rect(colour = "black",size = 0.6,fill = NA))
ggsave("Figures/heatmap_final.total_Paper.exc.sig.pdf", last_plot(),width = 6.5, height = 4)



## Heatmap
heatmap(
  as.matrix(bulk.tmp), Rowv=NA,
  Colv=as.dendrogram(hclust(dist(t(as.matrix(bulk.tmp)))))
)


## pHeatmap (basic)
out <- pheatmap(bulk.tmp, 
                show_rownames=F, cluster_cols=T, cluster_rows=T, scale="row",
                cex=1, clustering_distance_rows="euclidean", cex=1,
                clustering_distance_cols="euclidean", clustering_method="complete", border_color=FALSE)

### pHeatmap clustering with seed genes
my_gene_col = as.data.frame(c(exc.seed.down, exc.seed.up))
colnames(my_gene_col) = "gene"
my_gene_col <- my_gene_col %>% mutate(Exclusion = ifelse(test = gene %in% exc.seed.up, yes = "exc_up", no = "exc_down"))
rownames(my_gene_col) <- my_gene_col$gene
my_gene_col$gene <- NULL

my_sample_col <- as.data.frame(rownames(bulk.tmp))
colnames(my_sample_col) = "sample"
my_sample_col <- my_sample_col %>% mutate(Group = ifelse(test = sample %in% quant_tot.50, yes = "< 50th",
                                                         ifelse(test = sample %in% quant_tot.50.90, "50-90th", "> 90th")))
rownames(my_sample_col) <- my_sample_col$sample
my_sample_col$sample <- NULL

out <- pheatmap(bulk.tmp2, annotation_row = my_sample_col, annotation_col = my_gene_col, cluster_rows=F,
                show_rownames=F, legend_labels = NULL,
                treeheight_row = 80,
                treeheight_col = 80)

### pHeatmap clustering with paper genes
my_gene_col = as.data.frame(c(paper.up, paper.down))
colnames(my_gene_col) = "gene"
my_gene_col <- my_gene_col %>% mutate(Exclusion = ifelse(test = gene %in% paper.up, yes = "exc_up", no = "exc_down"))
rownames(my_gene_col) <- my_gene_col$gene
my_gene_col$gene <- NULL

my_sample_col <- as.data.frame(rownames(bulk2))
colnames(my_sample_col) = "sample"
my_sample_col <- my_sample_col %>% mutate(Group = ifelse(test = sample %in% quant_tot.50, yes = "< 50th",
                                                         ifelse(test = sample %in% quant_tot.50.90, "50-90th", "> 90th")))
rownames(my_sample_col) <- my_sample_col$sample
my_sample_col$sample <- NULL


## Paper version
bulk.up = bulk2[,colnames(bulk2) %in% paper.up]
bulk.down = bulk2[,colnames(bulk2) %in% paper.down]
## Combine bulk up and bulk down
bulk.tmp = cbind.data.frame(bulk.up, bulk.down)
bulk.tmp2= bulk.tmp %>% as.tibble() %>% mutate_if(is.numeric, ~ifelse(. > 2, 2, .)) %>% mutate_if(is.numeric, ~ifelse(. < -2, -2, .)) %>% as.data.frame
rownames(bulk.tmp2) = rownames(bulk.tmp)
out <- pheatmap(bulk.tmp2, annotation_row = my_sample_col, annotation_col = my_gene_col, cluster_rows=F,
                show_rownames=F, legend_labels = NULL,
                treeheight_row = 80,
                treeheight_col = 80)
ggsave("Figures/pheatmap_final.total_paper.exc.sig.pdf", out,width = 10, height = 5)


## pheatmap - Manually select exclusion program genes
bulk.up = bulk[,colnames(bulk) %in% exc.up[57*3+1:length(exc.up)]]
bulk.down = bulk[,colnames(bulk) %in% exc.down[41*3+1:length(exc.down)]]
bulk.tmp2 = cbind.data.frame(bulk.up, bulk.down)

my_gene_col = as.data.frame(c(colnames(bulk.up), colnames(bulk.down)))
colnames(my_gene_col) = "gene"
my_gene_col <- my_gene_col %>% mutate(Exclusion = ifelse(test = gene %in% exc.up, yes = "exc_up", no = "exc_down"))
rownames(my_gene_col) <- my_gene_col$gene
my_gene_col$gene <- NULL

my_sample_col <- as.data.frame(rownames(bulk.tmp2))
colnames(my_sample_col) = "sample"
my_sample_col <- my_sample_col %>% mutate(Group = ifelse(test = sample %in% quant_50, yes = "< 50th",
                                                         ifelse(test = sample %in% quant_50.90, "50-90th", "> 90th")))
rownames(my_sample_col) <- my_sample_col$sample
my_sample_col$sample <- NULL

out <- pheatmap(bulk.tmp2, annotation_row = my_sample_col, annotation_col = my_gene_col, cluster_rows=F,
                show_rownames=F, legend_labels = NULL,
                treeheight_row = 80,
                treeheight_col = 80)

ggsave("Figures/pheatmap(4)_final.total_exc.sig.pdf", out,width = 14, height = 5)


## pheatmap - Pearson correlation

d = limma::normalizeQuantiles(bulk)
a = as.data.frame(WGCNA::cor(t(d)))
a1 = a %>% select(ANGPT2) %>% mutate(ANGPT2 = abs(ANGPT2)) %>% arrange(desc(ANGPT2))
corr.list = a1 %>% filter(ANGPT2 >= 0.2) 
corr.list$group = ifelse(rownames(corr.list) %in% exc.up, "up", 
                         ifelse(rownames(corr.list) %in% exc.down, "down", "none"))


e = bulk[,colnames(d) %in% a3]

bulk.up = bulk2[,colnames(bulk2) %in% rownames(corr.list[corr.list$group == "up",])]
bulk.down = bulk2[,colnames(bulk2) %in% rownames(corr.list[corr.list$group == "down",])]
bulk.tmp = cbind.data.frame(bulk.up, bulk.down)
bulk.tmp2= bulk.tmp %>% as.tibble() %>% mutate_if(is.numeric, ~ifelse(. > 2, 2, .)) %>% mutate_if(is.numeric, ~ifelse(. < -2, -2, .)) %>% as.data.frame
rownames(bulk.tmp2) = rownames(bulk.tmp)

my_gene_col = as.data.frame(c(colnames(bulk.up), colnames(bulk.down)))
colnames(my_gene_col) = "gene"
my_gene_col <- my_gene_col %>% mutate(Exclusion = ifelse(test = gene %in% exc.up, yes = "exc_up", no = "exc_down"))
rownames(my_gene_col) <- my_gene_col$gene
my_gene_col$gene <- NULL

my_sample_col <- as.data.frame(rownames(bulk.tmp))
colnames(my_sample_col) = "sample"
my_sample_col <- my_sample_col %>% mutate(Group = ifelse(test = sample %in% quant_tot.50, yes = "< 50th",
                                                         ifelse(test = sample %in% quant_tot.50.90, "50-90th", "> 90th")))
rownames(my_sample_col) <- my_sample_col$sample
my_sample_col$sample <- NULL

out <- pheatmap(bulk.tmp2, annotation_row = my_sample_col, annotation_col = my_gene_col, cluster_rows=F,
                scale = "column",
                show_rownames=F, legend_labels = NULL,
                treeheight_row = 80,
                treeheight_col = 80)

ggsave("Figures/pheatmap(4)_final.total_exc.sig.pdf", out,width = 14, height = 5)






## Metastasis data heatmap
exc.seed.up = final.met$sig.final$exc.seed.up
exc.seed.down = final.met$sig.final$exc.seed.down
exc.up = final.met$sig.final$exc.up
exc.down = final.met$sig.final$exc.down

bulk.met.up = bulk.met2[,colnames(bulk.met2) %in% exc.seed.up]
#bulk.met.up = cbind.data.frame(rownames(bulk.met.up), bulk.met.up)
#colnames(bulk.met.up)[1] = "sample"
#bulk.met.up = bulk.met.up %>% gather(colnames(bulk.met.up)[2:ncol(bulk.met.up)], key = 'gene', value = 'exp')

bulk.met.down = bulk.met2[,colnames(bulk.met2) %in% exc.seed.down]
#bulk.met.down = cbind.data.frame(rownames(bulk.met.down), bulk.met.down)
#colnames(bulk.met.down)[1] = "sample"
#bulk.met.down = bulk.met.down %>% gather(colnames(bulk.met.down)[2:ncol(bulk.met.down)], key = 'gene', value = 'exp')

## Combine bulk.met up and bulk.met down
#bulk.met.tmp = rbind.data.frame(bulk.met.up, bulk.met.down)
bulk.met.tmp = cbind.data.frame(bulk.met.up, bulk.met.down)
#bulk.met.tmp = bulk.met.tmp[,2:ncol(bulk.met.tmp)]



bulk.met.up = bulk.met2[,colnames(bulk.met2) %in% paper.up]
bulk.met.up = cbind.data.frame(rownames(bulk.met.up), bulk.met.up)
colnames(bulk.met.up)[1] = "sample"
bulk.met.up = bulk.met.up %>% gather(colnames(bulk.met.up)[2:ncol(bulk.met.up)], key = 'gene', value = 'exp')

bulk.met.down = bulk.met2[,colnames(bulk.met2) %in% paper.down]
bulk.met.down = cbind.data.frame(rownames(bulk.met.down), bulk.met.down)
colnames(bulk.met.down)[1] = "sample"
bulk.met.down = bulk.met.down %>% gather(colnames(bulk.met.down)[2:ncol(bulk.met.down)], key = 'gene', value = 'exp')

## Combine bulk.met up and bulk.met down
bulk.met.tmp = rbind.data.frame(bulk.met.up, bulk.met.down)
bulk.met.tmp = cbind.data.frame(bulk.met.up, bulk.met.down)
#bulk.met.tmp = bulk.met.tmp[,2:ncol(bulk.met.tmp)]

gene.order = c(exc.seed.up, exc.seed.down)


## Paper figure style

b <- c(-3, -1.5, 0, 1.5, 3)
textcol <- "grey40"
ggplot(bulk.met.tmp,aes(x= factor(gene, level = gene.order) , y=factor(sample, level = ang.bulk), fill = exp)) +
  geom_tile(show.legend = F) +
  guides(fill=guide_legend(title="Gene Expression\nLevels.")) +
  labs(x="",y="",title="")+
  scale_fill_gradientn(limits = c(-3,3),
                       colours=c("#0B0BF1", "#898BF5", "white", "#EE8D8B", "#EA3627"),breaks=b, labels=format(b)) +
  #  scale_fill_gradient2(low="#3182bd", mid = "white", high="#E64B35", midpoint = 0) +
  theme_bw(base_size=10)+
  theme(axis.text.x=element_text(size=10,
                                 color = "black",
                                 face = "italic",
                                 hjust=1,
                                 angle=90),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        plot.margin=margin(0.7,0.4,0.1,0.2,"cm"),
        plot.title=element_text(colour=textcol,hjust=0,size=14,face="bold"),
        panel.border = element_rect(colour = "black",size = 0.6,fill = NA))

ggsave("Figures/heatmap_final.met_paper.exc.sig.pdf", last_plot(),width = 6.5, height = 4)

## Heatmap
heatmap(
  as.matrix(bulk.met.tmp), Rowv=NA,
  Colv=as.dendrogram(hclust(dist(t(as.matrix(bulk.met.tmp)))))
)

## pHeatmap
out <- pheatmap(bulk.met.tmp2, 
                show_rownames=F, cluster_cols=T, cluster_rows=T, scale="row",
                cex=1, clustering_distance_rows="euclidean", cex=1,
                clustering_distance_cols="euclidean", clustering_method="complete", border_color=FALSE)

### pHeatmap clustering with paper genes
my_gene_col = as.data.frame(c(paper.up, paper.down))
colnames(my_gene_col) = "gene"
my_gene_col <- my_gene_col %>% mutate(Exclusion = ifelse(test = gene %in% paper.up, yes = "exc_up", no = "exc_down"))
rownames(my_gene_col) <- my_gene_col$gene
my_gene_col$gene <- NULL

my_sample_col <- as.data.frame(rownames(bulk.met.tmp2))
colnames(my_sample_col) = "sample"
my_sample_col <- my_sample_col %>% mutate(Group = ifelse(test = sample %in% quant_met.50, yes = "< 50th",
                                                         ifelse(test = sample %in% quant_met.50.90, "50-90th", "> 90th")))
rownames(my_sample_col) <- my_sample_col$sample
my_sample_col$sample <- NULL

out <- pheatmap(bulk.met.tmp2, annotation_row = my_sample_col, annotation_col = my_gene_col, cluster_rows=F,
                show_rownames=F, legend_labels = NULL,
                treeheight_row = 80,
                treeheight_col = 80)
ggsave("Figures/pheatmap_met.total_paper.exc.sig.pdf", out,width = 10, height = 5)



### pHeatmap clustering

d = limma::normalizeQuantiles(bulk.met)
a = WGCNA::cor(t(d)) %>% data.frame()
a1 = a %>% select(ANGPT2) %>% mutate(ANGPT2 = abs(ANGPT2))
a3 = a2 %>% arrange(desc(ANGPT2)) %>% rownames()
a3 = a3[1:100]

library(pheatmap)


e = d[,colnames(d) %in% a3]
e1 = scale(e)
e1 = e1 %>% as.tibble() %>% mutate_if(is.numeric, ~ifelse(. > 2, 2, .)) %>% mutate_if(is.numeric, ~ifelse(. < -2, -2, .)) %>% as.matrix()

pheatmap(e, scale='column')
pheatmap(e1, scale='column')

a4 = a3 %>% mutate(exc = ifelse(rownames(a3) %in% exc.up, "up",
                                ifelse(rownames(a3) %in% exc.down, "down", "none")))
a5 = a4 %>% filter(ANGPT2 >= 0.25 & exc != "none") # up 9, down 17
my_gene_col = as.data.frame(rownames(a5))
colnames(my_gene_col) = "gene"
my_gene_col <- my_gene_col %>% mutate(Exclusion = ifelse(test = gene %in% exc.seed.up, yes = "exc_up", no = "exc_down"))
rownames(my_gene_col) <- my_gene_col$gene
my_gene_col$gene <- NULL

my_gene_col = as.data.frame(c(exc.seed.down, exc.seed.up))
colnames(my_gene_col) = "gene"
my_gene_col <- my_gene_col %>% mutate(Exclusion = ifelse(test = gene %in% exc.seed.up, yes = "exc_up", no = "exc_down"))
rownames(my_gene_col) <- my_gene_col$gene
my_gene_col$gene <- NULL
bulk.met.tmp2= bulk.met.tmp %>% as.tibble() %>% mutate_if(is.numeric, ~ifelse(. > 2, 2, .)) %>% mutate_if(is.numeric, ~ifelse(. < -2, -2, .)) %>% as.data.frame
rownames(bulk.met.tmp2) = rownames(bulk.met.tmp)

my_sample_col <- as.data.frame(rownames(bulk.met.tmp2))
colnames(my_sample_col) = "sample"
my_sample_col <- my_sample_col %>% mutate(Group = ifelse(test = sample %in% quant_met.50, yes = "< 50th",
                                                         ifelse(test = sample %in% quant_met.50.90, "50-90th", "> 90th")))
rownames(my_sample_col) <- my_sample_col$sample
my_sample_col$sample <- NULL

bulk.met.tmp2 = bulk.met.tmp2 %>% arrange(desc(ANGPT2))
out <- pheatmap(bulk.met.tmp2, annotation_row = my_sample_col, annotation_col = my_gene_col, cluster_rows=F,
                show_rownames=F, legend_labels = NULL,
                treeheight_row = 80,
                treeheight_col = 80)

ggsave("Figures/pheatmap(4)_final.total_exc.sig.pdf", out,width = 14, height = 5)



## Manually select exclusion program genes

exc.up = final.met$sig.final$exc.up
exc.down = final.met$sig.final$exc.down

d = limma::normalizeQuantiles(bulk.met)
a = as.data.frame(WGCNA::cor(t(d)))
a1 = a %>% select(ANGPT2) %>% mutate(ANGPT2 = abs(ANGPT2)) %>% arrange(desc(ANGPT2))
corr.list = a1 %>% filter(ANGPT2 >= 0.2) 
corr.list$group = ifelse(rownames(corr.list) %in% exc.up, "up", 
                         ifelse(rownames(corr.list) %in% exc.down, "down", "none"))


bulk.met.up = bulk.met2[,colnames(bulk.met2) %in% rownames(corr.list[corr.list$group == "up",])]
bulk.met.down = bulk.met2[,colnames(bulk.met2) %in% rownames(corr.list[corr.list$group == "down",])]
bulk.met.tmp = cbind.data.frame(bulk.met.up, bulk.met.down)
bulk.met.tmp2= bulk.met.tmp %>% as.tibble() %>% mutate_if(is.numeric, ~ifelse(. > 2, 2, .)) %>% mutate_if(is.numeric, ~ifelse(. < -2, -2, .)) %>% as.data.frame
rownames(bulk.met.tmp2) = rownames(bulk.met.tmp)


my_gene_col = as.data.frame(c(colnames(bulk.met.up), colnames(bulk.met.down)))
colnames(my_gene_col) = "gene"
my_gene_col <- my_gene_col %>% mutate(Exclusion = ifelse(test = gene %in% exc.up, yes = "exc_up", no = "exc_down"))
rownames(my_gene_col) <- my_gene_col$gene
my_gene_col$gene <- NULL

my_sample_col <- as.data.frame(rownames(bulk.met.tmp2))
colnames(my_sample_col) = "sample"
my_sample_col <- my_sample_col %>% mutate(Group = ifelse(test = sample %in% quant_met.50, yes = "< 50th",
                                                         ifelse(test = sample %in% quant_met.50.90, "50-90th", "> 90th")))
rownames(my_sample_col) <- my_sample_col$sample
my_sample_col$sample <- NULL

out <- pheatmap(bulk.met.tmp2, annotation_row = my_sample_col, annotation_col = my_gene_col, cluster_rows=F,
                show_rownames=F, legend_labels = NULL,
                treeheight_row = 80,
                treeheight_col = 80)

ggsave("Figures/pheatmap(4)_final.met_exc.sig.pdf", out,width = 14, height = 5)



## TCGA Metastasis

# ANGPT2 Correlated genes
d = readRDS('bulk.metastasis.RData')
a = WGCNA::cor(d)
a1 = as.data.frame(a)
a2 = a1 %>% select(ANGPT2)
a3 = abs(a2) %>% arrange(desc(ANGPT2))

highest = a3 %>% filter(ANGPT2 >= 0.5)
high = a3 %>% filter(ANGPT2 >= 0.4 & ANGPT2 < 0.5)
medium = a3 %>% filter(ANGPT2 >= 0.3 & ANGPT2 < 0.4)

# bulk data

bulk.met <- readRDS("Data/met_TCGA-SKCM_GDCdata.rds")
colnames(bulk.met) = gsub("-", ".", colnames(bulk.met)) 
colnames(bulk.met) = substr(colnames(bulk.met), 1, 15)
bulk.met = bulk[, colnames(bulk) %in% colnames(bulk.met)] # 366 obs.
bulk.met2 = scale(t(bulk.met), center = TRUE, scale = TRUE) %>% data.frame
quantile(bulk.met2$ANGPT2, c(0, 0.5, 0.9, 1))
bulk.met2$quant = ifelse(bulk.met2$ANGPT2 <= -0.08155082, "<50th", 
                         ifelse(bulk.met2$ANGPT2 >= 1.29747959, ">90th", "50-90th"))
quant_50 = bulk.met2 %>% filter(bulk.met2$quant == "<50th") %>% rownames(.)
quant_50.90 = bulk.met2 %>% filter(bulk.met2$quant == "50-90th") %>% rownames(.)
quant_90 = bulk.met2 %>% filter(bulk.met2$quant == ">90th") %>% rownames(.)
bulk.met = bulk.met %>% arrange(desc(ANGPT2))
bulk.met2 = bulk.met2 %>% arrange(desc(ANGPT2))





















