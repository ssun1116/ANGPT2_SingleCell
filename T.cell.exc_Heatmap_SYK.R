# load packages
library(ggplot2) # ggplot() for plotting
library(dplyr) # data reformatting
library(tidyr) # data reformatting
library(stringr) # string manipulation
library(pheatmap)
library(tidyverse)
library(WGCNA)
#setwd("~/Dropbox/ANGPT2_SingleCell")

#############################################
exc.up = c("FSCN1","TPM2", "H19", "BGN","HMGB1","EIF3H","DCAF13","LYPLA1","MEX3A","SOX4","TIMM50")
exc.down = c("PDE4DIP","PERP","ITGA3","NPC1","NAV2","LOC100126784","CTSD","IRF4","LGALS3","BHLHE41",
             "MMP14","TIMP2","IGSF8","ATP6V0C","IFI27L2","S100A13","S100A6","EMP1","LAMP2")

#############################################
## Exclusion - Total Data

# Total bulk sample
bulk = as.data.frame(read.table("Data/TCGA.SKCM.sampleMap_HiSeqV2", header = T)) # 475 obs.
rownames(bulk) = bulk$sample
bulk$sample <- NULL
bulk2 = scale(t(bulk), center = TRUE, scale = TRUE) %>% data.frame
quantile(bulk2$ANGPT2, c(0, 0.5, 0.9, 1))
bulk2$quant = ifelse(bulk2$ANGPT2 <= -0.07490557, "<50th", 
                     ifelse(bulk2$ANGPT2 >= 1.34707025, ">90th", "50-90th"))
quant_tot.50 = bulk2 %>% filter(bulk2$quant == "<50th") %>% rownames(.) # Sample names included in <50th
quant_tot.50.90 = bulk2 %>% filter(bulk2$quant == "50-90th") %>% rownames(.) # Sample names included in 50-90th
quant_tot.90 = bulk2 %>% filter(bulk2$quant == ">90th") %>% rownames(.) # Sample names included in >90th

# Use ANGPT2 expression level in ordering TCGA samples
bulk2 = bulk2 %>% data.frame %>% arrange(desc(ANGPT2)) 
ang.bulk = rev(rownames(bulk2)) # save sample order in string format

## Jerby Paper style Heatmap
bulk.up = bulk2[,colnames(bulk2) %in% exc.up] 
bulk.up = cbind.data.frame(rownames(bulk.up), bulk.up)
colnames(bulk.up)[1] = "sample"
bulk.up = bulk.up %>% gather(colnames(bulk.up)[2:ncol(bulk.up)], key = 'gene', value = 'exp')

bulk.down = bulk2[,colnames(bulk2) %in% exc.down]
bulk.down = cbind.data.frame(rownames(bulk.down), bulk.down)
colnames(bulk.down)[1] = "sample"
bulk.down = bulk.down %>% gather(colnames(bulk.down)[2:ncol(bulk.down)], key = 'gene', value = 'exp')

## Combine bulk up and bulk down
bulk.tmp = rbind.data.frame(bulk.up, bulk.down)
bulk.tmp2= bulk.tmp %>% as.tibble() %>% mutate_if(is.numeric, ~ifelse(. > 2, 2, .)) %>% mutate_if(is.numeric, ~ifelse(. < -2, -2, .)) %>% as.data.frame

gene.order = c(exc.up, exc.down)
b <- c(-2, -1, 0, 1, 2)
textcol <- "grey40"
ggplot(bulk.tmp2,aes(x= factor(gene, level = gene.order) , y=factor(sample, level = ang.bulk), fill = exp)) +
  geom_tile(show.legend = F) +
  guides(fill=guide_legend(title="Gene Expression\nLevels.")) +
  labs(x="",y="",title="")+
  scale_fill_gradientn(limits = c(-2,2),
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
ggsave("Figures/paper.heatmap_final.total_exc.sig_selected_1214.pdf", last_plot(),width = 6.5, height = 4)


## pHeatmap clustering with selected genes
my_gene_col = as.data.frame(c(exc.up, exc.down))
colnames(my_gene_col) = "gene"
my_gene_col <- my_gene_col %>% mutate(Exclusion = ifelse(test = gene %in% exc.up, yes = "exc_up", no = "exc_down"))
rownames(my_gene_col) <- my_gene_col$gene
my_gene_col$gene <- NULL

my_sample_col <- as.data.frame(rownames(bulk2))
colnames(my_sample_col) = "sample"
my_sample_col <- my_sample_col %>% mutate(Group = ifelse(test = sample %in% quant_tot.50, yes = "< 50th",
                                                         ifelse(test = sample %in% quant_tot.50.90, "50-90th", "> 90th")))
rownames(my_sample_col) <- my_sample_col$sample
my_sample_col$sample <- NULL


## Paper version
bulk.up = bulk2[,colnames(bulk2) %in% exc.up]
bulk.down = bulk2[,colnames(bulk2) %in% exc.down]

## Combine bulk up and bulk down
bulk.tmp = cbind.data.frame(bulk.up, bulk.down)
bulk.tmp2= bulk.tmp %>% as.tibble() %>% mutate_if(is.numeric, ~ifelse(. > 2, 2, .)) %>% mutate_if(is.numeric, ~ifelse(. < -2, -2, .)) %>% as.data.frame
rownames(bulk.tmp2) = rownames(bulk.tmp)
out <- pheatmap(bulk.tmp2, annotation_row = my_sample_col, annotation_col = my_gene_col, cluster_rows=F,
                show_rownames=F, legend_labels = NULL,
                treeheight_row = 80,
                treeheight_col = 80)
ggsave("Figures/pheatmap_final.total_paper.exc.sig_selected_1214.pdf", out,width = 10, height = 6)




#############################################
## Metastasis Data

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


# Use ANGPT2 in TCGA sample ordering 
bulk.met2 = bulk.met2 %>% data.frame %>% arrange(desc(ANGPT2)) 
ang.met.bulk = rev(rownames(bulk.met2)) # save sample order in string format


## Paper version Heatmap
bulk.met.up = bulk.met2[,colnames(bulk.met2) %in% exc.up]
bulk.met.up = cbind.data.frame(rownames(bulk.met.up), bulk.met.up)
colnames(bulk.met.up)[1] = "sample"
bulk.met.up = bulk.met.up %>% gather(colnames(bulk.met.up)[2:ncol(bulk.met.up)], key = 'gene', value = 'exp')

bulk.met.down = bulk.met2[,colnames(bulk.met2) %in% exc.down]
bulk.met.down = cbind.data.frame(rownames(bulk.met.down), bulk.met.down)
colnames(bulk.met.down)[1] = "sample"
bulk.met.down = bulk.met.down %>% gather(colnames(bulk.met.down)[2:ncol(bulk.met.down)], key = 'gene', value = 'exp')

## Combine bulk up and bulk down
bulk.met.tmp = rbind.data.frame(bulk.met.up, bulk.met.down)
bulk.met.tmp2= bulk.met.tmp %>% as.tibble() %>% mutate_if(is.numeric, ~ifelse(. > 2, 2, .)) %>% mutate_if(is.numeric, ~ifelse(. < -2, -2, .)) %>% as.data.frame

gene.order = c(exc.up, exc.down)
b <- c(-2, -1, 0, 1, 2)
textcol <- "grey40"
ggplot(bulk.met.tmp2,aes(x= factor(gene, level = gene.order) , y=factor(sample, level = ang.met.bulk), fill = exp)) +
  geom_tile(show.legend = F) +
  guides(fill=guide_legend(title="Gene Expression\nLevels.")) +
  labs(x="",y="",title="")+
  scale_fill_gradientn(limits = c(-2,2),
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
ggsave("Figures/paper.heatmap_final.met_exc.sig_selected_1214.pdf", last_plot(),width = 6.5, height = 4)


### pHeatmap clustering with paper genes
my_gene_col = as.data.frame(c(exc.up, exc.down))
colnames(my_gene_col) = "gene"
my_gene_col <- my_gene_col %>% mutate(Exclusion = ifelse(test = gene %in% exc.up, yes = "exc_up", no = "exc_down"))
rownames(my_gene_col) <- my_gene_col$gene
my_gene_col$gene <- NULL

my_sample_col <- as.data.frame(rownames(bulk.met2))
colnames(my_sample_col) = "sample"
my_sample_col <- my_sample_col %>% mutate(Group = ifelse(test = sample %in% quant_met.50, yes = "< 50th",
                                                         ifelse(test = sample %in% quant_met.50.90, "50-90th", "> 90th")))
rownames(my_sample_col) <- my_sample_col$sample
my_sample_col$sample <- NULL


## Paper version
bulk.met.up = bulk.met2[,colnames(bulk.met2) %in% exc.up]
bulk.met.down = bulk.met2[,colnames(bulk.met2) %in% exc.down]

## Combine bulk up and bulk down
bulk.met.tmp = cbind.data.frame(bulk.met.up, bulk.met.down)
bulk.met.tmp2= bulk.met.tmp %>% as.tibble() %>% mutate_if(is.numeric, ~ifelse(. > 2, 2, .)) %>% mutate_if(is.numeric, ~ifelse(. < -2, -2, .)) %>% as.data.frame
rownames(bulk.met.tmp2) = rownames(bulk.met.tmp)
out <- pheatmap(bulk.met.tmp2, annotation_row = my_sample_col, annotation_col = my_gene_col, cluster_rows=F,
                show_rownames=F, legend_labels = NULL,
                treeheight_row = 80,
                treeheight_col = 80)
ggsave("Figures/pheatmap_final.met_paper.exc.sig_selected_1214.pdf", out,width = 10, height = 6)




#############################################
## Supplementary. Primary data

bulk.pri <- readRDS("Data/pri_TCGA-SKCM_GDCdata.rds") 
colnames(bulk.pri) = gsub("-", ".", colnames(bulk.pri)) 
colnames(bulk.pri) = substr(colnames(bulk.pri), 1, 15)
bulk.pri = bulk[, colnames(bulk) %in% colnames(bulk.pri)] # 103 obs.
bulk.pri2 = scale(t(bulk.pri), center = TRUE, scale = TRUE) %>% data.frame
quantile(bulk.pri2$ANGPT2, c(0, 0.5, 0.9, 1))
bulk.pri2$quant = ifelse(bulk.pri2$ANGPT2 <= -0.157328, "<50th", 
                         ifelse(bulk.pri2$ANGPT2 >= 1.489399, ">90th", "50-90th"))
quant_pri.50 = bulk.pri2 %>% filter(bulk.pri2$quant == "<50th") %>% rownames(.)
quant_pri.50.90 = bulk.pri2 %>% filter(bulk.pri2$quant == "50-90th") %>% rownames(.)
quant_pri.90 = bulk.pri2 %>% filter(bulk.pri2$quant == ">90th") %>% rownames(.)


# Use ANGPT2 in TCGA sample ordering 
bulk.pri2 = bulk.pri2 %>% data.frame %>% arrange(desc(ANGPT2)) 
ang.pri.bulk = rev(rownames(bulk.pri2)) # save sample order in string format


## Paper version Heatmap
bulk.pri.up = bulk.pri2[,colnames(bulk.pri2) %in% exc.up]
bulk.pri.up = cbind.data.frame(rownames(bulk.pri.up), bulk.pri.up)
colnames(bulk.pri.up)[1] = "sample"
bulk.pri.up = bulk.pri.up %>% gather(colnames(bulk.pri.up)[2:ncol(bulk.pri.up)], key = 'gene', value = 'exp')

bulk.pri.down = bulk.pri2[,colnames(bulk.pri2) %in% exc.down]
bulk.pri.down = cbind.data.frame(rownames(bulk.pri.down), bulk.pri.down)
colnames(bulk.pri.down)[1] = "sample"
bulk.pri.down = bulk.pri.down %>% gather(colnames(bulk.pri.down)[2:ncol(bulk.pri.down)], key = 'gene', value = 'exp')

## Combine bulk up and bulk down
bulk.pri.tmp = rbind.data.frame(bulk.pri.up, bulk.pri.down)
bulk.pri.tmp2= bulk.pri.tmp %>% as.tibble() %>% mutate_if(is.numeric, ~ifelse(. > 2, 2, .)) %>% mutate_if(is.numeric, ~ifelse(. < -2, -2, .)) %>% as.data.frame

gene.order = c(exc.up, exc.down)
b <- c(-2, -1, 0, 1, 2)
textcol <- "grey40"
ggplot(bulk.pri.tmp2,aes(x= factor(gene, level = gene.order) , y=factor(sample, level = ang.pri.bulk), fill = exp)) +
  geom_tile(show.legend = F) +
  guides(fill=guide_legend(title="Gene Expression\nLevels.")) +
  labs(x="",y="",title="")+
  scale_fill_gradientn(limits = c(-2,2),
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
ggsave("Figures/paper.heatmap_final.pri_exc.sig_selected_1214.pdf", last_plot(),width = 6.5, height = 4)


### pHeatmap clustering with paper genes
my_gene_col = as.data.frame(c(exc.up, exc.down))
colnames(my_gene_col) = "gene"
my_gene_col <- my_gene_col %>% mutate(Exclusion = ifelse(test = gene %in% exc.up, yes = "exc_up", no = "exc_down"))
rownames(my_gene_col) <- my_gene_col$gene
my_gene_col$gene <- NULL

my_sample_col <- as.data.frame(rownames(bulk.pri2))
colnames(my_sample_col) = "sample"
my_sample_col <- my_sample_col %>% mutate(Group = ifelse(test = sample %in% quant_pri.50, yes = "< 50th",
                                                         ifelse(test = sample %in% quant_pri.50.90, "50-90th", "> 90th")))
rownames(my_sample_col) <- my_sample_col$sample
my_sample_col$sample <- NULL


## Paper version
bulk.pri.up = bulk.pri2[,colnames(bulk.pri2) %in% exc.up]
bulk.pri.down = bulk.pri2[,colnames(bulk.pri2) %in% exc.down]

## Combine bulk up and bulk down
bulk.pri.tmp = cbind.data.frame(bulk.pri.up, bulk.pri.down)
bulk.pri.tmp2= bulk.pri.tmp %>% as.tibble() %>% mutate_if(is.numeric, ~ifelse(. > 2, 2, .)) %>% mutate_if(is.numeric, ~ifelse(. < -2, -2, .)) %>% as.data.frame
rownames(bulk.pri.tmp2) = rownames(bulk.pri.tmp)
out <- pheatmap(bulk.pri.tmp2, annotation_row = my_sample_col, annotation_col = my_gene_col, cluster_rows=F,
                show_rownames=F, legend_labels = NULL,
                treeheight_row = 80,
                treeheight_col = 80)
ggsave("Figures/pheatmap_final.pri_paper.exc.sig_selected_1214.pdf", out,width = 10, height = 6)




















