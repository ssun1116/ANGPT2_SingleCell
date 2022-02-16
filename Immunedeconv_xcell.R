
## Step 1. Download the TCGA-SKCM data with TCGAbiolinks
library(BiocManager)
#BiocManager::install(c("TCGAbiolinks", "DT"))
library(TCGAbiolinks)
library(dplyr)
library(DT)

# query <- GDCquery(project = "TCGA-SKCM",
#                   data.category = "Transcriptome Profiling",
#                   data.type = "Gene Expression Quantification",
#                   workflow.type  = "HTSeq - FPKM",
#                   experimental.strategy = "RNA-Seq",
#                   legacy = FALSE)
# GDCdownload(query, method = "api", files.per.chunk = 10)
# data <- GDCprepare(query)
#saveRDS(data, "total_TCGA-SKCM_GDCdata.rds")
data = readRDS("total_TCGA-SKCM_GDCdata.rds")
data = data %>% filter()

## Step 2. Get expression matrix from TCGA data
#remotes::install_github("icbi-lab/immunedeconv", force = TRUE)
library(SummarizedExperiment)
matrix_data <- as.data.frame(assays(data)$`HTSeq - FPKM`)

## Step 3. Convert rownames from gene id to HGNC gene symbol
option_tx = F
if (option_tx){
  gm = rtracklayer::import('~/Dropbox/Resources/Gencode_hg38_v32/gencode.v32.primary_assembly.annotation.gtf.gz')
  gm1 = as.data.frame(gm) %>% filter(type=='transcript') %>% mutate(tss=ifelse(strand=='+', start, end)) %>% dplyr::select(transcript_id, gene_id, gene_name, tss, gene_type, transcript_type)
  gm1$transcript_id = do.call(rbind.data.frame, strsplit(gm1$transcript_id, '.', fixed = T))[[1]]
  gm1$gene_id = do.call(rbind.data.frame, strsplit(gm1$gene_id, '.', fixed = T))[[1]]
  saveRDS(gm1, '~/Dropbox/Resources/genes.Rdata')
} else{
  gm1 = readRDS('~/Dropbox/Resources/genes.Rdata')
}
gm1$gene_id = do.call(rbind.data.frame, strsplit(gm1$gene_id, '.', fixed = T))[[1]]
gene_coding = gm1 %>% filter(gene_type=='protein_coding') %>% select(gene_name, gene_id) %>% unique
gene_coding = gene_coding[!gene_coding$gene_name %in% c("MATR3", "PINX1", "TMSB15B"),] # duplicates in TCGA data : MATR3, PINX1, TMSB15B

d = as.data.frame(rownames(matrix_data))
colnames(d) = "gene_id" # 56602 obs.
d = merge(d, gene_coding, by = "gene_id") %>% unique

matrix_data = matrix_data %>% filter(rownames(.) %in% d$gene_id)
rownames(matrix_data) <- d$gene_name 


## Step 4.Start deconvolution

### QuanTIseq
#res_quantiseq = immunedeconv::deconvolute(matrix_data, "quantiseq")

### xCell
load(file = "xCell.data.rda")
res_xcell = immunedeconv::deconvolute(matrix_data, "xcell")

## ANGPT2
data_angpt2 <- matrix_data %>% filter(rownames(matrix_data) == "ANGPT2") %>% t %>% data.frame %>% arrange(ANGPT2)
#write.table(data_angpt2, "Tables/total.data_angpt2.ordered.txt")

angpt2.order <- rownames(data_angpt2)
res_xcell_angpt2.ordered = res_xcell[,angpt2.order]                                      
res_xcell_angpt2.ordered = cbind.data.frame(res_xcell$cell_type, res_xcell_angpt2.ordered)
colnames(res_xcell_angpt2.ordered)[1] <- "cell_type"
#write.table(res_xcell_angpt2.ordered, "Tables/total.res_xcell_angpt2.ordered.txt")

library(ggplot2)
library(tidyr)
library(cowplot)
library(RColorBrewer)
my.cols <- brewer.pal(12, "Paired")
my.cols1 <- my.cols[c(2,4,6,8,10,12)]
my.cols2 <- my.cols[c(1, 3, 5, 7, 9, 11)]

met.data = readRDS("met_TCGA-SKCM_GDCdata.rds")
met.matrix <- as.data.frame(assays(met.data)$`HTSeq - FPKM`)
met.list = colnames(met.matrix)
pri.data = readRDS("pri_TCGA-SKCM_GDCdata.rds")
pri.matrix <- as.data.frame(assays(pri.data)$`HTSeq - FPKM`)
pri.list = colnames(pri.matrix)

## ANGPT2
data_angpt2 <- matrix_data %>% filter(rownames(matrix_data) == "ANGPT2") %>% t %>% data.frame %>% arrange(ANGPT2)
#write.table(data_angpt2, "Tables/total.data_angpt2.ordered.txt")
data_angpt2$sample <- rownames(data_angpt2)

library(ggpubr)
for (i in 1:39){
  print(i)
  str <- res_xcell_angpt2.ordered$cell_type[i]
  d <- res_xcell_angpt2.ordered[i,] %>%
    gather(sample, score, -cell_type)
  d$type = ifelse(d$sample %in% met.list, "met", 
                  ifelse(d$sample %in% pri.list, "pri", NA))
  d = d[!is.na(d$type),] # 470 obs.
  samplelist <- unlist(d$sample)
  data_angpt2 = data_angpt2 %>% filter(data_angpt2$sample %in% samplelist)
  d = merge(d, data_angpt2, by = "sample")
  d = d %>% arrange(ANGPT2)
  
  p <- d %>% ggplot(aes(x= log2(ANGPT2), y=score, linetype= type, group = type)) +
    geom_point(size=1.5, aes(shape = type, color = type)) + 
    scale_shape_manual(values=c(16, 4))+
    scale_color_manual(values=c(my.cols2[i%%6 + 1],my.cols1[i%%6 + 1]))+
    scale_linetype_manual(values = c("solid", "dashed")) +
    # geom_smooth(aes(group = 1), size = 2,  color = my.cols[i%%6 + 1]) + ## loess function
    stat_smooth(method = 'lm', aes(group = type), size = 2, show.legend = F) + ## linear regression
    #    geom_vline(xintercept= c(247,425), color = "grey1", linetype="dashed", size = 0.5, alpha = 0.5) + #met
    #    geom_vline(xintercept= c(253,426), color = "grey1", linetype="dotted", size = 0.5, alpha = 0.5) +
    #    annotate("text", x = 130, y =  max(d$score), label = "\n<50th", size = 5, color = "grey1" ) +
    #    annotate("text", x = 130, y =  max(d$score), label = "\n\n\n(met.n = 179, pri.n = 51)", 
    #             size = 4.5, color = "grey1", alpha = 0.8, fontface = "italic") +
    #    annotate("text", x = 420, y = max(d$score) , label = "\n>90th", size = 5, color = "grey1" ) +
    #    annotate("text", x = 420, y = max(d$score), label = "\n\n\n\n(met.n = 36, \npri.n = 11)", 
    #             size = 4.5, color = "grey1", alpha = 0.8, fontface = "italic") +
    #scale_color_brewer(palette="Paired", guide= "none") +
    theme_bw() +
    theme(panel.background = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_text(color = "gray10", angle = 90, vjust = 0.5, hjust=1, size = 10),
          axis.title.x = element_text(vjust = -0.1),
          panel.grid.major.y = element_line(size = 0.4, linetype = "solid",color = "gray90"),
          panel.grid.major.x = element_blank(),
          legend.title = element_text(color = "gray10", size = 13),
          legend.text = element_text(color = "gray10", size = 13),
          axis.line = element_line(size = 0.5, linetype = "solid", color = "black"),
          axis.ticks = element_blank(),
          plot.title = element_text(color = "Black", size = 15),
          axis.title = element_text(color = "gray10", size = 14),
          #         legend.background = element_rect(fill="lightgrey", size=0.5, linetype="solid"),
          legend.position="bottom") +
    labs(x = expression(log[2]("ANGPT2 Expr.")),
         y = 'xCell score',
         title = str)
  p
  save_plot(paste("Figures/xcell_lm/lm_by.sample.", str, ".0920.png", sep = ""), p, base_height = 5.5, base_width = 7.5)
}

for (i in 1:39){
  print(i)
  str <- res_xcell_angpt2.ordered$cell_type[i]
  d <- res_xcell_angpt2.ordered[i,] %>%
    gather(sample, score, -cell_type)
  d$type = ifelse(d$sample %in% met.list, "met", 
                  ifelse(d$sample %in% pri.list, "pri", NA))
  d = d[!is.na(d$type),] # 470 obs.
  samplelist <- unlist(d$sample)
  data_angpt2 = data_angpt2 %>% filter(data_angpt2$sample %in% samplelist)
  d = merge(d, data_angpt2, by = "sample")
  d = d %>% arrange(ANGPT2)
  
  print(ggplot(d, aes(x= log2(ANGPT2), y=score)) +
          stat_smooth(method = 'lm', aes(linetype = type), color = "gray30", fill = "grey70", size = 3, alpha = 0.5, show.legend = T) + ## linear regression
          geom_point(aes(shape = type), color = "grey25", size=1.1, stroke =0.55, alpha = 0.7) + # Make dots bigger
          scale_shape_manual(values=c(5, 4))+
          geom_vline(xintercept= c(-0.5407935, 2.739935), color = c("darkred", "darkblue"), linetype="dashed", size = 1.75, alpha = 0.65) + #met
          geom_vline(xintercept= c(0.2470981, 2.800226), color = c("darkred", "darkblue"), linetype="solid", size = 1.75, alpha = 0.65) +
          annotate("text", x = -2.75, y = max(d$score), label = "\n<50th", size = 5, color = "grey1" ) +
          annotate("text", x = -2.75, y = 0.99*max(d$score), label = "\n\n\n\n(met.n = 179)\n(pri.n = 51)",
                   size = 4, color = "grey1", alpha = 0.8, fontface = "italic") +
          annotate("text", x = 5.5, y = max(d$score) , label = "\n>90th", size = 5, color = "grey1" ) +
          annotate("text", x = 5.5, y = 0.99*max(d$score), label = "\n\n\n\n(met.n = 36)\n(pri.n = 11)",
                   size = 4, color = "grey1", alpha = 0.8, fontface = "italic") +
          theme_bw() +
          theme(panel.background = element_blank(),
                axis.text.x = element_blank(),
                axis.text.y = element_text(color = "gray1", angle = 90, vjust = 0.5, hjust=1, size = 10),
                axis.title.x = element_text(vjust = -0.1, size = 13),
                panel.grid.major.y = element_line(size = 0.4, linetype = "solid",color = "gray90"),
                panel.grid.major.x = element_blank(),
                legend.title = element_text(color = "gray10", size = 13),
                legend.text = element_text(color = "gray10", size = 13),
                axis.line = element_line(size = 0.5, linetype = "solid", color = "black"),
                axis.ticks = element_blank(),
                plot.title = element_text(color = "Black", size = 16, face = "bold"),
                axis.title = element_text(color = "gray10", size = 14),
                legend.position="bottom") +
          labs(x = expression(log[2]("ANGPT2 Expr.")),
               y = 'xCell score',
               title = str) +
          guides(shape=guide_legend(override.aes=list(size=1), keywidth = 2, keyheight = 1))) 
  Sys.sleep(3)
  
  grid.ls(grid.force())    # To get the names of all the grobs in the ggplot
  # The edit - to set the size of the point in the legend to 4 mm
  grid.gedit("key-1-3-2.2-4-2-4", size = unit(3, "mm")) 
  grid.gedit("key-1-7-2.2-8-2-8", size = unit(3, "mm")) 
  g <- grid.grab()
  ggsave(plot = g, paste("Figures/xcell_lm/lm_by.sample.", str, ".0924.pdf", sep = ""), width = 12, height = 9)
  #save_plot(paste("Figures/xcell_lm/lm_by.sample.", str, ".0924.pdf", sep = ""), last_plot(), base_height = 4.5, base_width = 6)
}



#################################


## Step 1. Download the TCGA-SKCM data with TCGAbiolinks
library(BiocManager)
#BiocManager::install(c("TCGAbiolinks", "DT"))
library(TCGAbiolinks)
library(dplyr)
library(DT)

met.query <- GDCquery(project = "TCGA-SKCM",
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification",
                      workflow.type  = "HTSeq - FPKM",
                      experimental.strategy = "RNA-Seq",
                      sample.type = c("Metastatic"),
                      legacy = FALSE)

pri.query <- GDCquery(project = "TCGA-SKCM",
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification",
                      workflow.type  = "HTSeq - FPKM",
                      experimental.strategy = "RNA-Seq",
                      sample.type = c("Primary Tumor"),
                      legacy = FALSE)

GDCdownload(met.query, method = "api", files.per.chunk = 10)
met.data <- GDCprepare(met.query)
#saveRDS(met.data, "met_TCGA-SKCM_GDCdata.rds")
met.data = readRDS("met_TCGA-SKCM_GDCdata.rds")
GDCdownload(pri.query, method = "api", files.per.chunk = 10)
pri.data <- GDCprepare(pri.query)
#saveRDS(pri.data, "pri_TCGA-SKCM_GDCdata.rds")
pri.data = readRDS("pri_TCGA-SKCM_GDCdata.rds")

l <- list()
l[[1]] <- met.data
l[[2]] <- pri.data

for (j in 1:2){
  data <- l[[j]]
  ## Step 2. Get expression matrix from TCGA data
  #remotes::install_github("icbi-lab/immunedeconv", force = TRUE)
  library(SummarizedExperiment)
  matrix_data <- as.data.frame(assays(data)$`HTSeq - FPKM`)
  
  ## Step 3. Convert rownames from gene id to HGNC gene symbol
  option_tx = F
  if (option_tx){
    gm = rtracklayer::import('~/Dropbox/Resources/Gencode_hg38_v32/gencode.v32.primary_assembly.annotation.gtf.gz')
    gm1 = as.data.frame(gm) %>% filter(type=='transcript') %>% mutate(tss=ifelse(strand=='+', start, end)) %>% dplyr::select(transcript_id, gene_id, gene_name, tss, gene_type, transcript_type)
    gm1$transcript_id = do.call(rbind.data.frame, strsplit(gm1$transcript_id, '.', fixed = T))[[1]]
    gm1$gene_id = do.call(rbind.data.frame, strsplit(gm1$gene_id, '.', fixed = T))[[1]]
    saveRDS(gm1, '~/Dropbox/Resources/genes.Rdata')
  } else{
    gm1 = readRDS('~/Dropbox/Resources/genes.Rdata')
  }
  gm1$gene_id = do.call(rbind.data.frame, strsplit(gm1$gene_id, '.', fixed = T))[[1]]
  gene_coding = gm1 %>% filter(gene_type=='protein_coding') %>% select(gene_name, gene_id) %>% unique
  gene_coding = gene_coding[!gene_coding$gene_name %in% c("MATR3", "PINX1", "TMSB15B"),] # duplicates in TCGA data : MATR3, PINX1, TMSB15B
  
  d = as.data.frame(rownames(matrix_data))
  colnames(d) = "gene_id" # 56602 obs.
  d = merge(d, gene_coding, by = "gene_id") %>% unique
  
  matrix_data = matrix_data %>% filter(rownames(.) %in% d$gene_id)
  rownames(matrix_data) <- d$gene_name 
  
  ## Step 4.Start deconvolution
  
  ### QuanTIseq
  #res_quantiseq = immunedeconv::deconvolute(matrix_data, "quantiseq")
  ## xCell
  load(file = "xCell.data.rda")
  res_xcell = immunedeconv::deconvolute(matrix_data, "xcell")
  
  ## ANGPT2
  data_angpt2 <- matrix_data %>% filter(rownames(matrix_data) == "ANGPT2") %>% t %>% data.frame %>% arrange(ANGPT2)
  type = ifelse(j == 1, "met", "pri")
  # write.table(data_angpt2, paste("Tables/", type, ".data_angpt2.ordered.txt", sep = ""))
  
  angpt2.order <- rownames(data_angpt2)
  res_xcell_angpt2.ordered = res_xcell[,angpt2.order]
  res_xcell_angpt2.ordered = cbind.data.frame(res_xcell$cell_type, res_xcell_angpt2.ordered)
  colnames(res_xcell_angpt2.ordered)[1] <- "cell_type"
  type = ifelse(j == 1, "met", "pri")
  # write.table(res_xcell_angpt2.ordered, paste("Tables/", type, ".res_xcell_angpt2.ordered.txt", sep = ""))
  
  library(ggplot2)
  library(tidyr)
  library(cowplot)
  library(RColorBrewer)
  my.cols <- brewer.pal(12, "Paired")
  my.cols1 <- my.cols[c(2,4,6,8,10,12)]
  my.cols2 <- my.cols[c(1, 3, 5, 7, 9, 11)]
  
  for (i in 1:39){
    print(i)
    p <- res_xcell_angpt2.ordered[i,] %>%
      gather(sample, score, -cell_type) %>%
      ggplot(aes(x= factor(sample, level = angpt2.order), y=score), color=cell_type) +
      geom_point(size=2, color = my.cols2[i%%6 + 1]) + 
      #     geom_vline(xintercept= ifelse(j == 1, 179, 51), linetype=4, colour = "blue", size = 1, alpha = 0.7 ) +
      geom_vline(xintercept= ifelse(j == 1, 179, 51), color = "grey1", linetype="dashed", size = 0.5, alpha = 0.7) +
      geom_vline(xintercept= ifelse(j == 1, 331, 92), color = "grey1", linetype="dashed", size = 0.5, alpha = 0.7) +
      annotate("text", x = 14, y = 77, label = "<50th (n = 179)", size = 3, color = "gray10" ) +
      # geom_smooth(aes(group = 1), size = 2,  color = my.cols[i%%6 + 1]) + ## loess function
      stat_smooth(method = 'lm', aes(group = 1), size = 2,  color = my.cols1[i%%6 + 1]) + ## linear regression
      facet_wrap(~cell_type, scales="free_x", ncol=3) +
      scale_color_brewer(palette="Paired", guide= "none") +
      #     coord_flip() +
      theme_bw() +
      theme(axis.text.x = element_blank(),
            axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10),
            #panel.border = element_blank(),
            panel.grid.major = element_blank(),
            strip.text = element_text(size=14),
            axis.ticks = element_blank(),
            axis.title = element_text(size = 12)) +
      labs(x = 'ANGPT2 Expression level',
           y = 'xCell score')
    
    str <- res_xcell_angpt2.ordered$cell_type[i]
    
    save_plot(paste(type, str, "0830.png", sep = "."), p, base_height = 5, base_width = 5)
  }
}


## 0920.
l <- list()
l[[1]] <- met.data
l[[2]] <- pri.data

for (j in 1:2){
  data <- l[[j]]
  ## Step 2. Get expression matrix from TCGA data
  #remotes::install_github("icbi-lab/immunedeconv", force = TRUE)
  library(SummarizedExperiment)
  matrix_data <- as.data.frame(assays(data)$`HTSeq - FPKM`)
  
  ## Step 3. Convert rownames from gene id to HGNC gene symbol
  option_tx = F
  if (option_tx){
    gm = rtracklayer::import('~/Dropbox/Resources/Gencode_hg38_v32/gencode.v32.primary_assembly.annotation.gtf.gz')
    gm1 = as.data.frame(gm) %>% filter(type=='transcript') %>% mutate(tss=ifelse(strand=='+', start, end)) %>% dplyr::select(transcript_id, gene_id, gene_name, tss, gene_type, transcript_type)
    gm1$transcript_id = do.call(rbind.data.frame, strsplit(gm1$transcript_id, '.', fixed = T))[[1]]
    gm1$gene_id = do.call(rbind.data.frame, strsplit(gm1$gene_id, '.', fixed = T))[[1]]
    saveRDS(gm1, '~/Dropbox/Resources/genes.Rdata')
  } else{
    gm1 = readRDS('~/Dropbox/Resources/genes.Rdata')
  }
  gm1$gene_id = do.call(rbind.data.frame, strsplit(gm1$gene_id, '.', fixed = T))[[1]]
  gene_coding = gm1 %>% filter(gene_type=='protein_coding') %>% select(gene_name, gene_id) %>% unique
  gene_coding = gene_coding[!gene_coding$gene_name %in% c("MATR3", "PINX1", "TMSB15B"),] # duplicates in TCGA data : MATR3, PINX1, TMSB15B
  
  d = as.data.frame(rownames(matrix_data))
  colnames(d) = "gene_id" # 56602 obs.
  d = merge(d, gene_coding, by = "gene_id") %>% unique
  
  matrix_data = matrix_data %>% filter(rownames(.) %in% d$gene_id)
  rownames(matrix_data) <- d$gene_name 
  
  ## Step 4.Start deconvolution
  
  ### QuanTIseq
  #res_quantiseq = immunedeconv::deconvolute(matrix_data, "quantiseq")
  ## xCell
  load(file = "xCell.data.rda")
  res_xcell = immunedeconv::deconvolute(matrix_data, "xcell")
  res_xcell.t = t(res_xcell)
  colnames(res_xcell.t) <- res_xcell.t[1,]
  res_xcell.t = res_xcell.t[2:nrow(res_xcell.t),]
  
  ## ANGPT2
  data_angpt2 <- matrix_data %>% filter(rownames(matrix_data) == "ANGPT2") %>% t %>% data.frame %>% arrange(ANGPT2)
  type = ifelse(j == 1, "met", "pri")
  # write.table(data_angpt2, paste("Tables/", type, ".data_angpt2.ordered.txt", sep = ""))
  
  res_xcell.t_ordered = cbind.data.frame(data_angpt2, res_xcell.t) %>% arrange(ANGPT2)
  
  # angpt2.order <- rownames(data_angpt2)
  # res_xcell_angpt2.ordered = res_xcell[,angpt2.order]
  # res_xcell_angpt2.ordered = cbind.data.frame(res_xcell$cell_type, res_xcell_angpt2.ordered)
  # colnames(res_xcell_angpt2.ordered)[1] <- "cell_type"
  # type = ifelse(j == 1, "met", "pri")
  # write.table(res_xcell_angpt2.ordered, paste("Tables/", type, ".res_xcell_angpt2.ordered.txt", sep = ""))
  
  library(ggplot2)
  library(tidyr)
  library(cowplot)
  library(RColorBrewer)
  my.cols <- brewer.pal(12, "Paired")
  my.cols1 <- my.cols[c(2,4,6,8,10,12)]
  my.cols2 <- my.cols[c(1, 3, 5, 7, 9, 11)]
  
  d = res_xcell.t_ordered[,1:2] 
  d %>% ggplot(aes(x= log2(ANGPT2), y=d[,2])) +
    geom_point(size=2, color = my.cols2[i%%6 + 1]) + 
    #     geom_vline(xintercept= ifelse(j == 1, 179, 51), linetype=4, colour = "blue", size = 1, alpha = 0.7 ) +
    #    geom_vline(xintercept= ifelse(j == 1, 179, 51), color = "grey1", linetype="dashed", size = 0.5, alpha = 0.7) +
    #    geom_vline(xintercept= ifelse(j == 1, 331, 92), color = "grey1", linetype="dashed", size = 0.5, alpha = 0.7) +
    #    annotate("text", x = 14, y = 77, label = "<50th (n = 179)", size = 3, color = "gray10" ) +
    # geom_smooth(aes(group = 1), size = 2,  color = my.cols[i%%6 + 1]) + ## loess function
    stat_smooth(method = 'lm', aes(group = 1), size = 2,  color = my.cols1[i%%6 + 1]) + ## linear regression
    #    facet_wrap(~cell_type, scales="free_x", ncol=3) +
    scale_color_brewer(palette="Paired", guide= "none") +
    #     coord_flip() +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10),
          #panel.border = element_blank(),
          panel.grid.major = element_blank(),
          strip.text = element_text(size=14),
          axis.ticks = element_blank(),
          axis.title = element_text(size = 12)) +
    labs(x = 'ANGPT2 Expression level',
         y = 'xCell score')
  
  
  
  for (i in 1:39){
    print(i)
    p <- res_xcell_angpt2.ordered[i,] %>%
      gather(sample, score, -cell_type) %>%
      ggplot(aes(x= factor(sample, level = angpt2.order), y=score), color=cell_type) +
      geom_point(size=2, color = my.cols2[i%%6 + 1]) + 
      #     geom_vline(xintercept= ifelse(j == 1, 179, 51), linetype=4, colour = "blue", size = 1, alpha = 0.7 ) +
      geom_vline(xintercept= ifelse(j == 1, 179, 51), color = "grey1", linetype="dashed", size = 0.5, alpha = 0.7) +
      geom_vline(xintercept= ifelse(j == 1, 331, 92), color = "grey1", linetype="dashed", size = 0.5, alpha = 0.7) +
      annotate("text", x = 14, y = 77, label = "<50th (n = 179)", size = 3, color = "gray10" ) +
      # geom_smooth(aes(group = 1), size = 2,  color = my.cols[i%%6 + 1]) + ## loess function
      stat_smooth(method = 'lm', aes(group = 1), size = 2,  color = my.cols1[i%%6 + 1]) + ## linear regression
      facet_wrap(~cell_type, scales="free_x", ncol=3) +
      scale_color_brewer(palette="Paired", guide= "none") +
      #     coord_flip() +
      theme_bw() +
      theme(axis.text.x = element_blank(),
            axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10),
            #panel.border = element_blank(),
            panel.grid.major = element_blank(),
            strip.text = element_text(size=14),
            axis.ticks = element_blank(),
            axis.title = element_text(size = 12)) +
      labs(x = 'ANGPT2 Expression level',
           y = 'xCell score')
    
    str <- res_xcell_angpt2.ordered$cell_type[i]
    
    save_plot(paste(type, str, "0830.png", sep = "."), p, base_height = 5, base_width = 5)
  }
}

### lm
library(broom)
options(scipen = 999)

# total
total.res = read.table("Tables/total.res_xcell_angpt2.ordered.txt")
total.res = data.frame(t(total.res))
colnames(total.res) <- total.res[1,]
total.res = total.res[2:nrow(total.res),]
total.data = read.table("Tables/total.data_angpt2.ordered.txt")
total = cbind.data.frame(total.data, total.res)

str <- colnames(total)
str2 = gsub(" ", "_", str)
str2 = gsub("\\+", "", str2)
str2 = gsub("-", "_", str2)
colnames(total) <- str2
colnames(total)[6] = "T_cell_CD4_non_regulatory"
colnames(total)[37] ="T_cell_regs"


total[,2:ncol(total)] <- sapply(total[,2:ncol(total)], as.numeric)

lmfit <- lm(as.formula(paste(colnames(total)[1], "~", 
                             paste(colnames(total)[2:ncol(total)], collapse = " + "), sep = "")), data=total)


## met
met.res = read.table("Tables/met.res_xcell_angpt2.ordered.txt")
met.res = data.frame(t(met.res))
colnames(met.res) <- met.res[1,]
met.res = met.res[2:nrow(met.res),]
met.data = read.table("Tables/met.data_angpt2.ordered.txt")
met = cbind.data.frame(met.data, met.res)

str <- colnames(met)
str2 = gsub(" ", "_", str)
str2 = gsub("\\+", "", str2)
str2 = gsub("-", "_", str2)
colnames(met) <- str2
colnames(met)[6] = "T_cell_CD4_non_regulatory"
colnames(met)[37] ="T_cell_regs"
met[,2:ncol(met)] <- sapply(met[,2:ncol(met)], as.numeric)

met = met[,1:37] # Exclude score column

met.lmfit <- lm(as.formula(paste(colnames(met)[1], "~", 
                                 paste(colnames(met)[2:ncol(met)], collapse = " + "), sep = "")), data=met)
met.lmfit

met = met[,c(1, 2, 5, 7:9, 11:21, 23:37)] # individual samples
met.glmfit <- glm(as.formula(paste(colnames(met)[1], "~", 
                                   paste(colnames(met)[2:ncol(met)], collapse = " + "), sep = "")), data=met)
met.glmfit
tidy(met.glmfit)

## pri
pri.res = read.table("Tables/pri.res_xcell_angpt2.ordered.txt")
pri.res = data.frame(t(pri.res))
colnames(pri.res) <- pri.res[1,]
pri.res = pri.res[2:nrow(pri.res),]
pri.data = read.table("Tables/pri.data_angpt2.ordered.txt")
pri = cbind.data.frame(pri.data, pri.res)

str <- colnames(pri)
str2 = gsub(" ", "_", str)
str2 = gsub("\\+", "", str2)
str2 = gsub("-", "_", str2)
colnames(pri) <- str2
colnames(pri)[6] = "T_cell_CD4_non_regulatory"
colnames(pri)[37] ="T_cell_regs"
pri[,2:ncol(pri)] <- sapply(pri[,2:ncol(pri)], as.numeric)

pri = pri[,1:37] # Exclude score column

pri.lmfit <- lm(as.formula(paste(colnames(pri)[1], "~", 
                                 paste(colnames(pri)[2:ncol(pri)], collapse = " + "), sep = "")), data=pri)
pri.lmfit

pri = pri[,c(1, 2, 5, 7:9, 11:21, 23:37)] # Exclude score column
pri.lmfit <- lm(as.formula(paste(colnames(pri)[1], "~", 
                                 paste(colnames(pri)[2:ncol(pri)], collapse = " + "), sep = "")), data=pri)
pri.lmfit

## TCGA-SKCM


## Step 1. Download the TCGA-SKCM data with TCGAbiolinks
library(BiocManager)
#BiocManager::install(c("TCGAbiolinks", "DT"))
library(TCGAbiolinks)
library(dplyr)
library(DT)

query <- GDCquery(project = "TCGA-SKCM",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type  = "HTSeq - FPKM",
                  experimental.strategy = "RNA-Seq",
                  legacy = FALSE)
GDCdownload(query, method = "api", files.per.chunk = 10)
data <- GDCprepare(query)
#saveRDS(data, "total_TCGA-SKCM_GDCdata.rds")
data = readRDS("total_TCGA-SKCM_GDCdata.rds")


## Step 2. Get expression matrix from TCGA data
#remotes::install_github("icbi-lab/immunedeconv", force = TRUE)
library(SummarizedExperiment)
matrix_data <- as.data.frame(assays(data)$`HTSeq - FPKM`)


## Step 3. Convert rownames from gene id to HGNC gene symbol
option_tx = F
if (option_tx){
  gm = rtracklayer::import('~/Dropbox/Resources/Gencode_hg38_v32/gencode.v32.primary_assembly.annotation.gtf.gz')
  gm1 = as.data.frame(gm) %>% filter(type=='transcript') %>% mutate(tss=ifelse(strand=='+', start, end)) %>% dplyr::select(transcript_id, gene_id, gene_name, tss, gene_type, transcript_type)
  gm1$transcript_id = do.call(rbind.data.frame, strsplit(gm1$transcript_id, '.', fixed = T))[[1]]
  gm1$gene_id = do.call(rbind.data.frame, strsplit(gm1$gene_id, '.', fixed = T))[[1]]
  saveRDS(gm1, '~/Dropbox/Resources/genes.Rdata')
} else{
  gm1 = readRDS('~/Dropbox/Resources/genes.Rdata')
}
gm1$gene_id = do.call(rbind.data.frame, strsplit(gm1$gene_id, '.', fixed = T))[[1]]
gene_coding = gm1 %>% filter(gene_type=='protein_coding') %>% select(gene_name, gene_id) %>% unique
gene_coding = gene_coding[!gene_coding$gene_name %in% c("MATR3", "PINX1", "TMSB15B"),] # duplicates in TCGA data : MATR3, PINX1, TMSB15B

d = as.data.frame(rownames(matrix_data))
colnames(d) = "gene_id" # 56602 obs.
d = merge(d, gene_coding, by = "gene_id") %>% unique

matrix_data = matrix_data %>% filter(rownames(.) %in% d$gene_id)
rownames(matrix_data) <- d$gene_name 


## Step 4.Start deconvolution

# remotes::install_github("icbi-lab/immunedeconv", force = TRUE)
### QuanTIseq
#res_quantiseq = immunedeconv::deconvolute(matrix_data, "quantiseq")
### xCell
load(file = "Data/xCell.data.rda")
res_xcell = immunedeconv::deconvolute(matrix_data, "xcell")


## Step 5.Order deconvolution results by ANGPT2 exp.
data_angpt2 <- matrix_data %>% filter(rownames(matrix_data) == "ANGPT2") %>% t %>% data.frame %>% arrange(ANGPT2)
#write.table(data_angpt2, "Tables/total.data_angpt2.ordered.txt")
data_angpt2$sample <- rownames(data_angpt2)
angpt2.order <- rownames(data_angpt2)

res_xcell_angpt2.ordered = res_xcell[,angpt2.order]                                      
res_xcell_angpt2.ordered = cbind.data.frame(res_xcell$cell_type, res_xcell_angpt2.ordered)
colnames(res_xcell_angpt2.ordered)[1] <- "cell_type"
#write.table(res_xcell_angpt2.ordered, "Tables/total.res_xcell_angpt2.ordered.txt")


## Step 6.Discriminate primary tumor vs metastasis

# Download GDC data by sample (met, pri)
met.query <- GDCquery(project = "TCGA-SKCM",
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification",
                      workflow.type  = "HTSeq - FPKM",
                      experimental.strategy = "RNA-Seq",
                      sample.type = c("Metastatic"),
                      legacy = FALSE)

pri.query <- GDCquery(project = "TCGA-SKCM",
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification",
                      workflow.type  = "HTSeq - FPKM",
                      experimental.strategy = "RNA-Seq",
                      sample.type = c("Primary Tumor"),
                      legacy = FALSE)

GDCdownload(met.query, method = "api", files.per.chunk = 10)
met.data <- GDCprepare(met.query)
#saveRDS(met.data, "met_TCGA-SKCM_GDCdata.rds")
met.data = readRDS("Data/met_TCGA-SKCM_GDCdata.rds")
GDCdownload(pri.query, method = "api", files.per.chunk = 10)
pri.data <- GDCprepare(pri.query)
#saveRDS(pri.data, "pri_TCGA-SKCM_GDCdata.rds")
pri.data = readRDS("Data/pri_TCGA-SKCM_GDCdata.rds")

met.matrix <- as.data.frame(assays(met.data)$`HTSeq - FPKM`)
met.list = colnames(met.matrix)
pri.matrix <- as.data.frame(assays(pri.data)$`HTSeq - FPKM`)
pri.list = colnames(pri.matrix)


## Step 7. Visualization
library(tidyr)
library(ggplot)
for (i in 1:39){
  print(i)
  str <- res_xcell_angpt2.ordered$cell_type[i]
  d <- res_xcell_angpt2.ordered[i,] %>%
    gather(sample, score, -cell_type)
  d$type = ifelse(d$sample %in% met.list, "met", 
                  ifelse(d$sample %in% pri.list, "pri", NA))
  d = d[!is.na(d$type),] # 470 obs.
  samplelist <- unlist(d$sample)
  data_angpt2 = data_angpt2 %>% filter(data_angpt2$sample %in% samplelist)
  d = merge(d, data_angpt2, by = "sample")
  d = d %>% arrange(ANGPT2)
  
  print(ggplot(d, aes(x= log2(ANGPT2), y=score)) +
          stat_smooth(method = 'lm', aes(linetype = type), color = "gray30", fill = "grey70", size = 3, alpha = 0.5, show.legend = T) + ## linear regression
          geom_point(aes(shape = type), color = "grey25", size=1.1, stroke =0.55, alpha = 0.7) + # Make dots bigger
          scale_shape_manual(values=c(5, 4))+
          geom_vline(xintercept= c(-0.5407935, 2.739935), color = c("darkred", "darkblue"), linetype="dashed", size = 1.75, alpha = 0.65) + #met
          geom_vline(xintercept= c(0.2470981, 2.800226), color = c("darkred", "darkblue"), linetype="solid", size = 1.75, alpha = 0.65) +
          annotate("text", x = -2.75, y = max(d$score), label = "\n<50th", size = 5, color = "grey1" ) +
          annotate("text", x = -2.75, y = 0.99*max(d$score), label = "\n\n\n\n(met.n = 179)\n(pri.n = 51)",
                   size = 4, color = "grey1", alpha = 0.8, fontface = "italic") +
          annotate("text", x = 5.5, y = max(d$score) , label = "\n>90th", size = 5, color = "grey1" ) +
          annotate("text", x = 5.5, y = 0.99*max(d$score), label = "\n\n\n\n(met.n = 36)\n(pri.n = 11)",
                   size = 4, color = "grey1", alpha = 0.8, fontface = "italic") +
          theme_bw() +
          theme(panel.background = element_blank(),
                axis.text.x = element_blank(),
                axis.text.y = element_text(color = "gray1", angle = 90, vjust = 0.5, hjust=1, size = 10),
                axis.title.x = element_text(vjust = -0.1, size = 13),
                panel.grid.major.y = element_line(size = 0.4, linetype = "solid",color = "gray90"),
                panel.grid.major.x = element_blank(),
                legend.title = element_text(color = "gray10", size = 13),
                legend.text = element_text(color = "gray10", size = 13),
                axis.line = element_line(size = 0.5, linetype = "solid", color = "black"),
                axis.ticks = element_blank(),
                plot.title = element_text(color = "Black", size = 16, face = "bold"),
                axis.title = element_text(color = "gray10", size = 14),
                legend.position="bottom") +
          labs(x = expression(log[2]("ANGPT2 Expr.")),
               y = 'xCell score',
               title = str) +
          guides(shape=guide_legend(override.aes=list(size=1), keywidth = 2, keyheight = 1))) 
  Sys.sleep(3)
  
  grid.ls(grid.force())    # To get the names of all the grobs in the ggplot
  # The edit - to set the size of the point in the legend to 4 mm
  grid.gedit("key-1-3-2.2-4-2-4", size = unit(3, "mm")) 
  grid.gedit("key-1-7-2.2-8-2-8", size = unit(3, "mm")) 
  g <- grid.grab()
  #ggsave(plot = g, paste("Figures/xcell_lm/lm_by.sample.", str, ".0924.pdf", sep = ""), width = 12, height = 9)
  #save_plot(paste("Figures/xcell_lm/lm_by.sample.", str, ".0924.pdf", sep = ""), last_plot(), base_height = 4.5, base_width = 6)
}


## Step 8.Calculate lmfit value

library(broom)
options(scipen = 999)

# total
total.res = read.table("Tables/total.res_xcell_angpt2.ordered.txt")
total.res = data.frame(t(total.res))
colnames(total.res) <- total.res[1,]
total.res = total.res[2:nrow(total.res),]
total.data = read.table("Tables/total.data_angpt2.ordered.txt")
total = cbind.data.frame(total.data, total.res)

str <- colnames(total)
str2 = gsub(" ", "_", str)
str2 = gsub("\\+", "", str2)
str2 = gsub("-", "_", str2)
colnames(total) <- str2
colnames(total)[6] = "T_cell_CD4_non_regulatory"
colnames(total)[37] ="T_cell_regs"


total[,2:ncol(total)] <- sapply(total[,2:ncol(total)], as.numeric)

lmfit <- lm(as.formula(paste(colnames(total)[1], "~", 
                             paste(colnames(total)[2:ncol(total)], collapse = " + "), sep = "")), data=total)

## met
met.res = read.table("Tables/met.res_xcell_angpt2.ordered.txt")
met.res = data.frame(t(met.res))
colnames(met.res) <- met.res[1,]
met.res = met.res[2:nrow(met.res),]
met.data = read.table("Tables/met.data_angpt2.ordered.txt")
met = cbind.data.frame(met.data, met.res)

str <- colnames(met)
str2 = gsub(" ", "_", str)
str2 = gsub("\\+", "", str2)
str2 = gsub("-", "_", str2)
colnames(met) <- str2
colnames(met)[6] = "T_cell_CD4_non_regulatory"
colnames(met)[37] ="T_cell_regs"
met[,2:ncol(met)] <- sapply(met[,2:ncol(met)], as.numeric)

met = met[,1:37] # Exclude score column

met.lmfit <- lm(as.formula(paste(colnames(met)[1], "~", 
                                 paste(colnames(met)[2:ncol(met)], collapse = " + "), sep = "")), data=met)
met.lmfit

met = met[,c(1, 2, 5, 7:9, 11:21, 23:37)] # individual samples
met.glmfit <- glm(as.formula(paste(colnames(met)[1], "~", 
                                   paste(colnames(met)[2:ncol(met)], collapse = " + "), sep = "")), data=met)
met.glmfit
tidy(met.glmfit)

## pri
pri.res = read.table("Tables/pri.res_xcell_angpt2.ordered.txt")
pri.res = data.frame(t(pri.res))
colnames(pri.res) <- pri.res[1,]
pri.res = pri.res[2:nrow(pri.res),]
pri.data = read.table("Tables/pri.data_angpt2.ordered.txt")
pri = cbind.data.frame(pri.data, pri.res)

str <- colnames(pri)
str2 = gsub(" ", "_", str)
str2 = gsub("\\+", "", str2)
str2 = gsub("-", "_", str2)
colnames(pri) <- str2
colnames(pri)[6] = "T_cell_CD4_non_regulatory"
colnames(pri)[37] ="T_cell_regs"
pri[,2:ncol(pri)] <- sapply(pri[,2:ncol(pri)], as.numeric)

pri = pri[,1:37] # Exclude score column

pri.lmfit <- lm(as.formula(paste(colnames(pri)[1], "~", 
                                 paste(colnames(pri)[2:ncol(pri)], collapse = " + "), sep = "")), data=pri)
pri.lmfit

pri = pri[,c(1, 2, 5, 7:9, 11:21, 23:37)] # Exclude score column
pri.lmfit <- lm(as.formula(paste(colnames(pri)[1], "~", 
                                 paste(colnames(pri)[2:ncol(pri)], collapse = " + "), sep = "")), data=pri)
pri.lmfit



## Sampletype

# Analyze Immune deconvolution with tool "Immunedeconv"


## Step 1. Download the TCGA-SKCM data with TCGAbiolinks
library(BiocManager)
#BiocManager::install(c("TCGAbiolinks", "DT"))
library(TCGAbiolinks)
library(dplyr)
library(DT)

# Download GDC data by sample (met, pri)
met.query <- GDCquery(project = "TCGA-SKCM",
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification",
                      workflow.type  = "HTSeq - FPKM",
                      experimental.strategy = "RNA-Seq",
                      sample.type = c("Metastatic"),
                      legacy = FALSE)

pri.query <- GDCquery(project = "TCGA-SKCM",
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification",
                      workflow.type  = "HTSeq - FPKM",
                      experimental.strategy = "RNA-Seq",
                      sample.type = c("Primary Tumor"),
                      legacy = FALSE)

GDCdownload(met.query, method = "api", files.per.chunk = 10)
met.data <- GDCprepare(met.query)
#saveRDS(met.data, "met_TCGA-SKCM_GDCdata.rds")
met.data = readRDS("met_TCGA-SKCM_GDCdata.rds")
GDCdownload(pri.query, method = "api", files.per.chunk = 10)
pri.data <- GDCprepare(pri.query)
#saveRDS(pri.data, "pri_TCGA-SKCM_GDCdata.rds")
pri.data = readRDS("pri_TCGA-SKCM_GDCdata.rds")

l <- list()
l[[1]] <- met.data
l[[2]] <- pri.data

for (j in 1:2){
  data <- l[[j]]
  ## Step 2. Get expression matrix from TCGA data
  #remotes::install_github("icbi-lab/immunedeconv", force = TRUE)
  library(SummarizedExperiment)
  matrix_data <- as.data.frame(assays(data)$`HTSeq - FPKM`)
  
  ## Step 3. Convert rownames from gene id to HGNC gene symbol
  option_tx = F
  if (option_tx){
    gm = rtracklayer::import('~/Dropbox/Resources/Gencode_hg38_v32/gencode.v32.primary_assembly.annotation.gtf.gz')
    gm1 = as.data.frame(gm) %>% filter(type=='transcript') %>% mutate(tss=ifelse(strand=='+', start, end)) %>% dplyr::select(transcript_id, gene_id, gene_name, tss, gene_type, transcript_type)
    gm1$transcript_id = do.call(rbind.data.frame, strsplit(gm1$transcript_id, '.', fixed = T))[[1]]
    gm1$gene_id = do.call(rbind.data.frame, strsplit(gm1$gene_id, '.', fixed = T))[[1]]
    saveRDS(gm1, '~/Dropbox/Resources/genes.Rdata')
  } else{
    gm1 = readRDS('~/Dropbox/Resources/genes.Rdata')
  }
  gm1$gene_id = do.call(rbind.data.frame, strsplit(gm1$gene_id, '.', fixed = T))[[1]]
  gene_coding = gm1 %>% filter(gene_type=='protein_coding') %>% dplyr::select(gene_name, gene_id) %>% unique
  gene_coding = gene_coding[!gene_coding$gene_name %in% c("MATR3", "PINX1", "TMSB15B"),] # duplicates in TCGA data : MATR3, PINX1, TMSB15B
  
  d = as.data.frame(rownames(matrix_data))
  colnames(d) = "gene_id" # 56602 obs.
  d = merge(d, gene_coding, by = "gene_id") %>% unique
  
  matrix_data = matrix_data %>% filter(rownames(.) %in% d$gene_id)
  rownames(matrix_data) <- d$gene_name 
  
  ## Step 4.Start deconvolution
  
  ### QuanTIseq
  #res_quantiseq = immunedeconv::deconvolute(matrix_data, "quantiseq")
  ## xCell
  load(file = "Data/xCell.data.rda")
  res_xcell = immunedeconv::deconvolute(matrix_data, "xcell")
  
  ## ANGPT2
  data_angpt2 <- matrix_data %>% filter(rownames(matrix_data) == "ANGPT2") %>% t %>% data.frame %>% arrange(ANGPT2)
  type = ifelse(j == 1, "met", "pri")
  # write.table(data_angpt2, paste("Tables/", type, ".data_angpt2.ordered.txt", sep = ""))
  
  angpt2.order <- rownames(data_angpt2)
  res_xcell_angpt2.ordered = res_xcell[,angpt2.order]
  res_xcell_angpt2.ordered = cbind.data.frame(res_xcell$cell_type, res_xcell_angpt2.ordered)
  colnames(res_xcell_angpt2.ordered)[1] <- "cell_type"
  type = ifelse(j == 1, "met", "pri")
  # write.table(res_xcell_angpt2.ordered, paste("Tables/", type, ".res_xcell_angpt2.ordered.txt", sep = ""))
  
  ## ANGPT2 & ANGPT2 related genes
  df = read_excel("Tables/FindAllMarkers_endothelial.scobject_1110.xlsx")
  df = df %>% filter(cluster == "high") %>% arrange(desc(avg_log2FC))
  df = df[1:15,]
  data_angpt2 <- matrix_data %>% filter(rownames(matrix_data) %in% df$gene) %>% t %>% data.frame %>% arrange(ANGPT2)
  type = ifelse(j == 1, "met", "pri")
  # write.table(data_angpt2, paste("Tables/", type, ".data_angpt2.ordered.txt", sep = ""))
  
  data_angpt2$Median <- apply(data_angpt2[,1:13],1,median)
  
  angpt2.order <- rownames(data_angpt2)
  res_xcell_angpt2.ordered = res_xcell[,angpt2.order]
  res_xcell_angpt2.ordered = cbind.data.frame(res_xcell$cell_type, res_xcell_angpt2.ordered)
  colnames(res_xcell_angpt2.ordered)[1] <- "cell_type"
  type = ifelse(j == 1, "met", "pri")
  # write.table(res_xcell_angpt2.ordered, paste("Tables/", type, ".res_xcell_angpt2.ordered.txt", sep = ""))
  
  ## ANGPT2
  data_angpt2 <- matrix_data %>% filter(rownames(matrix_data) == "ANGPT2") %>% t %>% data.frame %>% arrange(ANGPT2)
  type = ifelse(j == 1, "met", "pri")
  # write.table(data_angpt2, paste("Tables/", type, ".data_angpt2.ordered.txt", sep = ""))
  
  res_xcell.t_ordered = cbind.data.frame(data_angpt2, res_xcell.t) %>% arrange(ANGPT2)
  
  # angpt2.order <- rownames(data_angpt2)
  # res_xcell_angpt2.ordered = res_xcell[,angpt2.order]
  # res_xcell_angpt2.ordered = cbind.data.frame(res_xcell$cell_type, res_xcell_angpt2.ordered)
  # colnames(res_xcell_angpt2.ordered)[1] <- "cell_type"
  # type = ifelse(j == 1, "met", "pri")
  # write.table(res_xcell_angpt2.ordered, paste("Tables/", type, ".res_xcell_angpt2.ordered.txt", sep = ""))
  
  library(ggplot2)
  library(tidyr)
  library(cowplot)
  library(RColorBrewer)
  my.cols <- brewer.pal(12, "Paired")
  my.cols1 <- my.cols[c(2,4,6,8,10,12)]
  my.cols2 <- my.cols[c(1, 3, 5, 7, 9, 11)]
  
  d = res_xcell.t_ordered[,1:2] 
  d %>% ggplot(aes(x= log2(ANGPT2), y=d[,2])) +
    geom_point(size=2, color = my.cols2[i%%6 + 1]) + 
    #     geom_vline(xintercept= ifelse(j == 1, 179, 51), linetype=4, colour = "blue", size = 1, alpha = 0.7 ) +
    #    geom_vline(xintercept= ifelse(j == 1, 179, 51), color = "grey1", linetype="dashed", size = 0.5, alpha = 0.7) +
    #    geom_vline(xintercept= ifelse(j == 1, 331, 92), color = "grey1", linetype="dashed", size = 0.5, alpha = 0.7) +
    #    annotate("text", x = 14, y = 77, label = "<50th (n = 179)", size = 3, color = "gray10" ) +
    # geom_smooth(aes(group = 1), size = 2,  color = my.cols[i%%6 + 1]) + ## loess function
    stat_smooth(method = 'lm', aes(group = 1), size = 2,  color = my.cols1[i%%6 + 1]) + ## linear regression
    #    facet_wrap(~cell_type, scales="free_x", ncol=3) +
    scale_color_brewer(palette="Paired", guide= "none") +
    #     coord_flip() +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10),
          #panel.border = element_blank(),
          panel.grid.major = element_blank(),
          strip.text = element_text(size=14),
          axis.ticks = element_blank(),
          axis.title = element_text(size = 12)) +
    labs(x = 'ANGPT2 Expression level',
         y = 'xCell score')
  
  
  
  for (i in 1:39){
    print(i)
    p <- res_xcell_angpt2.ordered[i,] %>%
      gather(sample, score, -cell_type) %>%
      ggplot(aes(x= factor(sample, level = angpt2.order), y=score), color=cell_type) +
      geom_point(size=2, color = my.cols2[i%%6 + 1]) + 
      #     geom_vline(xintercept= ifelse(j == 1, 179, 51), linetype=4, colour = "blue", size = 1, alpha = 0.7 ) +
      geom_vline(xintercept= ifelse(j == 1, 179, 51), color = "grey1", linetype="dashed", size = 0.5, alpha = 0.7) +
      geom_vline(xintercept= ifelse(j == 1, 331, 92), color = "grey1", linetype="dashed", size = 0.5, alpha = 0.7) +
      annotate("text", x = 14, y = 77, label = "<50th (n = 179)", size = 3, color = "gray10" ) +
      # geom_smooth(aes(group = 1), size = 2,  color = my.cols[i%%6 + 1]) + ## loess function
      stat_smooth(method = 'lm', aes(group = 1), size = 2,  color = my.cols1[i%%6 + 1]) + ## linear regression
      facet_wrap(~cell_type, scales="free_x", ncol=3) +
      scale_color_brewer(palette="Paired", guide= "none") +
      #     coord_flip() +
      theme_bw() +
      theme(axis.text.x = element_blank(),
            axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10),
            #panel.border = element_blank(),
            panel.grid.major = element_blank(),
            strip.text = element_text(size=14),
            axis.ticks = element_blank(),
            axis.title = element_text(size = 12)) +
      labs(x = 'ANGPT2 Expression level',
           y = 'xCell score')
    
    str <- res_xcell_angpt2.ordered$cell_type[i]
    
    save_plot(paste(type, str, "0830.png", sep = "."), p, base_height = 5, base_width = 5)
  }
}










## Step 1. Download the TCGA-SKCM data with TCGAbiolinks
library(BiocManager)
#BiocManager::install(c("TCGAbiolinks", "DT"))
library(TCGAbiolinks)
library(dplyr)
library(DT)

# query <- GDCquery(project = "TCGA-SKCM",
#                   data.category = "Transcriptome Profiling",
#                   data.type = "Gene Expression Quantification",
#                   workflow.type  = "HTSeq - FPKM",
#                   experimental.strategy = "RNA-Seq",
#                   legacy = FALSE)
# GDCdownload(query, method = "api", files.per.chunk = 10)
# data <- GDCprepare(query)
#saveRDS(data, "total_TCGA-SKCM_GDCdata.rds")
data = readRDS("total_TCGA-SKCM_GDCdata.rds")

## Step 2. Get expression matrix from TCGA data
#remotes::install_github("icbi-lab/immunedeconv", force = TRUE)
library(SummarizedExperiment)
matrix_data <- as.data.frame(assays(data)$`HTSeq - FPKM`)

## Step 3. Convert rownames from gene id to HGNC gene symbol
option_tx = F
if (option_tx){
  gm = rtracklayer::import('~/Dropbox/Resources/Gencode_hg38_v32/gencode.v32.primary_assembly.annotation.gtf.gz')
  gm1 = as.data.frame(gm) %>% filter(type=='transcript') %>% mutate(tss=ifelse(strand=='+', start, end)) %>% dplyr::select(transcript_id, gene_id, gene_name, tss, gene_type, transcript_type)
  gm1$transcript_id = do.call(rbind.data.frame, strsplit(gm1$transcript_id, '.', fixed = T))[[1]]
  gm1$gene_id = do.call(rbind.data.frame, strsplit(gm1$gene_id, '.', fixed = T))[[1]]
  saveRDS(gm1, '~/Dropbox/Resources/genes.Rdata')
} else{
  gm1 = readRDS('~/Dropbox/Resources/genes.Rdata')
}
gm1$gene_id = do.call(rbind.data.frame, strsplit(gm1$gene_id, '.', fixed = T))[[1]]
gene_coding = gm1 %>% filter(gene_type=='protein_coding') %>% select(gene_name, gene_id) %>% unique
gene_coding = gene_coding[!gene_coding$gene_name %in% c("MATR3", "PINX1", "TMSB15B"),] # duplicates in TCGA data : MATR3, PINX1, TMSB15B

d = as.data.frame(rownames(matrix_data))
colnames(d) = "gene_id" # 56602 obs.
d = merge(d, gene_coding, by = "gene_id") %>% unique

matrix_data = matrix_data %>% filter(rownames(.) %in% d$gene_id)
rownames(matrix_data) <- d$gene_name 


## Step 4.Start deconvolution
### QuanTIseq
#res_quantiseq = immunedeconv::deconvolute(matrix_data, "quantiseq")

### xCell
load(file = "xCell.data.rda")
res_xcell = immunedeconv::deconvolute(matrix_data, "xcell")
res_xcell.t = data.frame(t(res_xcell))
colnames(res_xcell.t) = res_xcell.t[1,]
res_xcell.t = res_xcell.t[2:nrow(res_xcell.t), ]


## Step 5. Need to discriminate primary tumor vs metastasis

met.data = readRDS("met_TCGA-SKCM_GDCdata.rds")
met.matrix <- as.data.frame(assays(met.data)$`HTSeq - FPKM`)
met.list = colnames(met.matrix)
pri.data = readRDS("pri_TCGA-SKCM_GDCdata.rds")
pri.matrix <- as.data.frame(assays(pri.data)$`HTSeq - FPKM`)
pri.list = colnames(pri.matrix)

res_xcell.t$sampletype = ifelse(rownames(res_xcell.t) %in% met.list, "met",
                                ifelse(rownames(res_xcell.t) %in% pri.list, "pri", NA))
res_xcell.t = res_xcell.t[!is.na(res_xcell.t$sampletype), ]
## ANGPT2
data_angpt2 <- matrix_data %>% filter(rownames(matrix_data) == "ANGPT2") %>% t %>% data.frame %>% arrange(ANGPT2)
data_angpt2 = data_angpt2 %>% filter(rownames(data_angpt2) %in% rownames(res_xcell.t))
angpt2.order <- rownames(data_angpt2)
res_xcell.ordered = res_xcell.t[angpt2.order,]  # 470 obs.
res_xcell.ordered = cbind.data.frame(data_angpt2, res_xcell.ordered)
write.table(res_xcell.ordered, "Tables/total.res_xcell.sampletype.ordered.txt")


## Step 6. Visualization

library(ggplot2)
library(tidyr)
library(cowplot)
library(RColorBrewer)
my.cols <- brewer.pal(12, "Paired")
my.cols1 <- my.cols[c(2,4,6,8,10,12)]
my.cols2 <- my.cols[c(1, 3, 5, 7, 9, 11)]

i = 1
res_xcell.ordered %>% ggplot(aes(x = ANGPT2, y = `Myeloid dendritic cell activated`)) +
  geom_point(size = 2, color = my.cols2[i%%6 + 1])
, shape = res_xcell.ordered$sampletype)
res_xcell.ordered = res_xcell.ordered %>% filter(ANGPT2 <= 100)

res_xcell.ordered[i,] %>%
  gather(sample, score, -cell_type) %>%
  ggplot(aes(x= factor(sample, level = angpt2.order), y=score), color=cell_type) +
  geom_point(size=2, color = my.cols2[i%%6 + 1]) +
  geom_hline(yintercept = quantile(ANG, knots), linetype = 'dashed') +
  # geom_smooth(aes(group = 1), size = 2,  color = my.cols[i%%6 + 1]) + ## loess function
  stat_smooth(method = 'lm', aes(group = 1), size = 2,  color = my.cols1[i%%6 + 1]) + ## linear regression
  facet_wrap(~cell_type, scales="free_x", ncol=3) +
  scale_color_brewer(palette="Paired", guide= "none") +
  coord_flip() +
  theme_bw() +
  stat_regline_equation(aes(label = ..eq.label..)) +
  stat_regline_equation(aes(label = ..rr.label..)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10),
        axis.text.y = element_blank(),
        #panel.border = element_blank(),
        panel.grid.major = element_blank(),
        strip.text = element_text(size=14),
        axis.ticks = element_blank(),
        axis.title = element_text(size = 12)) +
  labs(x = 'ANGPT2 Expression level',
       y = 'xCell score')

str <- res_xcell_angpt2.ordered$cell_type[i]
save_plot(paste("Figures/xcell_lm/lm_", str, ".png", sep = ""), p, base_height = 5, base_width = 5)




library(ggpubr)
for (i in 1:39){
  print(i)
  p <- res_xcell_angpt2.ordered[i,] %>%
    gather(sample, score, -cell_type) %>%
    ggplot(aes(x= factor(sample, level = angpt2.order), y=score), color=cell_type) +
    geom_point(size=2, color = my.cols2[i%%6 + 1]) +
    geom_hline(yintercept = quantile(ANG, knots), linetype = 'dashed') +
    # geom_smooth(aes(group = 1), size = 2,  color = my.cols[i%%6 + 1]) + ## loess function
    stat_smooth(method = 'lm', aes(group = 1), size = 2,  color = my.cols1[i%%6 + 1]) + ## linear regression
    facet_wrap(~cell_type, scales="free_x", ncol=3) +
    scale_color_brewer(palette="Paired", guide= "none") +
    coord_flip() +
    theme_bw() +
    stat_regline_equation(aes(label = ..eq.label..)) +
    stat_regline_equation(aes(label = ..rr.label..)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10),
          axis.text.y = element_blank(),
          #panel.border = element_blank(),
          panel.grid.major = element_blank(),
          strip.text = element_text(size=14),
          axis.ticks = element_blank(),
          axis.title = element_text(size = 12)) +
    labs(x = 'ANGPT2 Expression level',
         y = 'xCell score')
  
  str <- res_xcell_angpt2.ordered$cell_type[i]
  save_plot(paste("Figures/xcell_lm/lm_", str, ".png", sep = ""), p, base_height = 5, base_width = 5)
}

tmp = res_xcell_angpt2.ordered[i,] %>%  gather(sample, score, -cell_type) 
#################################


tmp = data_angpt2 %>% gather(sample, score, -name)

ggplot(tmp, aes(x= name, y=score, color=sample)) + geom_point(size=2)




## ANGPT2 related


# Analyze Immune deconvolution with tool "Immunedeconv"

## Step 1. Download the TCGA-SKCM data with TCGAbiolinks
library(BiocManager)
#BiocManager::install(c("TCGAbiolinks", "DT"))
library(TCGAbiolinks)
library(dplyr)
library(DT)

# Download GDC data by sample (met, pri)
met.query <- GDCquery(project = "TCGA-SKCM",
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification",
                      workflow.type  = "HTSeq - FPKM",
                      experimental.strategy = "RNA-Seq",
                      sample.type = c("Metastatic"),
                      legacy = FALSE)

pri.query <- GDCquery(project = "TCGA-SKCM",
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification",
                      workflow.type  = "HTSeq - FPKM",
                      experimental.strategy = "RNA-Seq",
                      sample.type = c("Primary Tumor"),
                      legacy = FALSE)

GDCdownload(met.query, method = "api", files.per.chunk = 10)
met.data <- GDCprepare(met.query)
#saveRDS(met.data, "met_TCGA-SKCM_GDCdata.rds")
met.data = readRDS("Data/met_TCGA-SKCM_GDCdata.rds")
GDCdownload(pri.query, method = "api", files.per.chunk = 10)
pri.data <- GDCprepare(pri.query)
#saveRDS(pri.data, "pri_TCGA-SKCM_GDCdata.rds")
pri.data = readRDS("Data/pri_TCGA-SKCM_GDCdata.rds")
data = readRDS("Data/total_TCGA-SKCM_GDCdata.rds")

l <- list()
l[[1]] <- met.data
l[[2]] <- pri.data

for (j in 1:2){
  data <- l[[j]]
  ## Step 2. Get expression matrix from TCGA data
  #remotes::install_github("icbi-lab/immunedeconv", force = TRUE)
  library(SummarizedExperiment)
  matrix_data <- as.data.frame(assays(data)$`HTSeq - FPKM`)
  bulk2 = read.csv("~/Downloads/TCGA.SKCM.sampleMap_HiSeqV2", sep = "\t")
  rownames(bulk2) = bulk2$sample
  bulk2$sample <- NULL
  
  library(stringi)
  met = colnames(bulk.met) %>% gsub("-", ".", .) %>% strtrim(., 15)
  met = intersect(met, colnames(bulk2))
  bulk.met = bulk2[met]
  
  pri = colnames(bulk.pri) %>% gsub("-", ".", .) %>% strtrim(., 15)
  pri = intersect(pri, colnames(bulk2))
  bulk.pri = bulk2[pri]
  
  ## Step 3. Convert rownames from gene id to HGNC gene symbol
  option_tx = F
  if (option_tx){
    gm = rtracklayer::import('~/Dropbox/Resources/Gencode_hg38_v32/gencode.v32.primary_assembly.annotation.gtf.gz')
    gm1 = as.data.frame(gm) %>% filter(type=='transcript') %>% mutate(tss=ifelse(strand=='+', start, end)) %>% dplyr::select(transcript_id, gene_id, gene_name, tss, gene_type, transcript_type)
    gm1$transcript_id = do.call(rbind.data.frame, strsplit(gm1$transcript_id, '.', fixed = T))[[1]]
    gm1$gene_id = do.call(rbind.data.frame, strsplit(gm1$gene_id, '.', fixed = T))[[1]]
    saveRDS(gm1, '~/Dropbox/Resources/genes.Rdata')
  } else{
    gm1 = readRDS('~/Dropbox/Resources/genes.Rdata')
  }
  gm1$gene_id = do.call(rbind.data.frame, strsplit(gm1$gene_id, '.', fixed = T))[[1]]
  gene_coding = gm1 %>% filter(gene_type=='protein_coding') %>% dplyr::select(gene_name, gene_id) %>% unique
  gene_coding = gene_coding[!gene_coding$gene_name %in% c("MATR3", "PINX1", "TMSB15B"),] # duplicates in TCGA data : MATR3, PINX1, TMSB15B
  
  d = as.data.frame(rownames(matrix_data))
  colnames(d) = "gene_id" # 56602 obs.
  d = merge(d, gene_coding, by = "gene_id") %>% unique
  
  matrix_data = matrix_data %>% filter(rownames(.) %in% d$gene_id)
  rownames(matrix_data) <- d$gene_name 
  
  ## Step 4.Start deconvolution
  
  ### QuanTIseq
  #res_quantiseq = immunedeconv::deconvolute(matrix_data, "quantiseq")
  ## xCell
  load(file = "Data/xCell.data.rda")
  res_xcell = immunedeconv::deconvolute(matrix_data, "xcell")
  
  ## ANGPT2 & ANGPT2 related genes
  df = read_excel("Tables/FindAllMarkers_endothelial.scobject_1110.xlsx")
  df = read_excel("Tables/FindAllMarkers_endothelial.threshold.scobject.xlsx")
  #df = df %>% filter(cluster == "high") %>% arrange(desc(avg_log2FC))
  #df = df[1:15,]
  data_angpt2 <- matrix_data %>% filter(rownames(matrix_data) %in% df$gene) %>% t %>% data.frame %>% arrange(ANGPT2)
  data_angpt2 = data_angpt2[1:471,]
  data_angpt2 <- bulk2 %>% filter(rownames( bulk2) %in% df$gene) %>% t %>% data.frame %>% arrange(ANGPT2)
  data_angpt2 <- bulk.met %>% filter(rownames( bulk.met) %in% df$gene) %>% t %>% data.frame %>% arrange(ANGPT2)
  data_angpt2 <- bulk.pri %>% filter(rownames( bulk.pri) %in% df$gene) %>% t %>% data.frame %>% arrange(ANGPT2)
  
  result = list()
  for (i in colnames(data_angpt2)){
    p <- ggscatter(data_angpt2, x = "ANGPT2", y = i, 
                   add = "reg.line", conf.int = TRUE, 
                   cor.coef = TRUE, cor.method = "pearson",
                   xlab = "ANGPT2", ylab = i)
    plot(p)
    res <- cor.test(data_angpt2$ANGPT2, data_angpt2[,i], 
                    method = "pearson")
    if(res$p.value < 0.001){
      if(abs(res$estimate) >= 0.25){
        ifelse(abs(res$estimate) >= 0.3,
               ifelse(abs(res$estimate) >= 0.4, result[["0.4"]]<-append(i, result[["0.4"]]), 
                      result[["0.3"]]<-append(i, result[["0.3"]])), result[["0.25"]]<-append(i, result[["0.25"]]))
      } 
      if(abs(res$estimate) >= 0.4){
        save_plot(paste('Figures/Pearson/strong.pri.ANGPT2_pearson_endo.marker_', i, '.1111.pdf', sep = ""), last_plot(), base_height = 5, base_width = 5.5)
      }
    } 
  }
  
  for (i in colnames(data_angpt2)){
    p <- ggscatter(data_angpt2, x = "ANGPT2", y = i, 
                   add = "reg.line", conf.int = TRUE, 
                   cor.coef = TRUE, cor.method = "pearson",
                   xlab = "ANGPT2", ylab = i)
    plot(p)
    res <- cor.test(data_angpt2$ANGPT2, data_angpt2[,i], 
                    method = "pearson")
    if(res$p.value < 0.001){
      save_plot(paste('Figures/Pearson/ANGPT2_pearson_0.001_endo.marker_', i, '.pdf', sep = ""), last_plot(), base_height = 5, base_width = 5.5)
    }
  }
  ggplot(data_angpt2, aes(y = ANGPT2)) + geom_boxplot()
  quantile(data_angpt2$ANGPT2, c(0, 0.25, 0.5, 0.75, 1))
  data_angpt2.quant = data_angpt2[data_angpt2$ANGPT2 <= 2.68759620,]
  
  for (i in colnames(data_angpt2.quant)){
    p <- ggscatter(data_angpt2.quant, x = "ANGPT2", y = i, 
                   add = "reg.line", conf.int = TRUE, 
                   cor.coef = TRUE, cor.method = "pearson",
                   xlab = "ANGPT2", ylab = i)
    plot(p)
    res <- cor.test(data_angpt2.quant$ANGPT2, data_angpt2.quant[,i], 
                    method = "pearson")
    if(res$p.value < 0.001 & abs(res$estimate) >= 0.25){
      save_plot(paste('Figures/Pearson/ANGPT2_pearson_quant.str_endo.marker_', i, '.pdf', sep = ""), last_plot(), base_height = 5, base_width = 5.5)
    }
  }
  
  
  
  data_angpt2.log = log(data_angpt2 + 1)
  for (i in colnames(data_angpt2.log)){
    p <- ggscatter(data_angpt2.log, x = "ANGPT2", y = i, 
                   add = "reg.line", conf.int = TRUE, 
                   cor.coef = TRUE, cor.method = "pearson",
                   xlab = "ANGPT2", ylab = i)
    plot(p)
    res <- cor.test(data_angpt2.log$ANGPT2, data_angpt2.log[,i], 
                    method = "pearson")
    if(res$p.value < 0.01){
      save_plot(paste('Figures/Pearson/ANGPT2_pearson_endo.marker_log.', i, '.pdf', sep = ""), last_plot(), base_height = 5, base_width = 5.5)
    }
  }
  
  
  
  type = ifelse(j == 1, "met", "pri")
  # write.table(data_angpt2, paste("Tables/", type, ".data_angpt2.ordered.txt", sep = ""))
  
  data_angpt2$Median <- apply(data_angpt2[,1:13],1,median)
  
  angpt2.order <- rownames(data_angpt2)
  res_xcell_angpt2.ordered = res_xcell[,angpt2.order]
  res_xcell_angpt2.ordered = cbind.data.frame(res_xcell$cell_type, res_xcell_angpt2.ordered)
  colnames(res_xcell_angpt2.ordered)[1] <- "cell_type"
  type = ifelse(j == 1, "met", "pri")
  # write.table(res_xcell_angpt2.ordered, paste("Tables/", type, ".res_xcell_angpt2.ordered.txt", sep = ""))
  
  ## Step 7. Visualization
  library(tidyr)
  library(ggplot)
  for (i in 1:39){
    print(i)
    str <- res_xcell_angpt2.ordered$cell_type[i]
    d <- res_xcell_angpt2.ordered[i,] %>%
      gather(sample, score, -cell_type)
    d$type = ifelse(d$sample %in% met.list, "met", 
                    ifelse(d$sample %in% pri.list, "pri", NA))
    d = d[!is.na(d$type),] # 470 obs.
    samplelist <- unlist(d$sample)
    data_angpt2 = data_angpt2 %>% filter(data_angpt2$sample %in% samplelist)
    d = merge(d, data_angpt2, by = "sample")
    d = d %>% arrange(ANGPT2)
    
    print(ggplot(d, aes(x= log2(ANGPT2), y=score)) +
            stat_smooth(method = 'lm', aes(linetype = type), color = "gray30", fill = "grey70", size = 3, alpha = 0.5, show.legend = T) + ## linear regression
            geom_point(aes(shape = type), color = "grey25", size=1.1, stroke =0.55, alpha = 0.7) + # Make dots bigger
            scale_shape_manual(values=c(5, 4))+
            geom_vline(xintercept= c(-0.5407935, 2.739935), color = c("darkred", "darkblue"), linetype="dashed", size = 1.75, alpha = 0.65) + #met
            geom_vline(xintercept= c(0.2470981, 2.800226), color = c("darkred", "darkblue"), linetype="solid", size = 1.75, alpha = 0.65) +
            annotate("text", x = -2.75, y = max(d$score), label = "\n<50th", size = 5, color = "grey1" ) +
            annotate("text", x = -2.75, y = 0.99*max(d$score), label = "\n\n\n\n(met.n = 179)\n(pri.n = 51)",
                     size = 4, color = "grey1", alpha = 0.8, fontface = "italic") +
            annotate("text", x = 5.5, y = max(d$score) , label = "\n>90th", size = 5, color = "grey1" ) +
            annotate("text", x = 5.5, y = 0.99*max(d$score), label = "\n\n\n\n(met.n = 36)\n(pri.n = 11)",
                     size = 4, color = "grey1", alpha = 0.8, fontface = "italic") +
            theme_bw() +
            theme(panel.background = element_blank(),
                  axis.text.x = element_blank(),
                  axis.text.y = element_text(color = "gray1", angle = 90, vjust = 0.5, hjust=1, size = 10),
                  axis.title.x = element_text(vjust = -0.1, size = 13),
                  panel.grid.major.y = element_line(size = 0.4, linetype = "solid",color = "gray90"),
                  panel.grid.major.x = element_blank(),
                  legend.title = element_text(color = "gray10", size = 13),
                  legend.text = element_text(color = "gray10", size = 13),
                  axis.line = element_line(size = 0.5, linetype = "solid", color = "black"),
                  axis.ticks = element_blank(),
                  plot.title = element_text(color = "Black", size = 16, face = "bold"),
                  axis.title = element_text(color = "gray10", size = 14),
                  legend.position="bottom") +
            labs(x = expression(log[2]("ANGPT2 Expr.")),
                 y = 'xCell score',
                 title = str) +
            guides(shape=guide_legend(override.aes=list(size=1), keywidth = 2, keyheight = 1))) 
    Sys.sleep(3)
    
    grid.ls(grid.force())    # To get the names of all the grobs in the ggplot
    # The edit - to set the size of the point in the legend to 4 mm
    grid.gedit("key-1-3-2.2-4-2-4", size = unit(3, "mm")) 
    grid.gedit("key-1-7-2.2-8-2-8", size = unit(3, "mm")) 
    g <- grid.grab()
    #ggsave(plot = g, paste("Figures/xcell_lm/lm_by.sample.", str, ".0924.pdf", sep = ""), width = 12, height = 9)
    #save_plot(paste("Figures/xcell_lm/lm_by.sample.", str, ".0924.pdf", sep = ""), last_plot(), base_height = 4.5, base_width = 6)
  }
}



## ggplot

## Step 1. Download the TCGA-SKCM data with TCGAbiolinks
library(BiocManager)
#BiocManager::install(c("TCGAbiolinks", "DT"))
library(TCGAbiolinks)
library(dplyr)
library(DT)

# query <- GDCquery(project = "TCGA-SKCM",
#                   data.category = "Transcriptome Profiling",
#                   data.type = "Gene Expression Quantification",
#                   workflow.type  = "HTSeq - FPKM",
#                   experimental.strategy = "RNA-Seq",
#                   legacy = FALSE)
# GDCdownload(query, method = "api", files.per.chunk = 10)
# data <- GDCprepare(query)
#saveRDS(data, "total_TCGA-SKCM_GDCdata.rds")
#data = readRDS("total_TCGA-SKCM_GDCdata.rds")
# 
# met.query <- GDCquery(project = "TCGA-SKCM",
#                       data.category = "Transcriptome Profiling",
#                       data.type = "Gene Expression Quantification",
#                       workflow.type  = "HTSeq - FPKM",
#                       experimental.strategy = "RNA-Seq",
#                       sample.type = c("Metastatic"),
#                       legacy = FALSE)
# pri.query <- GDCquery(project = "TCGA-SKCM",
#                       data.category = "Transcriptome Profiling",
#                       data.type = "Gene Expression Quantification",
#                       workflow.type  = "HTSeq - FPKM",
#                       experimental.strategy = "RNA-Seq",
#                       sample.type = c("Primary Tumor"),
#                       legacy = FALSE)
# 
# GDCdownload(met.query, method = "api", files.per.chunk = 10)
# met.data <- GDCprepare(met.query)
# #saveRDS(met.data, "met_TCGA-SKCM_GDCdata.rds")
# 
# GDCdownload(pri.query, method = "api", files.per.chunk = 10)
# pri.data <- GDCprepare(pri.query)
# #saveRDS(pri.data, "pri_TCGA-SKCM_GDCdata.rds")

met.data = readRDS("met_TCGA-SKCM_GDCdata.rds")
pri.data = readRDS("pri_TCGA-SKCM_GDCdata.rds")


## Step 2. Stratify Ang2 high and low (e.g. >90th (high) vs. <50th (low) expression)2
## Get expression matrix from TCGA data
library(SummarizedExperiment)
met.mtx <- as.data.frame(assays(met.data)$`HTSeq - FPKM`)
pri.mtx <- as.data.frame(assays(pri.data)$`HTSeq - FPKM`)

option_tx = F
if (option_tx){
  gm = rtracklayer::import('~/Dropbox/Resources/Gencode_hg38_v32/gencode.v32.primary_assembly.annotation.gtf.gz')
  gm1 = as.data.frame(gm) %>% filter(type=='transcript') %>% mutate(tss=ifelse(strand=='+', start, end)) %>% dplyr::select(transcript_id, gene_id, gene_name, tss, gene_type, transcript_type)
  gm1$transcript_id = do.call(rbind.data.frame, strsplit(gm1$transcript_id, '.', fixed = T))[[1]]
  gm1$gene_id = do.call(rbind.data.frame, strsplit(gm1$gene_id, '.', fixed = T))[[1]]
  saveRDS(gm1, '~/Dropbox/Resources/genes.Rdata')
} else{
  gm1 = readRDS('~/Dropbox/Resources/genes.Rdata')
}
gm1$gene_id = do.call(rbind.data.frame, strsplit(gm1$gene_id, '.', fixed = T))[[1]]
gene_coding = gm1 %>% filter(gene_type=='protein_coding') %>% select(gene_name, gene_id) %>% unique
gene_coding = gene_coding[!gene_coding$gene_name %in% c("MATR3", "PINX1", "TMSB15B"),] # duplicates in TCGA data : MATR3, PINX1, TMSB15B

met.d = as.data.frame(rownames(met.mtx))
colnames(met.d) = "gene_id" # 56602 obs.
met.d = merge(met.d, gene_coding, by = "gene_id") %>% unique
met.mtx = met.mtx %>% filter(rownames(.) %in% met.d$gene_id)
rownames(met.mtx) <- met.d$gene_name 

pri.d = as.data.frame(rownames(pri.mtx))
colnames(pri.d) = "gene_id" # 56602 obs.
pri.d = merge(pri.d, gene_coding, by = "gene_id") %>% unique
pri.mtx = pri.mtx %>% filter(rownames(.) %in% pri.d$gene_id)
rownames(pri.mtx) <- pri.d$gene_name 

met.mtx.t = as.data.frame(t(met.mtx)) %>% arrange(ANGPT2)
met.mtx.t$qnt <- ntile(met.mtx.t$ANGPT2, 10) 
met.mtx.t = met.mtx.t[c(1:179, 331:366), ] %>% select(ANGPT2, qnt)
table(met.mtx.t$qnt)
met.mtx.t$level = ifelse(met.mtx.t$qnt %in% c(9,10), "high", "low")
table(met.mtx.t$level)


pri.mtx.t = as.data.frame(t(pri.mtx)) %>% arrange(ANGPT2)
pri.mtx.t$qnt <- ntile(pri.mtx.t$ANGPT2, 10) 
pri.mtx.t = pri.mtx.t[c(1:51, 93:103), ] %>% select(ANGPT2, qnt)
table(pri.mtx.t$qnt)
pri.mtx.t$level = ifelse(pri.mtx.t$qnt %in% c(9,10), "high", "low")
table(pri.mtx.t$level)

# Stratify complete.


## Step 4.Start deconvolution

### QuanTIseq
#res_quantiseq = immunedeconv::deconvolute(matrix_data, "quantiseq")

### xCell
load(file = "xCell.data.rda")
res_met.xcell = immunedeconv::deconvolute(met.mtx, "xcell")
res_pri.xcell = immunedeconv::deconvolute(pri.mtx, "xcell")

res_met.xcell.t = data.frame(t(res_met.xcell))
colnames(res_met.xcell.t) = res_met.xcell.t[1,]
res_met.xcell.t = res_met.xcell.t[2:nrow(res_met.xcell.t),]
res_met.xcell.t = res_met.xcell.t[rownames(met.mtx.t),]
res_met.xcell = cbind.data.frame(met.mtx.t, res_met.xcell.t)

res_pri.xcell.t = data.frame(t(res_pri.xcell))
colnames(res_pri.xcell.t) = res_pri.xcell.t[1,]
res_pri.xcell.t = res_pri.xcell.t[2:nrow(res_pri.xcell.t),]
res_pri.xcell.t = res_pri.xcell.t[rownames(pri.mtx.t),]
res_pri.xcell = cbind.data.frame(pri.mtx.t, res_pri.xcell.t)

rm(list=setdiff(ls(), c("res_met.xcell", "res_pri.xcell")))

library(ggplot2)
library(tidyr)
library(cowplot)
library(RColorBrewer)
my.cols <- brewer.pal(12, "Paired")
my.cols1 <- my.cols[c(2,4,6,8,10,12)]
my.cols2 <- my.cols[c(1, 3, 5, 7, 9, 11)]

i = 1
res_met.xcell %>%
  ggplot(aes(x = reorder(`Myeloid dendritic cell activated`, ANGPT2), y = ANGPT2, shape = level)) +
  geom_point(size = 2, color = my.cols2[i%%6 + 1])


gather(sample, score, -cell_type) %>%
  ggplot(aes(x= factor(sample, level = angpt2.order), y=score), color=cell_type) +
  geom_point(size=2, color = my.cols2[i%%6 + 1]) +
  # geom_smooth(aes(group = 1), size = 2,  color = my.cols[i%%6 + 1]) + ## loess function
  stat_smooth(method = 'lm', aes(group = 1), size = 2,  color = my.cols1[i%%6 + 1]) + ## linear regression
  facet_wrap(~cell_type, scales="free_x", ncol=3) +
  scale_color_brewer(palette="Paired", guide= "none") +
  coord_flip() +
  theme_bw() +
  stat_regline_equation(aes(label = ..eq.label..)) +
  stat_regline_equation(aes(label = ..rr.label..)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10),
        axis.text.y = element_blank(),
        #panel.border = element_blank(),
        panel.grid.major = element_blank(),
        strip.text = element_text(size=14),
        axis.ticks = element_blank(),
        axis.title = element_text(size = 12)) +
  labs(x = 'ANGPT2 Expression level',
       y = 'xCell score')




































