
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

## ANGPT2
data_angpt2 <- matrix_data %>% filter(rownames(matrix_data) == "ANGPT2") %>% t %>% data.frame %>% arrange(ANGPT2)
angpt2.order <- rownames(data_angpt2)
res_xcell_angpt2.ordered = res_xcell[,angpt2.order]                                      
res_xcell_angpt2.ordered = cbind.data.frame(res_xcell$cell_type, res_xcell_angpt2.ordered)
colnames(res_xcell_angpt2.ordered)[1] <- "cell_type"

library(ggplot2)
library(tidyr)
library(cowplot)
library(RColorBrewer)
my.cols <- brewer.pal(12, "Paired")
my.cols1 <- my.cols[c(2,4,6,8,10,12)]
my.cols2 <- my.cols[c(1, 3, 5, 7, 9, 11)]

lm_eqn <- function(y,){
  m <- lm(y ~ x, df);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}
library(ggpubr)
for (i in 1:39){
  print(i)
  p <- res_xcell_angpt2.ordered[i,] %>%
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

  str <- res_xcell_angpt2.ordered$cell_type[i]
  save_plot(paste("Figures/xcell_lm/lm_", str, ".png", sep = ""), p, base_height = 5, base_width = 5)
}

tmp = res_xcell_angpt2.ordered[i,] %>%  gather(sample, score, -cell_type) 
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
saveRDS(met.data, "met_TCGA-SKCM_GDCdata.rds")
GDCdownload(pri.query, method = "api", files.per.chunk = 10)
pri.data <- GDCprepare(pri.query)
saveRDS(pri.data, "pri_TCGA-SKCM_GDCdata.rds")

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
  ### xCell
  load(file = "xCell.data.rda")
  res_xcell = immunedeconv::deconvolute(matrix_data, "xcell")
  
  ## ANGPT2
  data_angpt2 <- matrix_data %>% filter(rownames(matrix_data) == "ANGPT2") %>% t %>% data.frame %>% arrange(ANGPT2)
  angpt2.order <- rownames(data_angpt2)
  res_xcell_angpt2.ordered = res_xcell[,angpt2.order]                                      
  res_xcell_angpt2.ordered = cbind.data.frame(res_xcell$cell_type, res_xcell_angpt2.ordered)
  colnames(res_xcell_angpt2.ordered)[1] <- "cell_type"
  
  library(ggplot2)
  library(tidyr)
  library(cowplot)
  
  library(RColorBrewer)
  my.cols <- brewer.pal(12, "Paired")
  my.cols <- my.cols[c(1, 3, 5, 7, 9, 11)]
  
  for (i in 1:39){
    print(i)
    p <- res_xcell_angpt2.ordered[i,] %>%
      gather(sample, score, -cell_type) %>%
      ggplot(aes(x= factor(sample, level = angpt2.order), y=score), color=cell_type) +
      geom_point(size=2, color = my.cols2[i%%6 + 1]) +
      # geom_smooth(aes(group = 1), size = 2,  color = my.cols[i%%6 + 1]) + ## loess function
      stat_smooth(method = 'lm', aes(group = 1), size = 2,  color = my.cols1[i%%6 + 1]) + ## linear regression
      facet_wrap(~cell_type, scales="free_x", ncol=3) +
      scale_color_brewer(palette="Paired", guide= "none") +
      coord_flip() +
      theme_bw() +
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
    type = ifelse(j == 1, "met", "pri")
    save_plot(paste(type, str, "png", sep = "."), p, base_height = 5, base_width = 5)
  }
}
