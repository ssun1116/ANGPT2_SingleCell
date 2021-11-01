

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
met.data = readRDS("met_TCGA-SKCM_GDCdata.rds")
GDCdownload(pri.query, method = "api", files.per.chunk = 10)
pri.data <- GDCprepare(pri.query)
#saveRDS(pri.data, "pri_TCGA-SKCM_GDCdata.rds")
pri.data = readRDS("pri_TCGA-SKCM_GDCdata.rds")

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



