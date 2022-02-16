# T cell exclusion signature -> 1) Metastasis sample 2) Primary sample
# 1. Develop signature for each sample type with Malignant scRNA data
# 2. Compare the signature based on ANGPT2 expression on TCGA-SKCM samples
library(tidyverse)
library(Seurat)
setwd("ImmuneResistance/Code/")
source("ImmRes_source.R")
setwd("../../")
set.seed(1234)

bulk.met <- readRDS("Data/met_TCGA-SKCM_GDCdata.rds")
colnames(bulk.met) = gsub("-", ".", colnames(bulk.met)) 
colnames(bulk.met) = substr(colnames(bulk.met), 1, 15)
bulk.pri <- readRDS("Data/pri_TCGA-SKCM_GDCdata.rds")
colnames(bulk.pri) = gsub("-", ".", colnames(bulk.pri)) 
colnames(bulk.pri) = substr(colnames(bulk.pri), 1, 15)

### DATA -> should use xena data instead of GDC (RSEM-normalized RNA-seq = tpm)
bulk.xena = read.table("Data/TCGA.SKCM.sampleMap_HiSeqV2", header = T) # 475 obs.
rownames(bulk.xena) <- bulk.xena$sample
bulk.xena = bulk.xena[,2:ncol(bulk.xena)]

#bulk.xena.met = bulk.xena[, colnames(bulk.xena) %in% colnames(bulk.met)] # 366 obs.
#bulk.xena.pri = bulk.xena[, colnames(bulk.xena) %in% colnames(bulk.pri)] # 103 obs.


cell.cell.intr<-function(rB,r.sc,cellA.markers,cellB.markers,cellA.name = "malignant",
                         cellB.name = "T.cell",bulk.confounders = NULL,sc.confounders = NULL,
                         fileName = NULL, sample = NULL, pval = 0.1){
  # Characterizing malignant cells in tumors with low T cell infiltration.
  print(paste("Characterizing",cellA.name,"cells in tumors with low",cellB.name,"infiltration."))
  cellA.markers<-sort(intersect(rB$genes,cellA.markers))
  results<-list(cellA.markers = cellA.markers,
                cellB.markers = cellB.markers,
                bulk.confounders = bulk.confounders,
                sc.confounders = sc.confounders)
  
  print("1. Estimating cell type B abundance in bulk gene expression")
  results$bulk.cellB.abn <- round(get.OE.bulk(rB,gene.sign = list(cellB.markers),num.rounds = 1000),2)
  rownames(results$bulk.cellB.abn)<-rB$samples ## OK
  tmp = readRDS(paste("Results/Results_bulk.cellB.abn_", sample, "_1123.rds", sep = ""))
  results$bulk.cellB.abn = tmp

  print("2. Looking for genes correlated with cell type B abundance in bulk gene expression")
  if(is.null(bulk.confounders)){
    results$bulk.cor<-get.cor(t(rB$tpm),results$bulk.cellB.abn,method = "pearson")
  }else{
    # bulk.confounders<-as.matrix(bulk.confounders)
    # b<-!is.na(rowSums(bulk.confounders))
    # print(paste("Using",sum(b),"bulk samples (with full data)."))
    # bulk.confounders<-bulk.confounders[b,]
    # results$bulk.cor<-pcor.mat(t(rB$tpm[,b]),results$bulk.cellB.abn[b],bulk.confounders,method = "pearson")
  }
  results$bulk.cor<-add.onesided.p(results$bulk.cor)
  
  print("3. Getting the seed signatures")
  results$seed.sig<-get.top.elements(results$bulk.cor[cellA.markers,c("p.pos","p.neg")],no.elm = 20,min.cf = pval)
  names(results$seed.sig)<-c("att","exc")
  print(summary(results$seed.sig))
  saveRDS(results$seed.sig, paste("Results/Results_seed.sig_", sample, "_1123.rds", sep = ""))
  
  print("4. Testing if the seed signatures are anti correlated")
  results$sc.seed.scores<-round(get.OE.sc(r.sc,results$seed.sig,num.rounds = 1000),2)
  cor.plot(results$sc.seed.scores,main = "Seed overall expression",
           ylab = "Seed exclusion (up)",xlab = "Seed exclusion (down)")
  
  print("5. Expanding the seed signatures")
  r.sc$q<-cbind(r.sc$comp,log(r.sc$comp.reads))
  # if(!is.null(sc.confounders)){
  #   r.sc$q<-cbind(r.sc$q,sc.confounders)
  #   results$sc.seed.scores.rgr<-t(get.residuals(t(results$sc.seed.scores),sc.confounders))
  #   cor.plot(results$sc.seed.scores.rgr,main = "Seed residuals",
  #            ylab = "Seed exclusion (up)",xlab = "Seed exclusion (down)")
  # }
  results$att.sc<-pcor.mat(t(r.sc$tpm),results$sc.seed.scores[,1],r.sc$q,method = "spearman")
  results$exc.sc<-pcor.mat(t(r.sc$tpm),results$sc.seed.scores[,2],r.sc$q,method = "spearman")
  results$att.sc<-add.onesided.p(results$att.sc) ## Inf
  results$exc.sc<-add.onesided.p(results$exc.sc)
  
  print("6. Generating the final signatures")
  f<-function(s){
    names(s)<-gsub("p.","",names(s),fixed = T)
    s$exc<-intersect(s$att.neg,s$exc.pos)
    s$att<-intersect(s$att.pos,s$exc.neg) # No att data -> results X.
    return(s)
  }
  b.rp<-!startsWith(r.sc$genes,"RP")
  b.seed<-is.element(r.sc$genes,unlist(results$seed.sig)) # True : 20
  s1<-c(get.top.elements(m = results$att.sc[,c("p.pos","p.neg")],no.elm = 200,min.cf = 0.01,main = "att"),
        get.top.elements(m = results$exc.sc[,c("p.pos","p.neg")],no.elm = 200,min.cf = 0.01,main = "exc"))
  s2<-c(get.top.elements(m = results$att.sc[b.rp,c("p.pos","p.neg")],no.elm = 200,min.cf = 0.01,main = "att"),
        get.top.elements(m = results$exc.sc[b.rp,c("p.pos","p.neg")],no.elm = 200,min.cf = 0.01,main = "exc"))
  s3<-c(get.top.elements(m = results$att.sc[b.seed,c("p.pos","p.neg")],no.elm = 20,min.cf = 0.01,main = "att"),
        get.top.elements(m = results$exc.sc[b.seed,c("p.pos","p.neg")],no.elm = 20,min.cf = 0.01,main = "exc"))
  s1<-f(s1);s2<-f(s2);s3<-f(s3)
  results$sigA<-s1
  results$sig.no.rp<-s1
  results$sig<-lapply(names(s1),function(x) sort(unique(c(s1[[x]],s2[[x]],s3[[x]]))))
  names(results$sig)<-names(s1)
  print(summary(results$sig[c("exc","att")]))
  
  results$sig.final<-c(results$seed.sig[c("exc","att")],results$sig[c("exc","att")])
  names(results$sig.final)<-c("exc.seed.up","exc.seed.down","exc.up","exc.down")
  print(summary(results$sig.final))
  
  if(!is.null(fileName)){
    saveRDS(results,file = paste0("~/Dropbox/ANGPT2_SingleCell/Results/",fileName, "_", sample,"_1123.rds"))
  }
  return(results)
}


# bulk.fpkm = as.data.frame(bulk.met@assays@data@listData[["HTSeq - FPKM"]])
# 
# ### DATA -> should use xena data instead of GDC
# bulk.xena = read.table("Data/TCGA.SKCM.sampleMap_HiSeqV2_PANCAN", header = T)
# rownames(bulk.xena) = bulk.xena$sample
# bulk.xena.met = bulk.xena[,colnames(bulk.met)]
# 
# FPKMtoTPM <- function(x) {
#  return(exp(log(x) - log(sum(x)) + log(1e6)))
# }
# bulk.tpm = Matrix(as.matrix(bulk.fpkm %>% mutate_if(is.numeric, FPKMtoTPM)), sparse = TRUE)
# bulk.fpkm = Matrix(as.matrix(bulk.fpkm, sparse = TRUE))
# 
# colnames(bulk.tpm)= colnames(bulk.met)
# rownames(bulk.tpm) = bulk.met@rowRanges$external_gene_name
# colnames(bulk.fpkm)= colnames(bulk.met)
# rownames(bulk.fpkm) = bulk.met@rowRanges$external_gene_name


## Metastasis TCGA-SKCM
sc.tpm = read.csv("Data/GSE115978/GSE115978_tpm.csv") # dim 23686 7187
sc.meta = read.csv("ImmuneResistance/Data/tumors.mal_tsne_anno.txt", sep = "\t")
rownames(sc.tpm) = sc.tpm$X
sc.tpm = sc.tpm[, colnames(sc.tpm) %in% sc.meta$NAME] # dim 23686 1881

r.sc<-list(cells = colnames(sc.tpm),genes = rownames(sc.tpm),tpm = sc.tpm,
           comp = colSums(sc.tpm>0),comp.reads = length(colnames(sc.tpm)))

gene = intersect(rownames(bulk.xena.met), r.sc[["genes"]])
bulk.xena.met = bulk.xena.met[gene,]
rB.tpm<-list(samples = colnames(bulk.xena.met),genes = gene,tpm = bulk.xena.met)

load("ImmuneResistance/Results/CellTypes/cell.type.sig.RData")

results <-cell.cell.intr(rB.tpm,r.sc,
                        cellA.markers = cell.sig$Mal, # 389
                        cellB.markers = cell.sig$T.CD8, # 50
                        cellA.name = "malignant",
                        cellB.name = "T.CD8",
                        fileName = "FINAL_Tcell_Exclusion", 
                        sample = "Met", pval = 0.1)


## Primary TCGA-SKCM

r.sc<-list(cells = colnames(sc),genes = rownames(sc),tpm = sc@assays$RNA@counts,
           comp = colSums(sc@assays$RNA@counts>0),comp.reads = length(colnames(sc)))

gene = intersect(rownames(bulk.xena.pri), r.sc[["genes"]])
bulk.xena.pri = bulk.xena.pri[gene,]
rB.tpm<-list(samples = colnames(bulk.xena.pri),genes = gene,tpm = bulk.xena.pri)

load("ImmuneResistance/Results/CellTypes/cell.type.sig.RData")

results <-cell.cell.intr(rB.tpm,r.sc,
                         cellA.markers = cell.sig$Mal, # 389
                         cellB.markers = cell.sig$T.CD8, # 50
                         cellA.name = "malignant",
                         cellB.name = "T.CD8",
                         fileName = "FINAL_Tcell_Exclusion",
                         sample = "Pri", pval = 0.1)



## Total TCGA-SKCM

r.sc<-list(cells = colnames(sc),genes = rownames(sc),tpm = sc@assays$RNA@counts,
           comp = colSums(sc@assays$RNA@counts>0),comp.reads = length(colnames(sc)))

gene = intersect(rownames(bulk.xena), r.sc[["genes"]])
bulk.xena = bulk.xena[gene,2:ncol(bulk.xena)] 
rB.tpm<-list(samples = colnames(bulk.xena),genes = gene,tpm = bulk.xena)

results <-cell.cell.intr(rB.tpm,r.sc,
                         cellA.markers = cell.sig$Mal, # 389
                         cellB.markers = cell.sig$T.CD8, # 50
                         cellA.name = "malignant",
                         cellB.name = "T.CD8",
                         fileName = "FINAL_Tcell_Exclusion",
                         sample = "Total", pval = 0.1)
