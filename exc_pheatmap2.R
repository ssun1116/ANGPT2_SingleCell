library(tidyverse)
library(pheatmap)
library(WGCNA)
library(ggplot2) # ggplot() for plotting
library(dplyr) # data reformatting
library(tidyr) # data reformatting
library(stringr) # string manipulation
setwd("~/Dropbox/ANGPT2_SingleCell")

## TCGA Metastasis

# ANGPT2 Correlated genes
d = readRDS('bulk.metastasis.RData')
a = WGCNA::cor(d)
a1 = as.data.frame(a)
a2 = a1 %>% select(ANGPT2)