library(tidyverse)
library(pheatmap)
library(WGCNA)

d = readRDS('bulk.metastasis.RData')
a = WGCNA::cor(d)
a1 = as.data.frame(a)

a2 = a1 %>% select(ANGPT2)