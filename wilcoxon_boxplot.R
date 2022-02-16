# load packages
library(ggplot2) # ggplot() for plotting
library(dplyr) # data reformatting
library(ggpubr)
library(tidyverse)
library(cowplot)
#setwd("~/Dropbox/ANGPT2_SingleCell")

# Total bulk sample
bulk = as.data.frame(read.table("Data/TCGA.SKCM.sampleMap_HiSeqV2", header = T)) # 475 obs - first column : gene name
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

exc.up = c("FSCN1","TPM2", "H19", "BGN","HMGB1","EIF3H","DCAF13","LYPLA1","MEX3A","SOX4","TIMM50")
exc.down = c("PDE4DIP","PERP","ITGA3","NPC1","NAV2","LOC100126784","CTSD","IRF4","LGALS3","BHLHE41",
             "MMP14","TIMP2","IGSF8","ATP6V0C","IFI27L2","S100A13","S100A6","EMP1","LAMP2")

bulk.up = bulk2[,colnames(bulk2) %in% exc.up]
bulk.up = cbind.data.frame(rownames(bulk.up), bulk.up)
colnames(bulk.up)[1] = "sample"

bulk.down = bulk2[,colnames(bulk2) %in% exc.down]
bulk.down = cbind.data.frame(rownames(bulk.down), bulk.down)
colnames(bulk.down)[1] = "sample"

## Wilcox test
bulk.up$mean = rowMeans(bulk.up[,2:ncol(bulk.up)])
bulk.up$group = ifelse(bulk.up$sample %in% quant_tot.50, "<50th",
                       ifelse(bulk.up$sample %in% quant_tot.50.90, "50-90th", 
                              ifelse(bulk.up$sample %in% quant_tot.90, ">90th", NA)))

library(ggpubr)
my_comparisons <- list( c("<50th", "50-90th"), c("50-90th", ">90th"), c("<50th", ">90th"))
p1 <- ggboxplot(bulk.up, x = "group", y = "mean", color = "group", palette = "jco") +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")  +
  labs(title = "Exclusion-up", x = "ANGPT2 Groups", y = "Average z-score")

bulk.down$mean = rowMeans(bulk.down[,2:ncol(bulk.down)])
bulk.down$group = ifelse(bulk.down$sample %in% quant_tot.50, "<50th",
                         ifelse(bulk.down$sample %in% quant_tot.50.90, "50-90th", 
                                ifelse(bulk.down$sample %in% quant_tot.90, ">90th", NA)))

my_comparisons <- list( c("<50th", "50-90th"), c("50-90th", ">90th"), c("<50th", ">90th"))
p2 <- ggboxplot(bulk.down, x = "group", y = "mean", color = "group", palette = "jco") +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
  labs(title = "Exclusion-down", x = "ANGPT2 Groups", y = "Average z-score")

plot_row <- plot_grid(p1, p2)
# now add the title
title <- ggdraw() + 
  draw_label(
    "(L) T-cell exclusion pattern by ANGPT2 expression groups",
    fontface = 'bold',
    x = 0,
    hjust = 0
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 7)
  )
p <- plot_grid(
  title, plot_row,
  ncol = 1,
  # rel_heights values control vertical title margins
  rel_heights = c(0.1, 1)
)

ggsave(filename = "wilcoxon_boxplot.eps", plot = p, width = 8, height = 5)
ggsave(filename = "wilcoxon_boxplot.pdf", plot = p, width = 8, height = 5)
