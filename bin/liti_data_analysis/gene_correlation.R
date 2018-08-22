# Script to analyse whole genome P(Aff) scores and their correlations
#setwd('~/Projects/hog/')
library(MASS)
library(mixtools)
library(gplots)
library(tidyverse)
library(magrittr)
library(moments)
library(ggpubr)

#### Import Data ####
strains <- readRDS('data/Rdata/strain_meta.rds')
filtered_strains <- filter(strains, Ploidy == 2, Aneuploidies == 'euploid') %>% pull(`Standardized name`)

genes <- readRDS('data/Rdata/gene_meta_all.rds')
hog_genes <- readRDS('data/Rdata/gene_meta_hog.rds')

complexes <- readRDS('data/Rdata/complex_members.rds')

essential <- readRDS('data/Rdata/essential_genes.rds')
essential_hash <- structure(essential$essential, names=essential$locus)

growth <- readRDS('data/Rdata/growth_liti.rds') %>%
  filter(strain %in% filtered_strains)
  
probs <- readRDS('data/Rdata/paff_all_genes.rds') %>%
  filter(strain %in% filtered_strains, strain %in% growth$strain)

prob_mat <- readRDS('data/Rdata/paff_all_genes_mat.rds') %>%
  filter(strain %in% filtered_strains, strain %in% growth$strain)

cor_genes_melt <- readRDS('data/Rdata/all_gene_correlations.rds')
cor_growth <- readRDS('data/Rdata/gene_growth_correlations_matrix.rds')
cor_growth_melt <- readRDS('data/Rdata/gene_growth_correlations.rds')
gene_gene_cor <- readRDS('data/Rdata/all_gene_correlations_matrix.rds')

counts <- readRDS('data/Rdata/hog_gene_mut_counts.rds') %>%
  filter(strain %in% growth$strain)

#### Analysis ####
## Gene/Gene Correlation
pdf('figures/heatmaps/all_genes_paff_heatmap.pdf', width = 50, height = 50)
cols <- colorRampPalette(c("blue", "white","red"))(256)
heatmap.2(gene_gene_cor, symm = TRUE, revC = TRUE, col=cols,
          breaks=seq(-1,1,2/256), trace = "none")
dev.off()

## Gene/Growth Correlation
pdf('figures/heatmaps/all_genes_growth_heatmap.pdf', width = 50, height = 50)
cols <- colorRampPalette(c("blue", "white","red"))(256)
heatmap.2(cor_growth, col=cols, breaks=seq(-1,1,2/256), trace = "none", symkey = FALSE)
dev.off()

meanGrowthCor <- mean(cor_growth)

# normalise by subtraction
pdf('figures/heatmaps/all_genes_growth_heatmap_norm_sub.pdf', width = 50, height = 50)
cols <- colorRampPalette(c("blue", "white","red"))(256)
heatmap.2(cor_growth - meanGrowthCor, col=cols, breaks=seq(-1,1,2/256), trace = "none", symkey = FALSE)
dev.off()

# Normalise by division
pdf('figures/heatmaps/all_genes_growth_heatmap_norm_div.pdf', width = 50, height = 50)
cols <- colorRampPalette(c("blue", "white","red"))(256)
heatmap.2(cor_growth / meanGrowthCor, col=cols, trace = "none", symkey = FALSE)
dev.off()

# Hog Genes Only
cor_hog_counts <- cor(select(counts, -strain))
rownames(cor_hog_counts) <- structure(hog_genes$name, names=hog_genes$id)[rownames(cor_hog_counts)]

pdf('figures/heatmaps/hog_genes_mut_count_heatmap.pdf', width = 20, height = 20)
cols <- colorRampPalette(c("blue", "white","red"))(256)
heatmap.2(cor_hog_counts, symm = TRUE, revC = TRUE, col=cols,
          breaks=seq(-1,1,2/256), trace = "none")
dev.off()

# growth
cor_growth_hog <- cor_growth[rownames(cor_growth) %in% hog_genes$id,]
rownames(cor_growth_hog) <- structure(hog_genes$name, names=hog_genes$id)[rownames(cor_growth_hog)]
  
pdf('figures/heatmaps/hog_genes_growth_heatmap.pdf', width = 50, height = 50)
cols <- colorRampPalette(c("blue", "white","red"))(256)
heatmap.2(cor_growth_hog, breaks=seq(-1,1,2/256), 
          col=cols, trace = "none", symkey = FALSE, cexRow = 2.5, cexCol = 2.5, margins = c(24,15))
dev.off()

# Only select conditions
cor_growth_hog_osmotic <- cor_growth_hog[,c("ypdkcl2m","ypdnacl1m", "ypdnacl15m", "ypd14",
                                            "ypd40", "ypdcuso410mm", "ypdsodiummetaarsenite")]

pdf('figures/heatmaps/hog_genes_growth_osmotic_cons_heatmap.pdf', width = 15, height = 25)
cols <- colorRampPalette(c("blue", "white","red"))(256)
heatmap.2(cor_growth_hog_osmotic, breaks=seq(-1,1,2/256), 
          col=cols, trace = "none", symkey = FALSE, cexRow = 2.5, cexCol = 2.5, margins = c(30,15),
          labCol = c("KCl 2M", "NaCl 1M", "NaCl 1.5M", expression("14"~degree~"C"), expression("40"~degree~"C"), 
                                      expression("CuSO"[4]~10~"mM"), "Na Metaarsenite 2.5mM"))
dev.off()

# Against counts
cor_hog_growth_counts <- cor(select(counts, -strain), select(growth, -strain))
rownames(cor_hog_growth_counts) <- structure(hog_genes$name, names=hog_genes$id)[rownames(cor_hog_growth_counts)]

pdf('figures/heatmaps/hog_genes_mut_count_growth_heatmap.pdf', width = 20, height = 20)
cols <- colorRampPalette(c("blue", "white","red"))(256)
heatmap.2(cor_hog_growth_counts, col=cols, symkey=FALSE, cexRow = 2.5, cexCol = 2.5, margins = c(30,15),
          breaks=seq(-1,1,2/256), trace = "none")
dev.off()

