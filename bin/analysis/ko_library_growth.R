# Script performing analying ko library growth
setwd('~/Projects/hog/')
library(rlang)
library(tidyverse)
library(magrittr)
library(ggpubr)
library(gplots)
library(ggdendro)
library(dendextend)

#### Import and Process Data ####
genes <- readRDS('data/Rdata/gene_meta_all.rds')
essential <- readRDS('data/Rdata/essential_genes.rds')
essential_genes <- filter(essential, essential == 'E') %>% pull(locus)
non_essential_genes <- filter(essential, essential == 'NE') %>% pull(locus)

ko_growth <- read_tsv('data/raw/ko_scores.txt', col_names = TRUE) %>%
  filter(!gene == 'WT') %>%
  filter(!duplicated(.[,c('strain', 'condition', 'gene')])) %>%
  select(-position)

# Does not apply sorting because ko_scores tsv comes sorted
split_strains <- function(str, var){
  # Function to filter ko_growth to a given strain with a matrix of condition against 'var' (score or qvalue)
  ko <- filter(ko_growth, strain == str) %>%
    select(condition, gene, !!var) %>%
    spread(key = gene, value = !!var)
  return(ko)
}

strain_ko_scores <- lapply(unique(ko_growth$strain), split_strains, var='score')
names(strain_ko_scores) <- unique(ko_growth$strain)

strain_ko_qvalues <- lapply(unique(ko_growth$strain), split_strains, var='qvalue')
names(strain_ko_qvalues) <- unique(ko_growth$strain)

#### Analysis ####
## Dendograms
get_dend <- function(x){
  mat <- as.matrix(select(x, -condition))
  rownames(mat) <- x$condition
  return(as.dendrogram(hclust(dist(mat))))
}

strain_con_dends <- lapply(strain_ko_scores, get_dend)

p_dend_list <- lapply(strain_con_dends, ggdendrogram, rotate=TRUE)
p_dends <- ggarrange(plotlist = p_dend_list, ncol = 2, nrow = 2, labels = names(strain_con_dends))
ggsave('figures/ko_growth/strain_condition_dendograms.pdf', p_dends, width = 10, height = 10)

## Correlations
# 
get_cor <- function(x){
  mat <- as.matrix(select(x, -condition))
  rownames(mat) <- x$condition
  cors <- cor(t(mat), use = 'complete')
  cors[lower.tri(cors, diag = TRUE)] <- NA
  cors %<>% as_tibble(rownames='condition1') %>%
    gather(key = 'condition2', value = 'cor', -condition1) %>%
    filter(!is.na(cor))
  return(cors)
}

strain_con_cors <- bind_rows(lapply(strain_ko_scores, get_cor), .id = 'strain') %>%
  arrange(strain, condition1, condition2) %>%
  spread(strain, cor)

p_test <- ggplot(strain_con_cors, aes(x=S288C, y=UWOP)) +
  geom_point()





