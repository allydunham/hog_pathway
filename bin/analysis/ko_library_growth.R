# Script performing analying ko library growth
setwd('~/Projects/hog/')
library(rlang)
library(tidyverse)
library(magrittr)
library(ggpubr)
library(gplots)
library(ggdendro)
library(dendextend)
library(plotly)

#### Functions ####
## Turn tibble into a matrix with rownames from a column
tbl_to_matrix <- function(x, row){
  mat <- as.matrix(select(x, -row))
  rownames(mat) <- pull(x, !!row)
  return(mat)
}
  

## Split ko_growth tibble into subtable of a strain giving a variable spread over gene or condition
# Does not apply sorting because ko_scores tsv comes sorted
split_strains <- function(str, var, tbl, col='gene', row='condition'){
  # Function to filter ko_growth to a given strain with a matrix of condition against 'var' (score or qvalue)
  ko <- filter(tbl, strain == str) %>%
    select(!!row, !!col, !!var) %>%
    spread(key = col, value = !!var)
  return(ko)
}

## calculate condition ko growth profile dendogram from a tibble of strain scores
get_dend <- function(x){
  mat <- as.matrix(select(x, -condition))
  rownames(mat) <- x$condition
  return(as.dendrogram(hclust(dist(mat))))
}

## Calculate correlation between pairs of genes or conditions (based on profile according to the other)
# requires a tibble with the dependant variable in rows labeled by column 'var'
get_cor <- function(x, cor_meth = 'pearson', cor_use = 'pairwise', upper_tri=TRUE, var='condition'){
  mat <- as.matrix(select(x, -!!var))
  rownames(mat) <- pull(x, !!var)
  cors <- cor(t(mat), method = cor_meth, use = cor_use)
  if (upper_tri){
    cors[lower.tri(cors, diag = TRUE)] <- NA
  }
  var1 <- paste0(var,'1')
  var2 <- paste0(var,'2')
  cors %<>% as_tibble(rownames=var1) %>%
    gather(key = !!var2, value = 'cor', -!!var1) %>%
    drop_na()
  return(cors)
}

#### Import and Process Data ####
genes <- readRDS('data/Rdata/gene_meta_all.rds')
essential <- readRDS('data/Rdata/essential_genes.rds')
essential_genes <- filter(essential, essential == 'E') %>% pull(locus)
non_essential_genes <- filter(essential, essential == 'NE') %>% pull(locus)

ko_growth <- read_tsv('data/raw/ko_scores.txt', col_names = TRUE) %>%
  filter(!gene == 'WT') %>%
  filter(!duplicated(.[,c('strain', 'condition', 'gene')])) %>%
  select(-position) %>%
  mutate(condition = gsub('  ', ' ', condition)) %>% # Some conditions have double spaces in names
  mutate(name = if_else(is.na(name), gene, name))

strain_combs <- combn(c('S288C', 'UWOP', 'Y55', 'YPS'), 2)
strain_ko_scores <- sapply(unique(ko_growth$strain), split_strains, var='score', tbl=ko_growth, simplify = FALSE)
strain_ko_qvalues <- sapply(unique(ko_growth$strain), split_strains, var='qvalue', tbl=ko_growth, simplify = FALSE)

#### Analysis ####
## Dendograms
strain_con_dends <- lapply(strain_ko_scores, get_dend)

p_dend_list <- lapply(strain_con_dends, ggdendrogram, rotate=TRUE)
p_dends <- ggarrange(plotlist = p_dend_list, ncol = 2, nrow = 2, labels = names(strain_con_dends))
ggsave('figures/ko_growth/strain_condition_dendograms.pdf', p_dends, width = 10, height = 10)

# dendogram comparisons
# must filter Maltose 2%  (72H) since no data is available for Y55
strain_ko_scores_filtered <- lapply(unique(ko_growth$strain), split_strains, var='score', tbl=filter(ko_growth, !condition == 'Maltose 2% (72H)'))
names(strain_ko_scores_filtered) <- unique(ko_growth$strain)
strain_con_dends_filtered <- lapply(strain_ko_scores_filtered, get_dend)

pdf('figures/ko_growth/strain_condition_comparison_dendograms.pdf', width = 15, height = 20)
layout(matrix(1:(6*3), nrow = 3, byrow = TRUE), widths = rep(c(2,0.75,2), 2))
for (c in 1:dim(strain_combs)[2]){
  strains <- strain_combs[,c]
  tanglegram(untangle(strain_con_dends_filtered[[strains[1]]], strain_con_dends_filtered[[strains[2]]]),
             margin_inner = 15, common_subtrees_color_branches = TRUE,
             main_left = strains[1], main_right = strains[2], lwd=2, edge.lwd=2, just_one = FALSE)
}
layout(matrix(1, nrow = 1))
dev.off()

## Correlations
strain_con_cors <- bind_rows(lapply(strain_ko_scores, get_cor), .id = 'strain') %>%
  arrange(strain, condition1, condition2) %>%
  spread(strain, cor) %>%
  unite(col = 'pair', condition1, condition2, remove = FALSE)

p_strain_con_cor_cors <- mapply(function(x, y){
  ggplot(strain_con_cors, aes_string(x=x, y=y, text='pair')) + 
    geom_point(size=0.5) + 
    geom_abline(intercept = 0, slope = 1, colour='firebrick2') +
    geom_abline(intercept = -0.2, slope = 1, colour='firebrick2', linetype=2) +
    geom_abline(intercept = 0.2, slope = 1, colour='firebrick2', linetype=2)
  }, 
  strain_combs[1,], strain_combs[2,], SIMPLIFY = FALSE)
names(p_strain_con_cor_cors) <- paste(strain_combs[1,], strain_combs[2,], sep='_')

# removed rows correspond to pairs including Maltose 2%  (72H), since no data is available for Y55
p_strain_con_cor_cors_arr <- ggarrange(plotlist = p_strain_con_cor_cors, ncol = 3, nrow = 2, common.legend = TRUE)
ggsave('figures/ko_growth/strain_condition_condition_cors.pdf', p_strain_con_cor_cors_arr, width = 7, height = 5)

## Interactive plot
# Single
ggplotly(p_strain_con_cor_cors$UWOP_YPS, tooltip = c('pair'))

# Multiple
plotly_strain_cors <- subplot(sapply(p_strain_con_cor_cors, function(x){ggplotly(x, tooltip = c('pair'))}, simplify = FALSE),
                              nrows = 2, titleX = TRUE, titleY = TRUE, margin = c(0.05, 0.05, 0.07, 0.07))

# 3D Plots - show some joint
p_strain_cors_3d <- plot_ly(strain_con_cors, x=~S288C, y=~UWOP, z=~Y55, mode = 'markers', text=~pair) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'S288C'),
                      yaxis = list(title = 'UWOP'),
                      zaxis = list(title = 'Y55')))



## Per gene correlation
gene_ko_profile <- ko_growth %>%
#  mutate(score = if_else(qvalue < 0.9, score, 0)) %>%
  select(-qvalue, -gene) %>%
  rename(gene = name) %>%
  spread(key = strain, value = score)
  
gene_ko_profile_cor <- group_by(gene_ko_profile_cor, gene) %>%
  summarise(UWOP_S288C = cor(UWOP, S288C, use = 'na.or.complete'),
            UWOP_Y55 = cor(UWOP, Y55, use = 'na.or.complete'),
            UWOP_YPS = cor(UWOP, YPS, use = 'na.or.complete'),
            S288C_Y55 = cor(S288C, Y55, use = 'na.or.complete'),
            S288C_YPS = cor(S288C, YPS, use = 'na.or.complete'),
            YPS_Y55 = cor(YPS, Y55, use = 'na.or.complete')) %>%
  mutate(mean_cor = rowMeans(select(., UWOP_S288C, UWOP_Y55, UWOP_YPS, S288C_Y55, S288C_YPS, YPS_Y55), na.rm = TRUE), 
         var_cor = apply(select(., UWOP_S288C, UWOP_Y55, UWOP_YPS, S288C_Y55, S288C_YPS, YPS_Y55), 1, var, na.rm = TRUE))
  
p_gene_strain_cor_boxes <- ggplot(gather(gene_ko_profile_cor, key = 'strain_pair', value = 'cor', -gene), aes(x=strain_pair, y=cor, text=gene)) +
  geom_jitter()
ggplotly(p_gene_strain_cor_boxes, tooltip = c('gene'))

p_gene_mean_var_cor <- ggplot(gene_ko_profile_cor, aes(x=var_cor, y=mean_cor, label=gene)) +
  geom_point() +
  xlab('Variance of gene KO s-score correlation accross strains') +
  ylab('Mean gene KO s-score correlation accross strains')
ggplotly(p_gene_mean_var_cor)

# gene pair cor
# Only use genes which are significant in some condition

sig_genes <- read_tsv('data/raw/ko_scores.txt', col_names = TRUE) %>%
  filter(qvalue < 0.001) %>%
  filter(!gene == 'WT') %>%
  filter(!duplicated(.[,c('strain', 'condition', 'gene')])) %>%
  select(-position) %>%
  unite(comb, score, qvalue) %>%
  spread(key = strain, value = comb) %>%
  mutate(num_strains = (!is.na(S288C)) + (!is.na(UWOP)) + (!is.na(Y55)) + (!is.na(YPS))) %>%
  filter(num_strains > 2) %>%
  pull(gene) %>%
  unique()

strain_gene_cors <- bind_rows(lapply(sapply(unique(ko_growth$strain), split_strains, var='score', tbl=filter(ko_growth, gene %in% sig_genes), 
                                            col='condition', row='name', simplify = FALSE), 
                                     get_cor, var='name'), .id = 'strain') %>%
  rename(gene1 = name1, gene2 = name2) %>%
  arrange(strain, gene1, gene2) %>%
  spread(strain, cor) %>%
  unite(col = 'pair', gene1, gene2, remove = FALSE)

p_strain_gene_gene_cors <- mapply(function(x, y){
  ggplot(strain_gene_cors, aes_string(x=x, y=y, text='pair')) + 
    geom_point(size=0.5) + 
    geom_abline(intercept = 0, slope = 1, colour='firebrick2') +
    geom_abline(intercept = -0.2, slope = 1, colour='firebrick2', linetype=2) +
    geom_abline(intercept = 0.2, slope = 1, colour='firebrick2', linetype=2)
}, 
strain_combs[1,], strain_combs[2,], SIMPLIFY = FALSE)
names(p_strain_gene_gene_cors) <- paste(strain_combs[1,], strain_combs[2,], sep='_')

p_strain_gene_gene_cors_arr <- ggarrange(plotlist = p_strain_gene_gene_cors, ncol = 3, nrow = 2, common.legend = TRUE)
ggsave('figures/ko_growth/strain_gene_gene_cors.pdf', p_strain_gene_gene_cors_arr, width = 7, height = 5)

plotly_strain_gene_cors <- subplot(sapply(p_strain_gene_gene_cors, function(x){ggplotly(x, tooltip = c('pair'))}, simplify = FALSE),
                                   nrows = 2, titleX = TRUE, titleY = TRUE, margin = c(0.05, 0.05, 0.07, 0.07))

strain_gene_cors_all <- bind_rows(lapply(sapply(unique(ko_growth$strain), split_strains, var='score', tbl=ko_growth, 
                                            col='condition', row='name', simplify = FALSE), 
                                     get_cor, var='name'), .id = 'strain') %>%
  rename(gene1 = name1, gene2 = name2) %>%
  arrange(strain, gene1, gene2) %>%
  spread(strain, cor) %>%
  unite(col = 'pair', gene1, gene2, remove = FALSE) %>%
  mutate(mean = rowMeans(select(., S288C, UWOP, Y55, YPS), na.rm = TRUE))

strain_gene_cors_all_mat <- select(strain_gene_cors_all, gene1, gene2, S288C) %>%
  spread(key = 'gene2', value = 'S288C') %>%
  tbl_to_matrix(., 'gene1')

pdf('figures/ko_growth/s288c_gene_cor_heatmap.pdf', width = 20, height = 20)
heatmap.2(strain_gene_cors_all_mat, symm = TRUE, revC = TRUE)
dev.off()
