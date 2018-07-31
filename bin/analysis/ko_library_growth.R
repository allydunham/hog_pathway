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
library(Rtsne)

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

# Determine which level list a gene is in
# sig_list should have the most desirable category last if there are overlaps
get_sig <- function(item, sig_list){
  for (i in length(sig_list):1){
    if (item %in% sig_list[[i]]){
      return(i)
    }
  }
  return(0)
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

# Get gens sig in at least n conditions
sig_genes <- read_tsv('data/raw/ko_scores.txt', col_names = TRUE) %>%
  mutate(name = if_else(is.na(name), gene, name)) %>%
  filter(qvalue < 0.01) %>%
  filter(!gene == 'WT') %>%
  filter(!duplicated(.[,c('strain', 'condition', 'gene')])) %>%
  select(-position) %>%
  unite(comb, score, qvalue) %>%
  spread(key = strain, value = comb) %>%
  mutate(num_strains = (!is.na(S288C)) + (!is.na(UWOP)) + (!is.na(Y55)) + (!is.na(YPS)))

sig_genes_4 <- sig_genes %>%
  filter(num_strains > 3) %>%
  pull(name) %>%
  unique()

sig_genes_3 <- sig_genes %>%
  filter(num_strains > 2) %>%
  pull(name) %>%
  unique()

sig_genes_2 <- sig_genes %>%
  filter(num_strains > 1) %>%
  pull(name) %>%
  unique()

sig_genes_1 <- sig_genes %>%
  filter(num_strains > 0) %>%
  pull(name) %>%
  unique()

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
  spread(key = strain, value = score) %>%
  mutate(num_strains = rowSums(abs(select(., S288C, UWOP, YPS, Y55)) > 1, na.rm = TRUE))
  
gene_ko_profile_cor <- group_by(gene_ko_profile, gene) %>%
  summarise(UWOP_S288C = cor(UWOP, S288C, use = 'na.or.complete'),
            UWOP_Y55 = cor(UWOP, Y55, use = 'na.or.complete'),
            UWOP_YPS = cor(UWOP, YPS, use = 'na.or.complete'),
            S288C_Y55 = cor(S288C, Y55, use = 'na.or.complete'),
            S288C_YPS = cor(S288C, YPS, use = 'na.or.complete'),
            YPS_Y55 = cor(YPS, Y55, use = 'na.or.complete')) %>%
  mutate(mean_cor = rowMeans(select(., UWOP_S288C:YPS_Y55), na.rm = TRUE),
         max_cor = apply(select(., UWOP_S288C:YPS_Y55), 1, function(x){ifelse(all(is.na(x)), NA, unname(x)[which.max(abs(x))])}),
         abs_mean_cor = rowMeans(abs(select(., UWOP_S288C:YPS_Y55)), na.rm = TRUE),
         var_cor = apply(select(., UWOP_S288C:YPS_Y55), 1, var, na.rm = TRUE),
         concord = apply(sign(select(., UWOP_S288C:YPS_Y55)), 1, function(x){ifelse(all(is.na(x)), NA, max(x,na.rm=TRUE)==min(x,na.rm=TRUE))}),
         num_pos = rowSums(select(., UWOP_S288C:YPS_Y55) > 0, na.rm = TRUE),
         which_max = apply(select(., UWOP_S288C:YPS_Y55), 1, function(x){ifelse(all(is.na(x)), NA, names(x)[which.max(abs(x))])}),
         which_min = apply(select(., UWOP_S288C:YPS_Y55), 1, function(x){ifelse(all(is.na(x)), NA, names(x)[which.min(x)])}),
         sig_level = as.factor(sapply(gene, get_sig, sig_list=list(sig_genes_1, sig_genes_2, sig_genes_3, sig_genes_4))))
  
p_gene_strain_cor_boxes <- ggplot(gather(gene_ko_profile_cor, key = 'strain_pair', value = 'cor', UWOP_S288C:YPS_Y55), aes(x=strain_pair, y=cor, text=gene)) +
  geom_jitter()
ggplotly(p_gene_strain_cor_boxes, tooltip = c('gene'))

plot_ly(gene_ko_profile_cor, x=~var_cor, y=~mean_cor, color=~sig_level, text=~gene) %>%
  add_markers() %>%
  layout(title = "Summary of gene profiles correlations between strains",
         xaxis = list(title = 'Variance of gene KO s-score correlation accross strains'),
         yaxis = list(title = 'Mean gene KO s-score correlation accross strains'))

p_gene_mean_var_cor <- ggarrange(ggplot(gene_ko_profile_cor, aes(x=var_cor, y=mean_cor, label=gene, colour=sig_level)) +
                                   geom_point() +
                                   xlab('') +
                                   ylab('Mean Correlation') +
                                   guides(colour=guide_legend(title = 'Conditionally Significant Strains')), 
                                 ggplot(gene_ko_profile_cor, aes(x=var_cor, y=mean_cor, label=gene, colour=sig_level)) +
                                   geom_density_2d() +
                                   xlab('') +
                                   ylab('') +
                                   guides(colour=guide_legend(title = 'Conditionally Significant Strains')),
                                 ggplot(gene_ko_profile_cor, aes(x=var_cor, y=abs_mean_cor, label=gene, colour=sig_level)) +
                                   geom_point() +
                                   xlab('Correlation Variance') +
                                   ylab('Mean Absolute Correlation') +
                                   guides(colour=guide_legend(title = 'Conditionally Significant Strains')), 
                                 ggplot(gene_ko_profile_cor, aes(x=var_cor, y=abs_mean_cor, label=gene, colour=sig_level)) +
                                   geom_density_2d() +
                                   xlab('Correlation Variance') +
                                   ylab('') +
                                   guides(colour=guide_legend(title = 'Conditionally Significant Strains')),
                                 nrow = 2,
                                 ncol = 2,
                                 common.legend = TRUE,
                                 legend = 'top',
                                 align = 'hv')
ggsave('figures/ko_growth/gene_strain_impact_profile_correlations.pdf', p_gene_mean_var_cor, width = 9, height = 9)

p_gene_strain_cor_density <- ggplot(gather(gene_ko_profile_cor, key = 'strain_pair', value = 'cor', UWOP_S288C:YPS_Y55), aes(colour=strain_pair, x=cor)) +
  geom_density() +
  xlab('Strain-Strain Gene KO Profile Correlation') +
  ylab('Density') +
  guides(colour = guide_legend(title = 'Strain Pair'))
ggsave('figures/ko_growth/strain_strain_gene_cor_density.pdf', p_gene_strain_cor_density, width = 7, height = 5)

## Meaningful? Genes that show concordance tend to have higher absolute correlation values
p_gene_cor_concordance <- ggplot(filter(gene_ko_profile_cor, !is.na(concord)), aes(x=concord, y=abs_mean_cor, colour=concord)) +
  geom_boxplot() +
  stat_compare_means(comparisons = list(c('TRUE','FALSE'))) +
  stat_summary(geom = 'text', fun.data = function(x){return(c(y = -0.05, label = length(x)))}) +
  stat_summary(geom = 'text', aes(colour=concord), fun.data = function(x){return(c(y = mean(x), label = signif(mean(x), 3)))}) +
  guides(colour=FALSE) +
  xlab('Strain-Strain corelations have the same sign?') +
  ylab('Mean Absolute Correlation')
ggsave('figures/ko_growth/abs_gene_profile_correlation_concord_box.pdf', p_gene_cor_concordance, width = 7, height = 5)

p_gene_cor_num_pos <- ggplot(filter(gene_ko_profile_cor, !is.na(concord)), aes(x=num_pos, y=abs_mean_cor, group=cut_width(num_pos, width = 1))) +
  geom_boxplot() +
  stat_summary(geom = 'text', fun.data = function(x){return(c(y = -0.05, label = length(x)))}) +
  stat_summary(geom = 'text', colour='firebrick2', fun.data = function(x){return(c(y = mean(x), label = signif(mean(x), 3)))}) +
  guides(colour=FALSE) +
  xlab('Number of positive strain-strain correlations') +
  ylab('Mean Absolute Correlation')
ggsave('figures/ko_growth/abs_gene_profile_correlation_num_pos.pdf', p_gene_cor_num_pos, width = 7, height = 5)


plot_ly(data = gene_ko_profile_cor, x=~concord, y=~abs_mean_cor, text=~gene, color=~concord, type='box') %>%
  add_markers() %>%
  layout(title = "Relationship between correlation sign agreement and absolute correlation",
         xaxis = list(title = 'Gene Corrrelation Sign Agreement'),
         yaxis = list(title = 'Mean Absolute Correlation'))

high_cor_concordant_genes <- gene_ko_profile_cor %>%
  filter(concord == TRUE, abs_mean_cor > 0.3) %>%
  select(gene:var_cor) %>%
  mutate(id = structure(genes$id, names = genes$name)[gene])
write_tsv(high_cor_concordant_genes, 'data/high_strain_strain_cor_concord_genes.tsv')

# Filters tried - mean_cor > or < 0.3, in sig_genes_1/2/3
gene_ko_profile_cor_strong <- filter(gene_ko_profile_cor, sig_level == 3) %>% drop_na(-var_cor)

gene_strain_cor_tsne <- Rtsne(tbl_to_matrix(select(gene_ko_profile_cor_strong, gene, UWOP_S288C, UWOP_Y55, UWOP_YPS, S288C_Y55, S288C_YPS, YPS_Y55), 'gene'))

gene_strain_cor_tsne_df <- as_tibble(gene_strain_cor_tsne$Y) %>%
  rename(tsne1 = V1, tsne2 = V2) %>%
  mutate(gene = gene_ko_profile_cor_strong$gene, # These vars are mainly used to colour points
         mean_cor = gene_ko_profile_cor_strong$abs_mean_cor,
         var_cor = gene_ko_profile_cor_strong$var_cor,
         num_pos = gene_ko_profile_cor_strong$num_pos,
         concord = gene_ko_profile_cor_strong$concord,
         which_max = gene_ko_profile_cor_strong$which_max,
         which_min = gene_ko_profile_cor_strong$which_min,
         max_mag = gene_ko_profile_cor_strong$max_cor)
  
plot_ly(gene_strain_cor_tsne_df, x=~tsne1, y=~tsne2, color=~mean_cor, text=~gene, mode='markers') %>%
  add_markers() %>%
  layout(title = "tSNE of Gene Profile Correlation Between Strains",
         xaxis = list(title = 'tSNE 1'),
         yaxis = list(title = 'tSNE 2'))

# Only use genes which are significant in 3 conditions
strain_gene_cors <- bind_rows(lapply(sapply(unique(ko_growth$strain), split_strains, var='score', tbl=filter(ko_growth, name %in% sig_genes_3), 
                                            col='condition', row='name', simplify = FALSE), 
                                     get_cor, var='name'), .id = 'strain') %>%
  rename(gene1 = name1, gene2 = name2) %>%
  arrange(strain, gene1, gene2) %>%
  spread(strain, cor) %>%
  unite(col = 'pair', gene1, gene2, remove = FALSE)

p_strain_gene_gene_cors <- mapply(function(x, y){
  ggplot(strain_gene_cors, aes_string(x=x, y=y)) + 
    geom_bin2d(bins=40) + 
    geom_abline(intercept = 0, slope = 1, colour='firebrick2') +
    geom_abline(intercept = -0.2, slope = 1, colour='firebrick2', linetype=2) +
    geom_abline(intercept = 0.2, slope = 1, colour='firebrick2', linetype=2)
}, 
strain_combs[1,], strain_combs[2,], SIMPLIFY = FALSE)
names(p_strain_gene_gene_cors) <- paste(strain_combs[1,], strain_combs[2,], sep='_')

p_strain_gene_gene_cors_arr <- ggarrange(plotlist = p_strain_gene_gene_cors, ncol = 3, nrow = 2, common.legend = TRUE)
ggsave('figures/ko_growth/strain_sig_gene_gene_cors_dens.pdf', p_strain_gene_gene_cors_arr, width = 7, height = 5)

plotly_strain_gene_cors <- subplot(sapply(p_strain_gene_gene_cors, function(x){ggplotly(x, tooltip = c('pair'))}, simplify = FALSE),
                                   nrows = 2, titleX = TRUE, titleY = TRUE, margin = c(0.05, 0.05, 0.07, 0.07))

## Genes not required in one
odd_genes <- readRDS('data/Rdata/gene_ko_growth_full.rds') %>%
  mutate(name = if_else(is.na(name), gene, name)) %>%
  drop_na() %>%
  mutate(sig_strains = (S288C_qvalue < 0.05) + (UWOP_qvalue < 0.05) + (Y55_qvalue < 0.05) + (YPS_qvalue < 0.05)) %>%
  filter(sig_strains == 3)

