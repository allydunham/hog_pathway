# Script performing analying ko library growth
setwd('~/Projects/hog/')

## Load Packages and Import data
# Packages loaded before Tidyverse
library(gplots)
library(Rtsne)
library(EMT)
library(MASS)
library(e1071)

# Import data and load Tidyverse
source('bin/general_functions.R')
source('bin/ko_growth_analysis/ko_analysis_functions.R')
source('bin/ko_growth_analysis/load_ko_data.R')

# Other Packages
library(GGally)
library(ggpubr)
library(ggdendro)
library(dendextend)
library(plotly)
library(cowplot)
library(kSamples)

#######################################################
##################### Analysis ########################
#######################################################
######### By Condition ########
#### Condition Dendograms ####
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

#### Condition-Condition Correlations ####
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
# Multiple
plotly_strain_cors <- subplot(sapply(p_strain_con_cor_cors, function(x){ggplotly(x, tooltip = c('pair'))}, simplify = FALSE),
                              nrows = 2, titleX = TRUE, titleY = TRUE, margin = c(0.05, 0.05, 0.07, 0.07))

# 3D Plots - show some joint
p_strain_cors_3d <- plot_ly(strain_con_cors, x=~S288C, y=~UWOP, z=~Y55, mode = 'markers', text=~pair) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'S288C'),
                      yaxis = list(title = 'UWOP'),
                      zaxis = list(title = 'Y55')))

#### Conditions Clustering ####
gene_ko_profile <- ko_growth %>%
  #  mutate(score = if_else(qvalue < 0.9, score, 0)) %>%
  select(-qvalue, -gene) %>%
  rename(gene = name) %>%
  spread(key = strain, value = score) %>%
  mutate(num_strains = rowSums(abs(select(., S288C, UWOP, YPS, Y55)) > 1, na.rm = TRUE))

condition_profiles <- select(gene_ko_profile, -num_strains) %>%
  gather(key = 'strain', value = 'score', S288C:YPS) %>%
  spread(key = gene, value = score)

condition_profiles[is.na(condition_profiles)] <- 0

## Per strain condition tSNE
per_strain_con_tsne <- function(str, tbl){
  tbl <- filter(tbl, strain == !!str) %>% select(-strain) %>% tbl_to_matrix(., 'condition')
  tsne <- Rtsne(tbl, perplexity=10)$Y %>%
    as_tibble() %>%
    rename(tSNE1=V1, tSNE2=V2) %>%
    mutate(condition = rownames(tbl))
  return(tsne)
}
condition_tsnes_df <- sapply(c('S288C', 'UWOP', 'YPS', 'Y55'), per_strain_con_tsne,
                           tbl=select(condition_profiles, condition, strain, AAC1:ZWF1), simplify = FALSE) %>%
  bind_rows(.id='strain')

p_condition_tsne_by_strain <- ggplot(condition_tsnes_df, aes(x=tSNE1, y=tSNE2, colour=strain, text=condition)) +
  geom_point() +
  facet_wrap(~strain, scales = 'free') +
  guides(colour=FALSE)
ggplotly(p_condition_tsne_by_strain)

## Per strain condition PCA
per_strain_con_pca <- function(str, tbl){
  tbl <- filter(tbl, strain == !!str) %>% select(-strain) %>% tbl_to_matrix(., 'condition')
  pca <- prcomp(tbl)$x %>%
    as_tibble(rownames='condition')
  return(pca)
}
condition_pca_df <- sapply(c('S288C', 'UWOP', 'YPS', 'Y55'), per_strain_con_pca,
                           tbl=select(condition_profiles, condition, strain, AAC1:ZWF1), simplify = FALSE) %>%
  bind_rows(.id='strain')

p_condition_pca_by_strain <- ggplot(condition_pca_df, aes(x=PC5, y=PC6, colour=strain, text=condition)) +
  geom_point() +
  facet_wrap(~strain) +
  guides(colour=FALSE)
ggplotly(p_condition_pca_by_strain)

## Kmeans over all conditions
condition_profiles %<>% mutate(cluster = kmeans(select(condition_profiles, AAC1:ZWF1), centers = 8, nstart = 5)$cluster,
                               equiv_con = equiv_cons[condition])

p_con_clusters <- ggplot(condition_profiles, aes(x=cluster, fill=equiv_con, text=condition)) + geom_bar() + facet_wrap(~strain)
ggplotly(p_con_clusters)

## All strains condition gene profile tSNE
comb_tsne <- Rtsne(select(condition_profiles, condition, AAC1:ZWF1) %>% tbl_to_matrix(row = 'condition'))$Y %>%
  as_tibble() %>%
  rename(tSNE1=V1, tSNE2=V2) %>%
  mutate(strain = condition_profiles$strain,
         condition = condition_profiles$condition,
         equiv_con = condition_profiles$equiv_con,
         cluster = condition_profiles$cluster)

plot_ly(comb_tsne, x=~tSNE1, y=~tSNE2, color=~strain, text=~condition) %>% add_markers()

## All strains condition gene profile PCA
con_profile_pca <- prcomp(select(condition_profiles, condition, AAC1:ZWF1) %>% tbl_to_matrix(row = 'condition'))$x %>%
  as_tibble(rownames = 'condition') %>%
  mutate(strain = condition_profiles$strain,
         equiv_con = condition_profiles$equiv_con,
         cluster = condition_profiles$cluster)

subplot(plot_ly(con_profile_pca, x=~PC1, y=~PC2, color=~strain, text=~condition, type='scatter', mode='markers', legendgroup=~strain),
        plot_ly(con_profile_pca, x=~PC3, y=~PC4, color=~strain, text=~condition, type='scatter', mode='markers', legendgroup=~strain, showlegend=FALSE),
        plot_ly(con_profile_pca, x=~PC5, y=~PC6, color=~strain, text=~condition, type='scatter', mode='markers', legendgroup=~strain, showlegend=FALSE),
        plot_ly(con_profile_pca, x=~PC7, y=~PC8, color=~strain, text=~condition, type='scatter', mode='markers', legendgroup=~strain, showlegend=FALSE),
        nrows = 2, titleX = TRUE, titleY = TRUE)

#### Number of relavent genes for each strain in each condition ####
condition_sensitivity <- ko_growth %>%
  filter(qvalue < 0.05) %>%
  select(strain, condition, name, score) %>%
  group_by(condition) %>%
  spread(key = strain, value = score) %>%
  summarise(S288C_sig_gene_count = sum(!is.na(S288C)),
            UWOP_sig_gene_count = sum(!is.na(UWOP)),
            Y55_sig_gene_count = sum(!is.na(Y55)),
            YPS_sig_gene_count = sum(!is.na(YPS)))
write_tsv(condition_sensitivity, 'data/ko_condition_sensitivity.tsv')

p_con_sens_bar <- condition_sensitivity %>%
  rename(S288C = S288C_sig_gene_count,
         UWOP = UWOP_sig_gene_count,
         Y55 = Y55_sig_gene_count,
         YPS = YPS_sig_gene_count) %>%
  mutate(condition = factor(condition, levels = condition[order(.$S288C/(.$S288C + .$UWOP + .$Y55 + .$YPS))], ordered = TRUE)) %>%
  gather(key = 'strain', value = 'sig_count', S288C:YPS) %>%
  ggplot(., aes(x=condition, y=sig_count, fill=strain)) +
  geom_col(position='fill') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab('Relative Number of Conditionally Significant Genes') + 
  xlab('Condition') +
  guides(fill = guide_legend(title = 'Strain'))
ggsave('figures/ko_growth/relative_condition_sensitivity.pdf', p_con_sens_bar, width = 7, height = 5)

#### Proportional Overlap of conditionally significant genes between conditions in each strain ####
get_overlap <- function(con1, con2, strain){
  con1_genes <- filter(ko_growth_spread, condition == con1) %>% select(name, !!strain) %>% drop_na() %>% pull(name) %>% unique()
  con2_genes <- filter(ko_growth_spread, condition == con2) %>% select(name, !!strain) %>% drop_na() %>% pull(name) %>% unique()
  return(length(intersect(con1_genes, con2_genes))/((length(con1_genes) + length(con2_genes))/2))
}

get_overlap_pairs <- function(strain){
  return(mapply(get_overlap, MoreArgs = list(strain = strain), condition_combs[1,], condition_combs[2,]))
}

con_gene_overlaps <- sapply(strains, function(strain){tibble(condition1=condition_combs[1,],
                                                             condition2=condition_combs[2,],
                                                             overlap=get_overlap_pairs(strain))}, simplify = FALSE) %>%
  bind_rows(. ,.id='strain') %>%
  spread(key = strain, value = overlap) %>%
  unite(col = 'pair', remove = FALSE, condition1, condition2) %>%
  filter(!condition1 == condition2) %>%
  mutate(max = apply(select(., S288C:YPS), 1, max),
         min = apply(select(., S288C:YPS), 1, min),
         mean = apply(select(., S288C:YPS), 1, mean),
         diff = max - mean)

plotly_overlaps <- function(strain1, strain2){
  return(ggplotly(ggplot(con_gene_overlaps, aes_string(x=strain1, y=strain2, text='pair')) +
                    geom_point(size=0.5) + 
                    geom_abline(intercept = 0, slope = 1, colour='firebrick2') +
                    geom_abline(intercept = -0.2, slope = 1, colour='firebrick2', linetype=2) +
                    geom_abline(intercept = 0.2, slope = 1, colour='firebrick2', linetype=2)))
}

subplot(mapply(plotly_overlaps, strain_combs[1,], strain_combs[2,], SIMPLIFY = FALSE), nrows = 2, titleX = TRUE, titleY = TRUE)


######### By Gene ########
#### Gene-Gene Correlations ####
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

#### Per gene between strain correlation ####
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

#### Genes not required in one strain ####
odd_genes <- group_by(ko_growth, condition) %>%
  group_by(name, add = TRUE) %>%
  filter(qvalue < 0.01) %>%
  summarise(sig = n(),
            pos = sum(score > 0),
            neg = sum(score < 0),
            gene = first(gene),
            S288C_sig = 'S288C' %in% strain,
            UWOP_sig = 'UWOP' %in% strain,
            YPS_sig = 'YPS' %in% strain,
            Y55_sig = 'Y55' %in% strain) %>%
  filter(sig == 3, neg > 2) %>%
  left_join(., ko_growth_spread_all) %>%
  drop_na()

num_outliers <- dim(odd_genes)[1] - colSums(select(ungroup(odd_genes), S288C_sig:Y55_sig))
## Seems unlikely to have occured by chance
binom_tests <- sapply(num_outliers, binom.test, n=dim(odd_genes)[1], p=0.25)
multinom_test <- multinomial.test(num_outliers, rep(0.25, 4))


#### Genes with sig opposite effects ####
ko_growth_spread_diff <- filter(ko_growth_spread, apply(sign(select(ko_growth_spread,S288C:YPS)), 1,
                                                        function(x){!(min(x,na.rm = TRUE) == max(x, na.rm = TRUE))}))
write_tsv(ko_growth_spread_diff, 'data/ko_gene_opposite_growth.tsv')


#### Heatmaps of Genes of Interest ####
## Generic function to give heatmaps for genes in a set of cons
plot_con_gene_heatmaps <- function(tbl, genes, cons=NULL, strains=NULL, primary_strain='S288C', facet_cols=2, facet_rows=2, sig_level=0.01){
  if (is.null(cons)){
    cons <- unique(tbl$condition)
  }
  if (is.null(strains)){
    strains <- unique(tbl$strain)
  }
  if (is.null(primary_strain) | !primary_strain %in% strains){
    primary_strain <- strains[1]
  }
  
  
  tbl <- filter(tbl, name %in% genes, condition %in% cons)
  sig <- mutate(tbl, sig = ifelse(qvalue < sig_level, '*', '')) %>%
    select(strain, name, condition, sig)
  
  tbl <- select(tbl, strain, condition, name, score) %>%
    spread(key = condition, value = score) %>%
    gather(key = 'condition', value = 'score', -strain, -name) %>%
    left_join(., sig, by = c('strain', 'name', 'condition')) %>%
    mutate(sig = ifelse(is.na(sig), '', sig))
  
  limits <- c(min(tbl$score, na.rm = TRUE), max(tbl$score, na.rm = TRUE))
  
  if (limits[2] - limits[1] < 10){
    clrs <- c('firebrick2', 'white', 'cornflowerblue')
    vals <- scales::rescale(c(limits[1], 0, limits[2]))
  } else {
    clrs <- c('firebrick2', 'darkgoldenrod1', 'white', 'cornflowerblue', 'darkorchid')
    vals <- scales::rescale(c(limits[1], limits[1]/2, 0, limits[2]/2, limits[2]))
  }
  
  #### Plot heatmap faceted by strain
  ## Determine ordering based on primary strain
  mat <- split_strains(str=primary_strain, tbl=tbl, row='name', col='condition', var='score') %>%
    tbl_to_matrix(., row = 'name') %>%
    set_na(., 0)
  
  gene_dend <- as.dendrogram(hclust(dist(mat)))
  con_dend <- as.dendrogram(hclust(dist(t(mat))))
  gene_order <- rownames(mat)[order.dendrogram(gene_dend)]
  con_order <- colnames(mat)[order.dendrogram(con_dend)]
  
  # Add any gene/cons not in the primary strain in a random order
  gene_order <- c(gene_order, genes[!genes %in% gene_order])
  con_order <- c(con_order, cons[!cons %in% con_order])
  
  tbl %<>% mutate(name = factor(name, levels = gene_order, ordered = TRUE),
                  condition = factor(condition, levels = con_order, ordered = TRUE))
  
  p_prim_gene_dend <- ggdendrogram(gene_dend, rotate=TRUE)
  p_prim_con_dend <- ggdendrogram(con_dend, rotate=TRUE)
  
  ## plot heatmap
  p_strain_heatmaps <- ggplot(tbl, aes(x=condition, y=name, fill=score)) +
    geom_raster() +
    geom_text(aes(label=sig)) +
    facet_wrap(~strain, ncol=facet_cols, nrow=facet_rows) +
    xlab('Condition') +
    ylab('Gene') +
    scale_fill_gradientn(colours = clrs, limits=limits, na.value = 'black', values = vals) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    guides(fill=guide_colourbar(title = 'S-Score'))
  
  #### Plot general combined heatmap
  ## Determine ordering
  mat <- tbl %>% select(-sig) %>%
    spread(key = condition, value = score) %>%
    unite(col='name', strain, name) %>%
    tbl_to_matrix(., row = 'name') %>%
    set_na(., 0)
  
  gene_dend <- as.dendrogram(hclust(dist(mat)))
  con_dend <- as.dendrogram(hclust(dist(t(mat))))
  gene_order <- rownames(mat)[order.dendrogram(gene_dend)]
  con_order <- colnames(mat)[order.dendrogram(con_dend)]
  
  # Add any gene/cons not in the primary strain in a random order
  gene_order <- c(gene_order, genes[!genes %in% gene_order])
  con_order <- c(con_order, cons[!cons %in% con_order])
  
  tbl %<>% unite(col='unit', strain, name, remove = FALSE) %>%
    mutate(unit = factor(unit, levels = gene_order, ordered = TRUE),
           condition = factor(condition, levels = con_order, ordered = TRUE)) 
  
  p_all_gene_dend <- ggdendrogram(gene_dend, rotate=TRUE, labels = TRUE)
  p_all_con_dend <- ggdendrogram(con_dend, labels = TRUE)
  
  ## Plot heatmap
  p_all_heatmap <- ggplot(tbl, aes(x=condition, y=unit, fill=score)) +
    geom_raster() +
    geom_text(aes(label=sig)) +
    xlab('Condition') +
    ylab('Gene') +
    scale_fill_gradientn(colours = clrs, na.value = 'black', values = vals) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    guides(fill=guide_colourbar(title = 'S-Score'))
  
  return(list(strain_heatmap=p_strain_heatmaps,
              strain_gene_dend=p_prim_gene_dend,
              strain_condition_dend=p_prim_con_dend,
              all_heatmap=p_all_heatmap,
              all_gene_dend=p_all_gene_dend,
              all_condition_dend=p_all_con_dend))
}

## NaCl
nacl_cons <- grep('NaCl 0\\..M \\(', conditions, value = TRUE)
nacl_genes <- sapply(nacl_cons, get_sig_genes, growth_tbl = ko_growth, threshold=0.01) %>%
  unlist() %>%
  unique() %>%
  setdiff(., c('OPI9', 'ARV1', 'YOR345C')) # Excluded because not tested in UWOP

nacl_plots <- plot_con_gene_heatmaps(ko_growth, nacl_genes)
ggsave('figures/ko_growth/nacl_genes_strain_heatmaps.pdf', plot = nacl_plots$strain_heatmap, width = 13, height = 13)
ggsave('figures/ko_growth/nacl_genes_full_heatmaps.pdf', plot = nacl_plots$all_heatmap, width = 13, height = 15)
nacl_plots_sig <- plot_con_gene_heatmaps(ko_growth_sig, nacl_genes)
ggsave('figures/ko_growth/nacl_genes_strain_heatmaps_no_noise.pdf', plot = nacl_plots_sig$strain_heatmap, width = 13, height = 13)
ggsave('figures/ko_growth/nacl_genes_full_heatmaps_no_noise.pdf', plot = nacl_plots_sig$all_heatmap, width = 13, height = 15)

## Maltose/Glycerol
malt_cons <- c('Maltose 2% (48H)', 'Maltose 2% (72H)', 'Glycerol 2% (48H)', 'Glycerol 2% (72H)')
malt_genes <- sapply(malt_cons, get_sig_genes, growth_tbl = ko_growth, threshold=0.001) %>%
  unlist() %>%
  unique()

malt_plots <- plot_con_gene_heatmaps(ko_growth, malt_genes, primary_strain = 'UWOP')
ggsave('figures/ko_growth/malt_genes_strain_heatmaps.pdf', plot = malt_plots$strain_heatmap, width = 15, height = 17)
ggsave('figures/ko_growth/malt_genes_full_heatmaps.pdf', plot = malt_plots$all_heatmap, width = 15, height = 30)
malt_plots_sig <- plot_con_gene_heatmaps(ko_growth_sig, malt_genes, primary_strain = 'UWOP')
ggsave('figures/ko_growth/malt_genes_strain_heatmaps_no_noise.pdf', plot = malt_plots_sig$strain_heatmap, width = 15, height = 17)
ggsave('figures/ko_growth/malt_genes_full_heatmaps_no_noise.pdf', plot = malt_plots_sig$all_heatmap, width = 15, height = 30)

#### KO growth distances ####
ko_dists <- ko_growth_spread_all %>%
  mutate(S288C_UWOP = abs(S288C - UWOP),
         S288C_YPS = abs(S288C - YPS),
         S288C_Y55 = abs(S288C - Y55),
         Y55_UWOP = abs(Y55 - UWOP),
         YPS_UWOP = abs(YPS - UWOP),
         YPS_Y55 = abs(YPS - Y55))

p_dists <- ggplot(gather(ko_dists, key = 'pair', value = 'score_dist', S288C_UWOP:YPS_Y55), aes(x=pair, y=score_dist, colour=pair)) +
  geom_boxplot() +
  facet_wrap(~condition) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

p_dists_scatters <- ggpairs(ko_dists, columns = 3:8)

ko_dists_sum <- gather(ko_dists, key = 'pair', value = 'score_dist', S288C_UWOP:YPS_Y55) %>%
  group_by(condition, name) %>%
  summarise(mean_dist=mean(score_dist, na.rm = TRUE),
            max_dist=max(score_dist, na.rm = TRUE),
            min_dist=min(score_dist, na.rm = TRUE),
            num_large=sum(score_dist > 3, na.rm = TRUE),
            var_dist=var(score_dist, na.rm = TRUE),
            range=max_dist - min_dist,
            na_count=sum(is.na(score_dist)))

top_diff_genes <- filter(ko_dists_sum, na_count < 4, num_large > 2) %>%
  top_n(20, wt = mean_dist) %>%
  tbl_var_to_list(., 'name')

p_heat <- plot_con_gene_heatmaps(ko_growth, top_diff_genes$`39ºC (72H)`)
ggsave('figures/ko_growth/diff_heat_genes_strain_heatmap.pdf', p_heat$strain_heatmap, width = 12, height = 14)

p_caff <- plot_con_gene_heatmaps(ko_growth, top_diff_genes$`Caffeine 20mM (48H)`)


#### Gene set analysis ####
## Complexes
p_kornberg_complex <- plot_con_gene_heatmaps(ko_growth, unlist(filter(complexes, Complex == 'Kornberg\'s mediator (SRB) complex')$gene))
p_rpd3l_complex <- plot_con_gene_heatmaps(ko_growth, unlist(filter(complexes, Complex == 'Rpd3L complex')$gene))
p_ub_ligase_eradl_complex <- plot_con_gene_heatmaps(ko_growth, unlist(filter(complexes, Complex == 'ubiquitin ligase ERAD-L complex')$gene))

## Caffiene related sets based on tests performed
p_caff_dichol <- plot_con_gene_heatmaps(ko_growth, sets[['bp.dolichol_linked_oligosaccharide_biosynthetic_process(6)']], primary_strain = 'UWOP')
p_caff_swr1_comp <- plot_con_gene_heatmaps(ko_growth, unlist(filter(complexes, Complex == 'Swr1p complex')$gene))

## Gene Sets
osmotic_genes <- gene_sets_filt$bp[grep('osmo', names(gene_sets_filt$bp))] %>% unlist() %>% unique()
osmo_plots <- plot_con_gene_heatmaps(ko_growth, osmotic_genes, primary_strain = 'UWOP')
ggsave('figures/ko_growth/osmo_genes_strain_heatmaps.pdf', plot = osmo_plots$strain_heatmap, width = 20, height = 23)
ggsave('figures/ko_growth/osmo_genes_full_heatmaps.pdf', plot = osmo_plots$all_heatmap, width = 20, height = 40)
osmo_plots_sig <- plot_con_gene_heatmaps(ko_growth_sig, osmotic_genes, primary_strain = 'UWOP')
ggsave('figures/ko_growth/osmo_genes_strain_heatmaps_no_noise.pdf', plot = osmo_plots_sig$strain_heatmap, width = 20, height = 23)
ggsave('figures/ko_growth/osmo_genes_full_heatmaps_no_noise.pdf', plot = osmo_plots_sig$all_heatmap, width = 20, height = 40)

heat_genes <- gene_sets_filt$bp[grep('(temper|heat)', names(gene_sets_filt$bp))] %>% unlist %>% unique
heat_plots <- plot_con_gene_heatmaps(ko_growth, heat_genes)
ggsave('figures/ko_growth/heat_genes_strain_heatmaps.pdf', plot = heat_plots$strain_heatmap, width = 15, height = 17)
ggsave('figures/ko_growth/heat_genes_full_heatmaps.pdf', plot = heat_plots$all_heatmap, width = 15, height = 30)
heat_plots_sig <- plot_con_gene_heatmaps(ko_growth_sig, heat_genes)
ggsave('figures/ko_growth/heat_genes_strain_heatmaps_no_noise.pdf', plot = heat_plots_sig$strain_heatmap, width = 15, height = 17)
ggsave('figures/ko_growth/heat_genes_full_heatmaps_no_noise.pdf', plot = heat_plots_sig$all_heatmap, width = 15, height = 30)

aa_genes <- c(grep('amino_acid_biosynthetic',names(gene_sets_filt$bp),value = TRUE),
              "methionine_biosynthetic_process(7)",
              "glutamate_biosynthetic_process(8)",
              "arginine_biosynthetic_process(8)")
aa_genes <- gene_sets_filt$bp[aa_genes] %>% unlist() %>% unname() %>% unique()
aa_plots <- plot_con_gene_heatmaps(ko_growth, aa_genes)
ggsave('figures/ko_growth/aa_genes_strain_heatmaps.pdf', plot = aa_plots$strain_heatmap, width = 15, height = 23)


## Set Ks Tests - different to background in condition or strain-con pair
gene_set_test <- function(set, con, tbl){
  tbl %<>% filter(condition == con)
  if (sum(tbl$name %in% set) > 0){
    return(ks.test(filter(tbl, name %in% set)$score, filter(tbl, !name %in% set)$score)$p.value)
  }
  return(NA)
}

# set_ks_tests <- expand.grid(1:length(gene_sets_bp$geneset.names), conditions) %>%
#   as_tibble() %>%
#   rename(gene_set_id = Var1, condition = Var2) %>%
#   mutate(gene_set_name = gene_sets_bp$geneset.names[gene_set_id],
#          p.val = mapply(function(x,y){gene_set_test(gene_sets_bp$genesets[[x]],y,ko_growth)}, gene_set_id, condition),
#          p.adj = p.adjust(p.val,method = 'fdr'))
# saveRDS(set_ks_tests, 'data/Rdata/gene_set_condition_ks_tests.RDS')
set_ks_tests <- readRDS('data/Rdata/gene_set_condition_ks_tests.RDS')

set_ks_tests_mat <- select(set_ks_tests, -gene_set_name, -p.val) %>%
  mutate(p.adj = -log10(p.adj)) %>%
  spread(key = condition, value = p.adj) %>%
  tbl_to_matrix(., row = 'gene_set_id') %>%
  set_na(., 0) %>%
  set_inf(., 16)
heatmap.2(set_ks_tests_mat, col = colorRampPalette(colors = c('white','red'))(100), trace = 'none', margins = c(13,3))

## KS tests over all strain/condition pairs with all gene sets
ks_batches %<>% mutate(p.adj = p.adjust(p.value, method = 'fdr'))

p_pval_strain_density <- ggplot(ks_batches, aes(x=p.adj, colour=gene_set_group)) + geom_density() + facet_wrap(~strain)
p_pval_con_dist <- ggplot(ks_batches, aes(y=p.adj, x=condition)) + geom_boxplot() + theme(axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90))

p <- plot_con_gene_heatmaps(ko_growth, genes=gene_sets_filt$cc[c('ribosomal_subunit(4)')] %>% unlist() %>% unique())

# Complex T-tests - complexes different to background in strain-con pair
comp_t_test <- function(tbl, comp){
  ind <- tbl$name %in% comp
  # Return NAs if not enough genes in the set have been tested
  if (sum(ind) < 5){
    return(data_frame(estimate=NA, estimate1=NA, estimate2=NA,
                      statistic=NA, p.value=NA, parameter=NA,
                      conf.low=NA, conf.high=NA, method='Welch Two Sample t-test',
                      alternative='two.sided'))
  } else {
    return(tidy(t.test(tbl[ind,]$score, tbl[!ind,]$score)))
  }
}

do_comp_t_tests <- function(tbl, comps){
  return(
    bind_rows(
      sapply(comps, comp_t_test, tbl = tbl, simplify = FALSE),
      .id = 'complex'
    )
  )
}

complex_t_tests <- group_by(ko_growth, condition, strain) %>%
  do(do_comp_t_tests(., structure(complexes$gene, names=complexes$Complex)))%>%
  mutate(p.adj=p.adjust(p.value, method = 'fdr'))

## Compare each strain to the rest in each con
comp_strain_t_test <- function(tbl, str){
  ind <- tbl$strain == str
  # Return NAs if not enough genes in the set have been tested
  if (sum(ind) < 3 | sum(!ind) < 3){
    return(data_frame(estimate=NA, estimate1=NA, estimate2=NA,
                      statistic=NA, p.value=NA, parameter=NA,
                      conf.low=NA, conf.high=NA, method='Welch Two Sample t-test',
                      alternative='two.sided'))
  } else {
    return(tidy(t.test(tbl[ind,]$score, tbl[!ind,]$score)))
  }
}

do_comp_strain_t_tests <- function(tbl, comps, strains){
  return(
    bind_rows(
      sapply(strains, function(str){
        bind_rows(
          sapply(comps,
                 function(comp){comp_strain_t_test(tbl=filter(tbl, name %in% comp), str=str)},
                 simplify = FALSE),
          .id = 'complex'
          )
        },
        simplify = FALSE
      ),
      .id = 'strain'
    )
  )
}

complex_strain_t_tests <- group_by(ko_growth, condition) %>%
  do(do_comp_strain_t_tests(., structure(complexes$gene, names=complexes$Complex), strains)) %>%
  mutate(p.adj=p.adjust(p.value, method = 'fdr'))

set_strain_t_tests <- group_by(ko_growth, condition) %>%
  do(do_comp_strain_t_tests(., unlist(gene_sets_filt, recursive = FALSE), strains)) %>%
  mutate(p.adj=p.adjust(p.value, method = 'fdr'))

## Strain difference Anderson Darling tests
strain_diff_ad_tests <- function(tbl){
  str_sizes <- table(tbl$strain)
  strs <- names(str_sizes)[str_sizes > 4]
  # Return NAs if not enough genes in the set have been tested
  if (length(strs) < 2){
    return(data_frame(strains=0, S288C=str_sizes['S288C'], UWOP=str_sizes['UWOP'],
                      Y55=str_sizes['Y55'], YPS=str_sizes['YPS'], sample_total = sum(str_sizes),
                      ties=NA, sig=NA, ad=NA, t.ad=NA, asymp.p.val=NA, method='Anderson-Darling'))
  } else {
    t <- ad.test(lapply(strs, function(str){pull(filter(tbl, strain == str), score)}), method = 'asymptotic')
    return(data_frame(strains=length(strs), S288C=str_sizes['S288C'], UWOP=str_sizes['UWOP'],
                      Y55=str_sizes['Y55'], YPS=str_sizes['YPS'], sample_total = sum(str_sizes),
                      ties=t$n.ties, sig=t$sig, ad=t$ad[1,1], t.ad=t$ad[1,2], asymp.p.val=t$ad[1,3], method='Anderson-Darling'))
  }
}

do_strain_diff_ad_tests <- function(tbl, sets){
  return(
    bind_rows(
      sapply(sets,
             function(set){strain_diff_ad_tests(tbl=filter(tbl, name %in% set))},
             simplify = FALSE),
      .id = 'gene_set')
  )
}

complex_strain_diff_ad_tests <- group_by(ko_growth, condition) %>%
  do(do_strain_diff_ad_tests(., structure(complexes$gene, names=complexes$Complex))) %>%
  mutate(p.adj = p.adjust(asymp.p.val, method = 'fdr'))

growth_diff_cons <- c('Paraquat (48H)', 'Caffeine 20mM (48H)', 'Caffeine 15mM (48H)', 'Acetic acid (48H)', 'Maltose 2% (48H)', 'Maltose 2% (72H)')

# gene_set_strain_diff_ad_tests <- group_by(ko_growth, condition) %>%
#   do(do_strain_diff_ad_tests(., sets)) %>%
#   mutate(p.adj = p.adjust(asymp.p.val, method = 'fdr'))
# saveRDS(gene_set_strain_diff_ad_tests, 'data/Rdata/gene_set_ad_tests.RDS')
gene_set_strain_diff_ad_tests <- readRDS('data/Rdata/gene_set_ad_tests.RDS') %>%
  left_join(., set_meta)
View(filter(gene_set_strain_diff_ad_tests,
            p.adj < 0.05,
            condition %in% growth_diff_cons,
            set_size < 20))

p <- plot_con_gene_heatmaps(ko_growth, sets[['mf.iron_ion_transmembrane_transporter_activity(10)']], primary_strain = 'UWOP')
p <- plot_con_gene_heatmaps(ko_growth, sets[['bp.cytogamy(4)']], primary_strain = 'UWOP')

pathway_strain_diff_ad_tests <- group_by(ko_growth, condition) %>%
  do(do_strain_diff_ad_tests(., pathways)) %>%
  mutate(p.adj = p.adjust(asymp.p.val, method = 'fdr')) %>%
  left_join(., pathway_meta)

View(filter(pathway_strain_diff_ad_tests,
            p.adj < 0.05,
            set_size < 30))

p <- plot_con_gene_heatmaps(ko_growth,
                            pathways[['Biosynthesis of the N-glycan precursor (dolichol lipid-linked oligosaccharide, LLO) and transfer to a nascent protein']])

pathway_strain_diff_ad_tests_filtered_cons <- filter(ko_growth, condition %in% growth_diff_cons) %>%
  do(do_strain_diff_ad_tests(., pathways)) %>%
  mutate(p.adj = p.adjust(asymp.p.val, method = 'fdr')) %>%
  left_join(., pathway_meta)

## Use gene switching probabilities
gene_sig_summary <- filter(ko_growth, condition %in% growth_diff_cons) %>%
  mutate(signed_qvalue = sign(score) * qvalue) %>%
  select(strain, condition, name, signed_qvalue) %>%
  spread(key = 'strain', value = 'signed_qvalue')

switch_prob_summary_set <- filter(switch_probs, name %in% set_genes, condition %in% growth_diff_cons) %>%
  unite(col = 'pair', strain1, strain2) %>%
  select(condition, gene, name, pair, qval) %>%
  spread(key = pair, value = qval) %>%
  left_join(., gene_sig_summary, by=c("condition", "name")) %>%
  left_join(., bind_rows(sapply(sets, function(x){data_frame(name=x)}, simplify=FALSE), .id='gene_set'), by='name') %>%
  group_by(condition, gene_set) %>%
  summarise(S288C_switches = sum(Y55_S288C < 0.001, na.rm = TRUE) + sum(YPS_S288C < 0.001, na.rm = TRUE) + sum(UWOP_S288C < 0.001, na.rm = TRUE),
            UWOP_switches = sum(UWOP_Y55 < 0.001, na.rm = TRUE) + sum(YPS_UWOP < 0.001, na.rm = TRUE) + sum(UWOP_S288C < 0.001, na.rm = TRUE),
            Y55_switches = sum(Y55_S288C < 0.001, na.rm = TRUE) + sum(YPS_Y55 < 0.001, na.rm = TRUE) + sum(UWOP_Y55 < 0.001, na.rm = TRUE),
            YPS_switches = sum(YPS_S288C < 0.001, na.rm = TRUE) + sum(YPS_UWOP < 0.001, na.rm = TRUE) + sum(YPS_Y55 < 0.001, na.rm = TRUE),
            switches_per_strain = (S288C_switches + UWOP_switches + YPS_switches + Y55_switches)/4,
            S288C_non_sigs = sum(S288C > 0.001 & (abs(UWOP) < 0.001 | abs(YPS) < 0.001 | abs(Y55) < 0.001), na.rm = TRUE)) %>%
  left_join(., set_meta, by='gene_set') %>%
  mutate_at(vars(contains('switches'), 'S288C_non_sigs'), .funs = funs(./set_size))

## Identified sets have a strongly recuring theme of ALG/DIE genes
p <- plot_con_gene_heatmaps(ko_growth, sets[['mf.glucosyltransferase_activity(6)']], primary_strain = 'UWOP')
p <- plot_con_gene_heatmaps(ko_growth, sets[['bp.dolichol_linked_oligosaccharide_biosynthetic_process(6)']], primary_strain = 'UWOP')

p_dichol <- plot_con_gene_heatmaps(ko_growth,
                                   unique(c(sets[['bp.dolichol_linked_oligosaccharide_biosynthetic_process(6)']],
                                            sets[['mf.glucosyltransferase_activity(6)']])),
                                   primary_strain = 'UWOP')
ggsave('figures/ko_growth/dichol_genes_strain_heatmaps.pdf', plot = p_dichol$strain_heatmap, width = 15, height = 17)


## Method based on filtering switches first
switch_prob_summary <- filter(switch_probs, name %in% set_genes) %>%
  mutate_at(vars(contains('phenotype')), as.logical) %>%
  filter(phenotype1 | phenotype2,
         !sign(scores1) == sign(scores2)) %>%
  left_join(., bind_rows(sapply(sets, function(x){data_frame(name=x)}, simplify=FALSE), .id='gene_set'), by='name') %>%
  group_by(condition, gene_set) %>%
  summarise(switches = length(unique(name)),
            S288C_switches = sum(strain1 == 'S288C' | strain2 == 'S288C')) %>%
  left_join(., set_meta, by = 'gene_set') %>%
  mutate_at(vars(contains('switches')), .funs = funs(./set_size))

# Again bp.dolichol_linked_oligosaccharide_biosynthetic_process(6)
## For cyclohexamide, when looking for paraquat (both contain same gene set)
p <- plot_con_gene_heatmaps(ko_growth, sets[['cc.Sin3_type_complex(5)']])
p <- plot_con_gene_heatmaps(ko_growth, sets[['bp.negative_regulation_of_chromatin_silencing_at_silent_mating_type_cassette(6)']])

switch_prob_path_summary <- filter(switch_probs, name %in% unique(unlist(pathways))) %>%
  mutate_at(vars(contains('phenotype')), as.logical) %>%
  filter(phenotype1 | phenotype2, !sign(scores1) == sign(scores2)) %>%
  left_join(., bind_rows(sapply(pathways, function(x){data_frame(name=x)}, simplify=FALSE), .id='gene_set'), by='name') %>%
  group_by(condition, gene_set) %>%
  summarise(switches = length(unique(name)), S288C_switches = sum(strain2 == 'S288C')) %>%
  left_join(., pathway_meta, by = 'gene_set') %>%
  mutate_at(vars(contains('switches')), .funs = funs(./set_size))

# Again dolichol related appears to be best match
p_dichol_path <- plot_con_gene_heatmaps(ko_growth, pathways[['Biosynthesis of the N-glycan precursor (dolichol lipid-linked oligosaccharide, LLO) and transfer to a nascent protein']])
ggsave('figures/ko_growth/dichol_path_genes_strain_heatmaps.pdf', plot = p_dichol_path$strain_heatmap, width = 15, height = 17)

## Compare distribution of score diffs in S288C pairs vs other pairs
strain_diffs_test <- function(tbl){
  ind <- tbl$strain2 == 'S288C'
  # Return NAs if not enough genes in the set have been tested
  if (sum(ind) < 3 | sum(!ind) < 3){
    return(data_frame(statistic=NA, p.value=NA, method='Two-sample Kolmogorov-Smirnov test', alternative='two-sided'))
  } else {
    return(tidy(ks.test(tbl[ind,]$sub, tbl[!ind,]$sub)))
  }
}

do_strain_diffs_tests <- function(tbl, sets){
  return(
    bind_rows(
      sapply(sets,
             function(set){strain_diffs_test(tbl=filter(tbl, name %in% set))},
             simplify = FALSE),
      .id = 'gene_set')
  )
}

gene_set_strain_diffs_tests_condition <- group_by(switch_probs, condition) %>%
  do(do_strain_diffs_tests(., sets)) %>%
  mutate(p.adj = p.adjust(p.value, method = 'fdr')) %>%
  left_join(., set_meta, by='gene_set')

# # Not particularly useful
# gene_set_strain_diffs_tests_condition_filtered <- filter(switch_probs, name %in% set_genes) %>%
#   mutate_at(vars(contains('phenotype')), as.logical) %>%
#   filter(phenotype1 | phenotype2, !sign(scores1) == sign(scores2)) %>%
#   group_by(condition) %>%
#   do(do_strain_diffs_tests(., sets)) %>%
#   mutate(p.adj = p.adjust(p.value, method = 'fdr'))  %>%
#   left_join(., set_meta, by='gene_set')

gene_set_strain_diffs_tests <- do(switch_probs, do_strain_diffs_tests(., sets)) %>%
  mutate(p.adj = p.adjust(p.value, method = 'fdr'))  %>%
  left_join(., set_meta, by='gene_set')

gene_set_strain_diffs_tests_filt <- filter(switch_probs, condition %in% growth_diff_cons) %>%
  do(do_strain_diffs_tests(., sets)) %>%
  mutate(p.adj = p.adjust(p.value, method = 'fdr')) %>%
  left_join(., set_meta, by='gene_set')

# Identifies many groups that relate to NaCl but apparantly not to conditions of interest
p <- plot_con_gene_heatmaps(ko_growth, sets[['bp.cellular_response_to_salt_stress(6)']], primary_strain = 'UWOP')
p <- plot_con_gene_heatmaps(ko_growth, sets[['bp.regulation_of_ion_transport(5)']], primary_strain = 'UWOP')

pathway_strain_diffs_tests_condition <- group_by(switch_probs, condition) %>%
  do(do_strain_diffs_tests(., pathways)) %>%
  mutate(p.adj = p.adjust(p.value, method = 'fdr')) %>%
  left_join(., pathway_meta, by='gene_set')

#### Compare Gene sets vs random set ####
# Compare known gene sets to random sets of genes in various forms
# Determine group sizes for random samples
set_sizes <- sapply(sets, length)
set_size_exp_fit <- fitdistr(set_sizes, 'exponential')

hist(set_sizes, freq = FALSE)
curve(dexp(x, rate = set_size_exp_fit$estimate), from = 0, to = 700, add = TRUE)

# functions
range_span <- function(x, na.rm=TRUE){
  r <- range(x)
  if (is.infinite(r[1])){
    return(NA)
  } else{
    return(range(x, na.rm=na.rm)[2] - range(x, na.rm=na.rm)[1])
  }
}

summarise_gene_set <- function(tbl){
  return(summarise(tbl,
                   var=var(score, na.rm = TRUE),
                   range=range_span(score), mean=mean(score, na.rm = TRUE),
                   median=median(score, na.rm = TRUE),
                   skew=skewness(score, na.rm = TRUE),
                   max=max(score, na.rm = TRUE),
                   min=min(score, na.rm = TRUE),
                   sig_mean=mean(score[qvalue < 0.01], na.rm = TRUE),
                   abs_mean=mean(abs(score), na.rm = TRUE)))
}

apply_gene_score_fun <- function(genes, tbl){
  tbl <- filter(tbl, name %in% genes)
  
  none <- summarise_gene_set(tbl)
  strain <- summarise_gene_set(group_by(tbl, strain))
  con <- summarise_gene_set(group_by(tbl, condition))
  strain_con <- summarise_gene_set(group_by(tbl, strain, condition))
  
  return(mutate(bind_rows(list(none=none, strain=strain, condition=con, strain_condition=strain_con), .id = 'subgroup'), set_size=length(genes)))
}

random_gene_sets <- replicate(1000, sample(genes, max(rexp(1, rate = set_size_exp_fit$estimate), 5)))
names(random_gene_sets) <- paste0('randGroup', 1:length(random_gene_sets))
set_vars <- bind_rows(
  list(
    random_genes = bind_rows(
                    sapply(random_gene_sets,
                      apply_gene_score_fun,
                      tbl=ko_growth,
                      simplify = FALSE),
                    .id = 'gene_set'),
    gene_sets = bind_rows(
                  sapply(
                    sets,
                    apply_gene_score_fun,
                    tbl=ko_growth,
                    simplify = FALSE),
                  .id = 'gene_set')
  ),
  .id = 'type'
) %>%
  mutate(subgroup = factor(subgroup, levels = c('none', 'strain', 'condition', 'strain_condition'), ordered = TRUE),
         type = factor(type, levels = c('random_genes', 'gene_sets'), ordered = TRUE))

p_gene_set_boxes <- sapply(c('var', 'range', 'mean', 'median', 'skew', 'min', 'max', 'sig_mean', 'abs_mean'),
                           function(x){ggplot(set_vars, aes_string(x='subgroup', y=x, colour='type')) +
                                        geom_boxplot() +
                                        theme(axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45))},
                           simplify = FALSE)
p_gene_set_boxes_arr <- ggarrange(plotlist = p_gene_set_boxes, common.legend = TRUE)
ggsave('figures/ko_growth/gene_set_stats_boxes.pdf', p_gene_set_boxes_arr, width=14, height = 14)

p_gene_set_scatters <- sapply(list(mean_med=c('mean','median'), min_max=c('min', 'max'), size_var=c('set_size', 'var'), sigmean_mean=c('mean','sig_mean')),
                              function(x){ggplot(set_vars, aes_string(x=x[1], y=x[2], colour='subgroup')) + 
                                            geom_point() +
                                            facet_wrap(facets = vars(type))},
                              simplify = FALSE)
p_gene_set_scatters_arr <- ggarrange(plotlist = p_gene_set_scatters, common.legend = TRUE)
ggsave('figures/ko_growth/gene_set_stats_scatters.jpg', p_gene_set_scatters_arr, width = 14, height = 14)

set_mat <- set_vars %>%
  filter(type == 'gene_sets', subgroup == 'condition', set_size > 10) %>%
  select(gene_set, abs_mean, condition) %>%
  spread(key = condition, value = abs_mean) %>%
  tbl_to_matrix(., row = 'gene_set')

pdf('figures/ko_growth/gene_set_abs_mean_sscore_heatmap.pdf', width = 20, height = 20)
heatmap.2(set_mat, trace = 'none', col = colorRampPalette(c('white','red'))(100), margins = c(13,5), labRow = '')
dev.off()

set_con_summary <- set_vars %>% 
  filter(type == 'gene_sets', subgroup=='condition') %>%
  group_by(gene_set) %>%
  summarise(ind = which.max(abs_mean),
            condition=condition[ind],
            max_abs_mean=abs_mean[ind],
            max_mean = mean[ind],
            max_var = var[ind],
            enrich_abs_mean=max_abs_mean/mean(abs_mean),
            enrich_mean=max_mean/mean(mean),
            enrich_var=max_var/mean(var))

set_str_con_summary <- set_vars %>% 
  filter(type == 'gene_sets', subgroup=='strain_condition') %>%
  group_by(gene_set, strain) %>%
  summarise(ind = which.max(abs_mean),
            condition=condition[ind],
            max_abs_mean=abs_mean[ind],
            max_mean = mean[ind],
            max_var = var[ind],
            enrich_abs_mean=max_abs_mean/mean(abs_mean),
            enrich_mean=max_mean/mean(mean),
            enrich_var=max_var/mean(var))

## Bigger differentiation between gene sets effect in most impactful strain than average in S288C
p_mean_enrich <- ggplot(set_str_con_summary, aes(x=strain, y=enrich_abs_mean)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90)) +
  geom_hline(yintercept = 1)
p_mean_enrich_per_con <- ggplot(set_str_con_summary, aes(x=condition, y=enrich_abs_mean, colour=strain)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90)) +
  geom_hline(yintercept = 1)
ggsave('figures/ko_growth/abs_mean_enrichment_of_most_impacted_con.pdf',
       ggarrange(p_mean_enrich, p_mean_enrich_per_con,
                 nrow = 1, ncol = 2, widths = c(1,4),align = 'hv'),
       width = 20, height = 10)

## Test number of sig genes in real vs random sets
gene_set_lengths <- sapply(sets, length)
set_strain_con_num_sig <- data_frame(gene_set = rep(names(sets), gene_set_lengths),
                                   name = unname(unlist(sets)),
                                   set_size = rep(gene_set_lengths, gene_set_lengths)) %>%
  left_join(., ko_growth, by='name') %>%
  mutate(abs_score = abs(score)) %>%
  group_by(strain, condition, gene_set) %>%
  summarise(set_size=first(set_size),
            n_sig = sum(qvalue < 0.01, na.rm = TRUE),
            prop_sig = n_sig/set_size)


## Determine expected number of significant genes for a set of size n in given strain/condition
rand_sets <- read_tsv('data/rand_gene_set_expected_sig_counts.tsv', col_names = TRUE)

rand_strain_set_mean_props <- rand_sets %>%
  mutate(mean_prop = mean_sig/set_size) %>%
  group_by(strain, condition) %>%
  summarise(mean = mean(mean_prop), var=var(mean_prop))

N <- length(genes)
strain_set_props <- ko_growth %>%
  group_by(strain, condition) %>%
  summarise(expected_prop = sum(qvalue<0.01)/N)

func_set_sizes <- data_frame(func_size=sapply(sets, function(x){sum(x %in% genes)}), gene_set = names(sets))

# Binomial test?
set_strain_con_num_sig %<>% left_join(., strain_set_props, by = c('strain', 'condition')) %>%
  left_join(., func_set_sizes, by = 'gene_set') %>%
  mutate(expected_n_sig = func_size * expected_prop,
         prop_enriched = prop_sig/expected_prop)

# Dichols are most enriched in caffeine as found in other analyses
p <- plot_con_gene_heatmaps(ko_growth, gene_sets_bp_filt[['arginine_biosynthetic_process(8)']])

