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

p_gene_strain_cor_boxes <- ggplot(gather(gene_ko_profile_cor, key = 'strain_pair', value = 'cor', UWOP_S288C:YPS_Y55), aes(x=strain_pair, y=cor)) +
  geom_boxplot()
ggsave('figures/ko_growth/gene_cors_between_strains.pdf', p_gene_strain_cor_boxes, width = 8, height = 6)

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

p_gene_sig_vs_cor_boxes <- ggplot(gene_ko_profile_cor, aes(x=sig_level, y=mean_cor)) +
  geom_boxplot() +
  xlab('Number of strains with q.value < 0.01') +
  ylab('Mean Correlation Between Strains') +
  stat_summary(geom = 'text', fun.data = function(x){return(c(y = -0.55, label = length(x)))})
ggsave('figures/ko_growth/mean_gene_cor_by_sig.pdf', p_gene_sig_vs_cor_boxes, width = 8, height = 5)

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


