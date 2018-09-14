# Script performing analysis of gene sets (pathways/sets/complexes) impact during knockouts
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
#### 1 Manually looking for interesting complexes/sets ####
## Gene Sets
man_osmotic_genes <- gene_sets_filt$bp[grep('osmo', names(gene_sets_filt$bp))] %>% unlist() %>% unique()

man_heat_genes <- gene_sets_filt$bp[grep('(temper|heat)', names(gene_sets_filt$bp))] %>% unlist %>% unique

aa_genes <- c(grep('amino_acid_biosynthetic',names(gene_sets_filt$bp),value = TRUE),
              "methionine_biosynthetic_process(7)",
              "glutamate_biosynthetic_process(8)",
              "arginine_biosynthetic_process(8)")
aa_genes <- gene_sets_filt$bp[aa_genes] %>% unlist() %>% unname() %>% unique()

mal_genes <- grep('MAL', genes, value = TRUE)
########

#### 2 Sets determined by conditional impact ####
## NaCl
nacl_cons <- grep('NaCl 0\\..M \\(', conditions, value = TRUE)
nacl_genes <- sapply(nacl_cons, get_sig_genes, growth_tbl = ko_growth, threshold=0.01) %>%
  unlist() %>%
  unique() %>%
  setdiff(., c('OPI9', 'ARV1', 'YOR345C')) # Excluded because not tested in UWOP

# No gene set contains more than a few of the NaCl sig genes
nacl_sets <- data_frame(set = names(sets), count = sapply(sets, function(x){sum(x %in% nacl_genes)}), len = gene_set_lengths) %>%
  mutate(rel_count = count/len)

## Maltose/Glycerol
malt_cons <- c('Maltose 2% (48H)', 'Maltose 2% (72H)', 'Glycerol 2% (48H)', 'Glycerol 2% (72H)')
malt_genes <- sapply(malt_cons, get_sig_genes, growth_tbl = ko_growth, threshold=0.001) %>%
  unlist() %>%
  unique()

########

#### 3 Sets based on KO growth distances ####
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
########

#### 4 KS tests on set S-Score within each strain/con background ####
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
#heatmap.2(set_ks_tests_mat, col = colorRampPalette(colors = c('white','red'))(100), trace = 'none', margins = c(13,3))

## KS tests over all strain/condition pairs with all gene sets
ks_batches %<>% mutate(p.adj = p.adjust(p.value, method = 'fdr')) %>%
  mutate(gene_set = paste(gene_set_group, gene_set, sep = '.')) %>%
  select(-gene_set_group) %>%
  left_join(., set_meta, by = "gene_set")

p_pval_strain_density <- ggplot(ks_batches, aes(x=p.adj, colour=gene_set_group)) + geom_density() + facet_wrap(~strain)
p_pval_con_dist <- ggplot(ks_batches, aes(y=p.adj, x=condition)) + geom_boxplot() + theme(axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90))
########

#### 5 T-tests on difference to background in each strain/con pair ####
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
########

#### 6 T-tests comparing each strain to the others ####
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
########

#### 7 Anderson-Darling tests (difference between strains in a gene set in a condition) ####
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

## Complexes
complex_strain_diff_ad_tests <- group_by(ko_growth, condition) %>%
  do(do_strain_diff_ad_tests(., structure(complexes$gene, names=complexes$Complex))) %>%
  mutate(p.adj = p.adjust(asymp.p.val, method = 'fdr'))

## Gene Sets
# gene_set_strain_diff_ad_tests <- group_by(ko_growth, condition) %>%
#   do(do_strain_diff_ad_tests(., sets)) %>%
#   mutate(p.adj = p.adjust(asymp.p.val, method = 'fdr'))
# saveRDS(gene_set_strain_diff_ad_tests, 'data/Rdata/gene_set_ad_tests.RDS')
gene_set_strain_diff_ad_tests <- readRDS('data/Rdata/gene_set_ad_tests.RDS') %>%
  left_join(., set_meta)

## Pathways
pathway_strain_diff_ad_tests <- group_by(ko_growth, condition) %>%
  do(do_strain_diff_ad_tests(., pathways)) %>%
  mutate(p.adj = p.adjust(asymp.p.val, method = 'fdr')) %>%
  left_join(., pathway_meta)

pathway_strain_diff_ad_tests_filtered_cons <- filter(ko_growth, condition %in% growth_diff_cons) %>%
  do(do_strain_diff_ad_tests(., pathways)) %>%
  mutate(p.adj = p.adjust(asymp.p.val, method = 'fdr')) %>%
  left_join(., pathway_meta)
########

#### 8 Gene Switching Probabilities ####
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

## Method based on filtering switches first - anomolies seen in relation to stars on heatmaps from different sig threshold
switch_prob_summary <- filter(switch_probs, name %in% set_genes) %>%
  mutate_at(vars(contains('phenotype')), as.logical) %>%
  filter(xor(phenotype1, phenotype2) | (!(sign(scores1) == sign(scores2)) & qval< 0.01)) %>%
  left_join(., bind_rows(sapply(sets, function(x){data_frame(name=x)}, simplify=FALSE), .id='gene_set'), by='name') %>%
  group_by(condition, gene_set) %>%
  summarise(switches = length(unique(name)),
            S288C_switches = sum(strain1 == 'S288C' | strain2 == 'S288C')) %>%
  left_join(., set_meta, by = 'gene_set') %>%
  mutate_at(vars(contains('switches')), .funs = funs(./set_size))

# Again bp.dolichol_linked_oligosaccharide_biosynthetic_process(6)

# # Method restricted to S288C switches
# switch_prob_summary <- filter(switch_probs, name %in% set_genes, condition %in% c("Cyclohexamide (48H)", "39ºC (48H)", "39ºC (72H)")) %>%
#   mutate_at(vars(contains('phenotype')), as.logical) %>%
#   filter(xor(phenotype1, phenotype2) | (!(sign(scores1) == sign(scores2)) & qval< 0.01),
#          strain2 == 'S288C' & phenotype2) %>%
#   left_join(., bind_rows(sapply(sets, function(x){data_frame(name=x)}, simplify=FALSE), .id='gene_set'), by='name') %>%
#   group_by(condition, gene_set) %>%
#   summarise(switches = length(unique(name))) %>%
#   left_join(., set_meta, by = 'gene_set') %>%
#   mutate(switches = switches/set_size)

switch_prob_path_summary <- filter(switch_probs, name %in% unique(unlist(pathways))) %>%
  mutate_at(vars(contains('phenotype')), as.logical) %>%
  filter(phenotype1 | phenotype2, !sign(scores1) == sign(scores2)) %>%
  left_join(., bind_rows(sapply(pathways, function(x){data_frame(name=x)}, simplify=FALSE), .id='gene_set'), by='name') %>%
  group_by(condition, gene_set) %>%
  summarise(switches = length(unique(name)), S288C_switches = sum(strain2 == 'S288C')) %>%
  left_join(., pathway_meta, by = 'gene_set') %>%
  mutate_at(vars(contains('switches')), .funs = funs(./set_size))

# Again dolichol related appears to be best match
########

#### 9 KS tests on gene set gene's strain differences in each condition ####
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

pathway_strain_diffs_tests_condition <- group_by(switch_probs, condition) %>%
  do(do_strain_diffs_tests(., pathways)) %>%
  mutate(p.adj = p.adjust(p.value, method = 'fdr')) %>%
  left_join(., pathway_meta, by='gene_set')
########

## Reorganise plots into lists?
#### 10 Compare Gene sets vs random set ####
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
########

#### Interesting Sets Identified (from previous sections and manual) ####
# Identified by the section number they were drawn from
## Heat
p_manual_heat <- plot_con_gene_heatmaps(ko_growth, man_heat_genes) #1
ggsave('figures/ko_growth/heat_genes_strain_heatmaps.pdf', plot = p_manual_heat$strain_heatmap, width = 15, height = 17)

p_kornberg_complex <- plot_con_gene_heatmaps(ko_growth, unlist(filter(complexes, Complex == 'Kornberg\'s mediator (SRB) complex')$gene)) #1

p_heat <- plot_con_gene_heatmaps(ko_growth, top_diff_genes$`39ºC (72H)`) #3
ggsave('figures/ko_growth/diff_heat_genes_strain_heatmap.pdf', p_heat$strain_heatmap, width = 12, height = 14)

p_tubulin_heat <- plot_con_gene_heatmaps(ko_growth, sets[['bp.tubulin_complex_assembly(6)']]) #8

## Caffiene
p_caff_diff <- plot_con_gene_heatmaps(ko_growth, top_diff_genes$`Caffeine 20mM (48H)`) #3
p_caff_swr1_comp <- plot_con_gene_heatmaps(ko_growth, unlist(filter(complexes, Complex == 'Swr1p complex')$gene)) #10

p_dichol <- plot_con_gene_heatmaps(ko_growth,
                                   unique(c(sets[['bp.dolichol_linked_oligosaccharide_biosynthetic_process(6)']],
                                            sets[['mf.glucosyltransferase_activity(6)']])),
                                   primary_strain = 'UWOP') # 10, 8, 7
ggsave('figures/ko_growth/dichol_genes_strain_heatmaps.pdf', plot = p_dichol$strain_heatmap, width = 15, height = 17)

p_dichol_path <- plot_con_gene_heatmaps(ko_growth, pathways[['Biosynthesis of the N-glycan precursor (dolichol lipid-linked oligosaccharide, LLO) and transfer to a nascent protein']]) #8, 7
ggsave('figures/ko_growth/dichol_path_genes_strain_heatmaps.pdf', plot = p_dichol_path$strain_heatmap, width = 15, height = 17)

## NaCl
p_manual_osmo <- plot_con_gene_heatmaps(ko_growth, man_osmotic_genes, primary_strain = 'UWOP') #1
ggsave('figures/ko_growth/osmo_genes_strain_heatmaps.pdf', plot = p_manual_osmo$strain_heatmap, width = 20, height = 23)

p_resp_to_salt <- plot_con_gene_heatmaps(ko_growth, sets[['bp.cellular_response_to_salt_stress(6)']], primary_strain = 'UWOP') #9
p_ion_transport_reg <- plot_con_gene_heatmaps(ko_growth, sets[['bp.regulation_of_ion_transport(5)']], primary_strain = 'UWOP') #9
nacl_plots <- plot_con_gene_heatmaps(ko_growth, nacl_genes) #2
ggsave('figures/ko_growth/nacl_genes_strain_heatmaps.pdf', plot = nacl_plots$strain_heatmap, width = 13, height = 13)

## Cyclohexamide
p_rpd3l_complex <- plot_con_gene_heatmaps(ko_growth, unlist(filter(complexes, Complex == 'Rpd3L complex')$gene)) #1
# For cyclohexamide, when looking for paraquat (both contain same gene set)
p_sin3_comp <- plot_con_gene_heatmaps(ko_growth, sets[['cc.Sin3_type_complex(5)']]) #8
p_neg_chrom_silence <- plot_con_gene_heatmaps(ko_growth, sets[['bp.negative_regulation_of_chromatin_silencing_at_silent_mating_type_cassette(6)']]) #8
p_chrom_silence <- plot_con_gene_heatmaps(ko_growth, sets[['bp.regulation_of_chromatin_silencing_at_rDNA(6)']]) #8

# Also heat and caff
p_his_ubi <- plot_con_gene_heatmaps(ko_growth, sets[['bp.histone_ubiquitination(5)']]) #8

## Maltose/Glycerol
p_mal_genes <- plot_con_gene_heatmaps(ko_growth, mal_genes) #1, manually identified mal gene example
p_rib_subunit <- plot_con_gene_heatmaps(ko_growth, genes=gene_sets_filt$cc[c('ribosomal_subunit(4)')] %>% unlist() %>% unique()) #4
malt_plots <- plot_con_gene_heatmaps(ko_growth, malt_genes, primary_strain = 'UWOP') #2
ggsave('figures/ko_growth/malt_genes_strain_heatmaps.pdf', plot = malt_plots$strain_heatmap, width = 15, height = 17)

##AA Starvation
p_aa <- plot_con_gene_heatmaps(ko_growth, aa_genes) #1
ggsave('figures/ko_growth/aa_genes_strain_heatmaps.pdf', plot = p_aa$strain_heatmap, width = 15, height = 23)
########

#### Test Genes of interest for mutations ####
# Identify variants carried by the ko strains
uwop_muts <- geno$mut_id[geno$SACE_GAT > 0]
y55_muts <- geno$mut_id[geno$ADA > 0]

## Look at whole dichol gene set (was significant in caffiene)
imp_gens <- filter(imp, name %in% unique(c(sets[['bp.dolichol_linked_oligosaccharide_biosynthetic_process(6)']],
                                           sets[['mf.glucosyltransferase_activity(6)']])))
imp_strs <- filter(imp_gens, mut_id %in% c(uwop_muts, y55_muts))

p <- plot_ly(filter(imp_gens, type=='nonsynonymous'), x=~freq, y=~sift_score, text=~id)

## Look only at genes with phenotypes
imp_gens_restricted <- filter(imp, name %in% c('TPS1', 'ALG3', 'DIE2', 'ALG5', 'ALG8', 'ALG6'))
imp_strs_restricted <- filter(imp_gens_restricted, mut_id %in% c(uwop_muts, y55_muts))
# Single non-synonomous variant found in dichol genes with phenotypes (chrVII:947885_G/A), only carried in Y55

# Look at frequency vs sift score in nonsynonymous variants in phenotype genes - single high freq, high sift variant - YGR227W M61V
p <- plot_ly(filter(imp_gens_restricted, type=='nonsynonymous'), x=~freq, y=~sift_score, text=~id)

## Look only at genes with phenotypes
imp_strs_both <- filter(imp_gens, mut_id %in% intersect(uwop_muts, y55_muts))

# Look at frequency vs sift score in nonsynonymous variants in phenotype genes - single high freq, high sift variant - YGR227W M61V
p <- plot_ly(filter(imp_gens_both, type=='nonsynonymous'), x=~freq, y=~sift_score, text=~id)
########

#### Test Gene Sets for RNA seq enrichment ####
rna_seq_sets <- data_frame(gene_set = rep(names(sets), gene_set_lengths), name = unlist(sets)) %>%
  left_join(., rna_seq, by='name') %>%
  group_by(gene_set, strain) %>%
  summarise(all = sum(padj < 0.01),
            up = sum(padj < 0.01 & log2FoldChange > 0),
            down = sum(padj < 0.01 & log2FoldChange < 0)) %>%
  drop_na() %>%
  gather(key = 'type', value = 'count', all, up, down) %>%
  unite(col = 'unite', strain, type, sep = '_') %>%
  spread(key = 'unite', value = 'count') %>%
  mutate(total_all = UWOPS87_2421_all + Y55_all + YPS606_all,
         total_up = UWOPS87_2421_up + Y55_all + YPS606_up,
         total_down = UWOPS87_2421_down + Y55_down + YPS606_down) %>%
  left_join(., data_frame(gene_set=names(gene_set_lengths), set_size=gene_set_lengths), by='gene_set') %>%
  mutate_at(.vars = vars(-gene_set, set_size), funs(norm = ./set_size)) %>%
  mutate_at(.vars = vars(matches('total_.+_norm')), funs(./3))
########

#### RNA Seq Comparison in Interesting Sets ####
## Heat
p_manual_heat_expr <- plot_rna_seq_fold_change(rna_seq, gene_names = man_heat_genes) #1

p_kornberg_complex_expr <- plot_rna_seq_fold_change(rna_seq, gene_names = unlist(filter(complexes, Complex == 'Kornberg\'s mediator (SRB) complex')$gene)) #1

p_heat_expr <- plot_rna_seq_fold_change(rna_seq, gene_names = top_diff_genes$`39ºC (72H)`) #3

p_tubulin_heat_expr <- plot_rna_seq_fold_change(rna_seq, gene_names = sets[['bp.tubulin_complex_assembly(6)']]) #8

## Caffiene
p_caff_diff_expr <- plot_rna_seq_fold_change(rna_seq, gene_names = top_diff_genes$`Caffeine 20mM (48H)`) #3
p_caff_swr1_comp_expr <- plot_rna_seq_fold_change(rna_seq, gene_names = unlist(filter(complexes, Complex == 'Swr1p complex')$gene)) #10

p_dichol_expr <- plot_rna_seq_fold_change(rna_seq, gene_names = unique(c(sets[['bp.dolichol_linked_oligosaccharide_biosynthetic_process(6)']],
                                                                         sets[['mf.glucosyltransferase_activity(6)']]))) # 10, 8, 7
ggsave('figures/ko_growth/dichol_genes_expression.pdf', plot = p_dichol_expr, width = 8, height = 8)
p_dichol_expr_caff <- plot_rna_seq_fold_change(rna_seq_caff, gene_names = unique(c(sets[['bp.dolichol_linked_oligosaccharide_biosynthetic_process(6)']],
                                                                         sets[['mf.glucosyltransferase_activity(6)']]))) # 10, 8, 7
ggsave('figures/ko_growth/dichol_genes_expression_caff.pdf', plot = p_dichol_expr_caff, width = 8, height = 8)

p_dichol_path_expr <- plot_rna_seq_fold_change(rna_seq, gene_names =  pathways[['Biosynthesis of the N-glycan precursor (dolichol lipid-linked oligosaccharide, LLO) and transfer to a nascent protein']]) #8, 7

## NaCl
p_manual_osmo_expr <- plot_rna_seq_fold_change(rna_seq, gene_names =  man_osmotic_genes) #1

p_resp_to_salt_expr <- plot_rna_seq_fold_change(rna_seq, gene_names = sets[['bp.cellular_response_to_salt_stress(6)']]) #9
p_ion_transport_reg_expr <- plot_rna_seq_fold_change(rna_seq, gene_names = sets[['bp.regulation_of_ion_transport(5)']]) #9
p_nacl_expr <- plot_rna_seq_fold_change(rna_seq, gene_names = nacl_genes) #2

## Cyclohexamide
p_rpd3l_complex_expr <- plot_rna_seq_fold_change(rna_seq, gene_names = unlist(filter(complexes, Complex == 'Rpd3L complex')$gene)) #1
# For cyclohexamide, when looking for paraquat (both contain same gene set)
p_sin3_comp_expr <- plot_rna_seq_fold_change(rna_seq, gene_names = sets[['cc.Sin3_type_complex(5)']]) #8
p_neg_chrom_silence_expr <- plot_rna_seq_fold_change(rna_seq, gene_names = sets[['bp.negative_regulation_of_chromatin_silencing_at_silent_mating_type_cassette(6)']]) #8
p_chrom_silence_expr <- plot_rna_seq_fold_change(rna_seq, gene_names = sets[['bp.regulation_of_chromatin_silencing_at_rDNA(6)']]) #8

# Also heat and caff
p_his_ubi_expr <- plot_rna_seq_fold_change(rna_seq, gene_names = sets[['bp.histone_ubiquitination(5)']]) #8

## Maltose/Glycerol
p_mal_genes_expr <- plot_rna_seq_fold_change(rna_seq, gene_names = mal_genes) #1, manually identified mal gene example
p_rib_subunit_expr <- plot_rna_seq_fold_change(rna_seq, gene_names = gene_sets_filt$cc[c('ribosomal_subunit(4)')] %>% unlist() %>% unique()) #4
p_malt_expr <- plot_rna_seq_fold_change(rna_seq, gene_names = malt_genes) #2

##AA Starvation
p_aa_expr <- plot_rna_seq_fold_change(rna_seq, gene_names = aa_genes) #1
########
