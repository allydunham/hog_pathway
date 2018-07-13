# Script performing analyses on bede's growth data across all strains
setwd('~/Projects/hog/')
library(rlang)
library(tidyverse)
library(magrittr)
library(ggpubr)

figure_root <- 'figures/bede_growth_all/'
ko_thresh <- 0.95

#### Import Data ####
strains <- readRDS('data/Rdata/strain_meta.rds')
distance_to_ref <- structure(strains$`Total number of SNPs`, names=strains$`Standardized name`)

filtered_strains <- filter(strains, Ploidy <= 2, Aneuploidies == 'euploid') %>%
  pull(`Standardized name`) %>%
  setdiff(. , c("AMH", "BAG", "BAH", "BAL", "CEG", "CEI"))

genes <- readRDS('data/Rdata/gene_meta_all.rds')

essential <- readRDS('data/Rdata/essential_genes.rds')
essential_genes <- filter(essential, essential == 'E') %>% pull(locus)
non_essential_genes <- filter(essential, essential == 'NE') %>% pull(locus)

growth <- read_tsv('data/raw/yeasts_liti_fixed.tsv', col_names = TRUE, col_types = cols(strain = col_character())) %>%
  filter(info %in% filtered_strains) %>% 
  filter(subset == 'liti') %>%
  rename(strain_id = strain, strain=info)


probs <- readRDS('data/Rdata/paff_all_genes.rds') %>%
  filter(strain %in% growth$strain) %>%
  mutate(ko = p_aff > 0.95) %>%
  left_join(., select(growth, strain, condition, score, qvalue), by='strain')

path <- readRDS('data/Rdata/hog_path_probs.rds') %>%
  filter(strain %in% growth$strain) %>%
  left_join(., select(growth, strain, condition, score, qvalue), by='strain')

ko_growth <- readRDS('data/Rdata/gene_ko_growth.rds') %>%
  filter(gene %in% unique(probs$gene))

ko_growth_full <- readRDS('data/Rdata/gene_ko_growth_full.rds') %>%
  filter(gene %in% unique(probs$gene))


#### Analysis ####
sig_genes_strong <- sapply(unique(growth$condition), function(x){filter(ko_growth, condition == x, num_strains >= 2) %>% pull(gene)})
sig_genes_strong <- sig_genes_strong[sapply(sig_genes_strong, length) > 0]

# Generic function for gene ko growth box plots - taken from ko_gene_associations.R
plot_ko_growth_box <- function(condition, gene_names=NULL, gene_ids=NULL, prob_tbl=probs, ylab=NULL, xlab=NULL, nrow=NULL, ncol=NULL,
                               ylim=NULL){
  if (!is.null(gene_ids)){
    gene_names <- structure(genes$name, names=genes$id)[gene_ids]
    gene_names[is.na(gene_names)] <- gene_ids[is.na(gene_names)]
    gene_ids <- structure(names(gene_names), names=unname(gene_names))
  } else if (!is.null(gene_names)){
    gene_ids <- structure(genes$id, names=genes$name)[str_to_upper(gene_names)]
    gene_names <- structure(names(gene_ids), names=unname(gene_ids))
  } else {
    warning('One of "gene_names" or "gene_ids" must be given')
  }
  
  probs_filter <- filter(prob_tbl, gene %in% gene_ids, condition==!!condition) %>%
    select(strain, gene, score, ko) %>%
    spread(key = 'gene', value = 'ko') %>%
    rename_at(.vars = vars(!!gene_ids), .funs = function(x){gene_names[x]}) %>%
    gather(key = 'gene', value = ko, -strain, -score) %>%
    mutate(ko = factor(ifelse(ko, 'True', 'False'), levels = c('False', 'True'))) %>%
    mutate(gene = factor(gene, levels = gene_names))
  
  if (is.null(ylab)){
    ylab = paste0('S-Score in ', condition)
  }
  if (is.null(xlab)){
    xlab = paste0('P(Aff) > ', ko_thresh)
  }
  if (is.null(ylim)){
    r <- range(probs_filter[,'score'], na.rm = TRUE)
    err <- (r[2] - r[1]) * 0.3
    ylim = c(r[1] - err, r[2] + err)
  }
  
  p <- ggplot(probs_filter, aes(x=ko, y=score)) +
    geom_boxplot(aes(fill=gene), varwidth = TRUE, alpha=0.4) +
    facet_wrap(~gene, nrow = nrow, ncol = ncol) +
    stat_compare_means(comparisons = list(c('False','True'))) +
    stat_summary(geom = 'text', fun.data = function(x){return(c(y = ylim[1], label = length(x)))}) +
    stat_summary(geom = 'text', aes(colour=gene), fun.data = function(x){return(c(y = ylim[1] + 0.5*err, label = signif(mean(x), 3)))}) +
    xlab(xlab) +
    ylab(ylab) +
    ylim(ylim[1], ylim[2]) +
    guides(fill=FALSE, colour=FALSE)
  return(p)
}

assoc_plots <- lapply(names(sig_genes_strong), function(x){plot_ko_growth_box(x, gene_ids = sig_genes_strong[[x]])})
cons <- gsub('^X', '', gsub('\\.+', '_', make.names(names(sig_genes_strong))))
for (i in 1:length(assoc_plots)){ggsave(paste0(figure_root, 'assocs/' , cons[i], 'ko_association_plot.pdf'), assoc_plots[[i]],
                                                width = 10, height = 10)}

## Modified function with no gene faceting
plot_ko_growth_box_no_facet <- function(condition, gene_names=NULL, gene_ids=NULL, prob_tbl=probs, ylab=NULL, xlab=NULL, nrow=NULL, ncol=NULL,
                               ylim=NULL){
  if (!is.null(gene_ids)){
    gene_names <- structure(genes$name, names=genes$id)[gene_ids]
    gene_names[is.na(gene_names)] <- gene_ids[is.na(gene_names)]
    gene_ids <- structure(names(gene_names), names=unname(gene_names))
  } else if (!is.null(gene_names)){
    gene_ids <- structure(genes$id, names=genes$name)[str_to_upper(gene_names)]
    gene_names <- structure(names(gene_ids), names=unname(gene_ids))
  } else {
    warning('One of "gene_names" or "gene_ids" must be given')
  }
  
  probs_filter <- filter(prob_tbl, gene %in% gene_ids, condition==!!condition) %>%
    select(strain, gene, score, ko) %>%
    spread(key = 'gene', value = 'ko') %>%
    mutate(ko = rowSums(select(., -strain, -score)) > 0) %>%
    mutate(ko_count = rowSums(select(., -strain, -score))) %>%
    select(strain, score, ko, ko_count) %>%
    mutate(ko = factor(ifelse(ko, 'True', 'False'), levels = c('False', 'True')))
  
  if (is.null(ylab)){
    ylab = paste0('S-Score in ', condition)
  }
  if (is.null(xlab)){
    xlab = paste0('P(Aff) > ', ko_thresh, ' in any conditionally significant gene')
  }
  if (is.null(ylim)){
    r <- range(probs_filter[,'score'], na.rm = TRUE)
    err <- (r[2] - r[1]) * 0.3
    ylim = c(r[1] - err, r[2] + err)
  }
  
  p_box <- ggplot(probs_filter, aes(x=ko, y=score, fill=ko)) +
    geom_boxplot(varwidth = TRUE, alpha=0.4) +
    stat_compare_means(comparisons = list(c('False','True'))) +
    stat_summary(geom = 'text', fun.data = function(x){return(c(y = ylim[1], label = length(x)))}) +
    stat_summary(geom = 'text', colour='firebrick2', fun.data = function(x){return(c(y = ylim[1] + 0.5*err, label = signif(mean(x), 3)))}) +
    xlab(xlab) +
    ylab(ylab) +
    ylim(ylim[1], ylim[2]) +
    guides(fill=FALSE, colour=FALSE) + 
    scale_fill_manual(values = c(True='cornflowerblue', False='firebrick2'))
  
  fit <- summary(lm(score~ko_count, data=probs_filter))
  p_count <- ggplot(probs_filter, aes(x=ko_count, y=score)) +
    geom_boxplot(aes(group=cut_width(ko_count, 1, center = 0)), fill='cornflowerblue') + 
    geom_smooth(method = 'lm') +
    ggtitle(paste0('Adj. R squared = ', signif(fit$adj.r.squared, 4), ', p = ', signif(pf(fit$fstatistic[1],fit$fstatistic[2],fit$fstatistic[3],lower.tail = FALSE), 3))) +
    ylab('') +
    xlab(paste0('Number of conditionally significant genes where P(Aff) > ', ko_thresh))
  return(ggarrange(p_box, p_count))
}
assoc_plots_no_facet <- lapply(names(sig_genes_strong), function(x){plot_ko_growth_box_no_facet(x, gene_ids = sig_genes_strong[[x]])})
for (i in 1:length(assoc_plots_no_facet)){ggsave(paste0(figure_root, 'assocs/' , cons[i], 'ko_association_plot_all_genes.pdf'), 
                                                 assoc_plots_no_facet[[i]],width = 10, height = 6)}


## Generic function to plot volcano style plots per gene
plot_ko_growth_volcano <- function(condition, gene_names=NULL, gene_ids=NULL, prob_tbl=probs, ylab=NULL, xlab=NULL, nrow=NULL, ncol=NULL){
  if (!is.null(gene_ids)){
    gene_names <- structure(genes$name, names=genes$id)[gene_ids]
    gene_names[is.na(gene_names)] <- gene_ids[is.na(gene_names)]
    gene_ids <- structure(names(gene_names), names=unname(gene_names))
  } else if (!is.null(gene_names)){
    gene_ids <- structure(genes$id, names=genes$name)[str_to_upper(gene_names)]
    gene_names <- structure(names(gene_ids), names=unname(gene_ids))
  } else {
    warning('One of "gene_names" or "gene_ids" must be given')
  }
  
  probs_filter <- filter(prob_tbl, gene %in% gene_ids, condition==!!condition) %>%
    select(strain, gene, ko, score, qvalue) %>%
    spread(key = 'gene', value = 'ko') %>%
    rename_at(.vars = vars(!!gene_ids), .funs = function(x){gene_names[x]}) %>%
    gather(key = 'gene', value = 'ko', -strain, -score, -qvalue) %>%
    mutate(ko := factor(ifelse(ko, 'True', 'False'), levels = c('False', 'True'))) %>%
    mutate(gene = factor(gene, levels = gene_names))
  
  if (is.null(ylab)){
    ylab = 'Q-Value'
  }
  if (is.null(xlab)){
    xlab = paste0('S-Score in ', condition)
  }
  
  p <- ggplot(probs_filter, aes(x=score, y=qvalue)) +
    geom_line() +
    geom_point(data = filter(probs_filter, ko == 'True'), aes(colour=ko)) +
    facet_wrap(~gene, nrow = nrow, ncol = ncol) +
    xlab(xlab) +
    ylab(ylab) +
    scale_color_manual(values = c(True='red',False='black')) +
    guides(colour=FALSE) +
    ggtitle(condition) +
    theme(plot.title = element_text(hjust = 0.5)) +
    geom_hline(yintercept = 0.05, linetype = 'dashed', colour = 'blue')
  return(p)
}
assoc_volcano_plots <- lapply(names(sig_genes_strong), function(x){plot_ko_growth_volcano(x, gene_ids = sig_genes_strong[[x]])})

for (i in 1:length(assoc_volcano_plots)){ggsave(paste0(figure_root, 'volcanos/' , cons[i], 'volcano_ko_plot.pdf'), assoc_volcano_plots[[i]],
                                                width = 10, height = 10)}

## Path analysis
fit4 <- lm(score ~ hog_probability, filter(path, condition == "NaCl 0.4M (72H)"))
fit6 <- lm(score ~ hog_probability, filter(path, condition == "NaCl 0.6M (72H)"))
fitS <- lm(score ~ hog_probability, filter(path, condition == "Sorbitol 1M (48H)"))

p_path_prob <- ggplot(filter(path, condition %in% c("NaCl 0.4M (72H)", "NaCl 0.6M (72H)", "Sorbitol 1M (48H)")),
                      aes(x=hog_probability, y=score)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  facet_wrap(~ condition) +
  xlab('P(Hog)') +
  ylab('S-Score')
ggsave(paste0(figure_root, 'path_probability.pdf'), p_path_prob, width = 7, height = 7)


## General linear model test
probs_ace <- filter(probs, condition == 'Anaerobic growth (48H)', gene %in% sig_genes_strong$`Anaerobic growth (48H)`) %>%
  select(-ko) %>%
  spread(key = 'gene', value = 'p_aff')

fit <- lm(score ~ .^2, data = select(probs_ace, -strain, -condition, -qvalue))

probs_ace %<>% mutate(pred_score = predict(fit, newdata = probs_ace)) %>%
  mutate(sig = qvalue < 0.05)

p_lm <- ggplot(probs_ace, aes(x=score, y=pred_score, colour=sig)) +
  geom_point() +
  ylim(-5, 5) +
  xlim(-5,5) +
  geom_abline(slope = 1, intercept = 0)
