# Script performing analyses on bede's growth data across all strains
setwd('~/Projects/hog/')
library(rlang)
library(tidyverse)
library(magrittr)
library(ggpubr)
library(gplots)

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
  mutate(ko = p_aff > ko_thresh) %>%
  left_join(., select(growth, strain, condition, score, qvalue), by='strain')

path <- readRDS('data/Rdata/hog_path_probs.rds') %>%
  filter(strain %in% growth$strain) %>%
  left_join(., select(growth, strain, condition, score, qvalue), by='strain')

ko_growth <- readRDS('data/Rdata/gene_ko_growth.rds') %>%
  filter(gene %in% unique(probs$gene))

ko_growth_full <- readRDS('data/Rdata/gene_ko_growth_full.rds') %>%
  filter(gene %in% unique(probs$gene))


#### Analysis ####
sig_genes <- sapply(unique(growth$condition), function(x){filter(ko_growth, condition == x, num_strains >= 1) %>% pull(gene)})
sig_genes <- sig_genes[sapply(sig_genes, length) > 0]

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

## Use fisher method to test variants
fisher_method <- function(condition, prob_tbl=NULL, gene_ids=NULL, gene_names=NULL, thresh = ko_thresh){
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
  
  probs_filter <- filter(prob_tbl, gene %in% gene_ids, condition==!!condition, p_aff > thresh) %>%
    select(strain, gene, score, qvalue) %>%
    mutate(gene = gene_names[gene]) %>%
    group_by(gene) %>%
    summarise(k = n(),
              sum_pos = sum(score > 0),
              sum_neg = sum(score < 0),
              fisher_X = -2*sum(log(qvalue))) %>%
    mutate(fisher_p = pchisq(fisher_X, df = 2*k, lower.tail = FALSE)) %>%
    mutate(fisher_p_adj = p.adjust(p, method = 'fdr')) %>%
    mutate(condition = condition)
  
  return(probs_filter)
}

fisher_thresh <- ko_thresh
fisher_combined_probs <- bind_rows(lapply(names(sig_genes_strong), function(x){fisher_method(x, gene_ids = sig_genes_strong[[x]], prob_tbl = probs, thresh = fisher_thresh)}))
write_tsv(fisher_combined_probs, paste0('data/bede_condition_sig_gene_affect_growth_impacts_', fisher_thresh, '.tsv'))

## Using of FDR adjusted p-values before combining is a problem because it potentially violates the assumption of uniform distribution over [0,1]
## that means the value follows a chi squared distribution

## Path analysis
# prob Model
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

# All conditions vs hog
hog_prob_lm_test <- function(con, tbl){
  fit <- summary(lm(score ~ hog_probability, data = filter(tbl, condition == con)))
  return(c(p_val = pf(fit$fstatistic[1], fit$fstatistic[2], fit$fstatistic[3], lower.tail = FALSE), R_sqr = fit$r.squared, Adj_R_sqr=fit$adj.r.squared))
}

hog_sig_cons <- c("NaCl 0.4M (48H)", "NaCl 0.4M (72H)", "NaCl 0.6M (48H)", "NaCl 0.6M (72H)", "Sorbitol 1M (48H)")
hog_semi_sig_cons <- c("NaCl 0.4M + 39ºC (48H)", "NaCl 0.4M + 39ºC (72H)", "NaCl 0.6M + 39ºC (48H)", "NaCl 0.6M + 39ºC (72H)")
hog_prob_tests <- sapply(unique(path$condition), hog_prob_lm_test, tbl=path) %>%
  t(.) %>%
  as_tibble(rownames='condition') %>%
  rename(p_val = p_val.value) %>%
  mutate(hog_sig = if_else(condition %in% hog_sig_cons, 'Signif', if_else(condition %in% hog_semi_sig_cons, 'Partial', 'None'))) %>%
  mutate(log_p = -log10(p_val))

p_hog_prob_conditional_sig_p <- ggplot(hog_prob_tests, aes(x=condition, y=log_p, fill=hog_sig)) +
  geom_col() +
  geom_hline(yintercept = -log10(0.01), linetype='dashed') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  guides(fill=guide_legend(title = 'HOG Significance')) +
  xlab('Condition') +
  ylab(expr('-log'[10]~'p'))
  
ggsave(paste0(figure_root, 'hog_path_prob_condition_cor_p_vals.pdf'), p_hog_prob_conditional_sig_p, width = 8, height = 5)


# Number of KO'd OS sig genes
prob_nacl <- filter(probs,
                    condition %in% c("NaCl 0.4M (72H)", "NaCl 0.6M (72H)", "Sorbitol 1M (48H)"),
                    gene %in% sig_genes$`NaCl 0.6M (72H)`) %>%
  select(strain, gene, ko, condition, score, qvalue) %>%
  spread(key = gene, value = ko) %>%
  mutate(ko_sum = rowSums(select(., -strain, -condition, -score,-qvalue)))

fit4 <- lm(score ~ ko_sum, filter(prob_nacl, condition == "NaCl 0.4M (72H)"))
fit6 <- lm(score ~ ko_sum, filter(prob_nacl, condition == "NaCl 0.6M (72H)"))
fitS <- lm(score ~ ko_sum, filter(prob_nacl, condition == "Sorbitol 1M (48H)"))

p_nacl_sig_gene_count_growth_cor <- ggplot(prob_nacl, aes(x=ko_sum, y=score)) +
  geom_point() +
  geom_smooth(method='lm') +
  facet_wrap(~condition)
ggsave(paste0(figure_root, 'nacl_sig_ko_count_growth.pdf'), p_nacl_sig_gene_count_growth_cor, width = 7, height = 5)

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


## Intersection of significant genes
cols <- colorRampPalette(c('white','blue','red'))
pdf('figures/bede_growth_condition_strong_sig_gene_heatmap.pdf', width = 10, height = 10)
sig_gene_overlap_strong <- sapply(sig_genes_strong, function(x){sapply(sig_genes_strong, function(y){length(intersect(x,y))/length(x)})})
width <- dim(sig_gene_overlap_strong)[1]
den <- as.dendrogram(hclust(dist(sig_gene_overlap_strong)))
heatmap.2(t(sig_gene_overlap_strong), Rowv = den, Colv = den, revC = TRUE, trace = 'none', breaks = width + 1, margins = c(15,15),
          colsep = 1:width, rowsep = 1:width, sepwidth = c(0.025,0.025), sepcolor = 'black', col = cols(width))
dev.off()

pdf('figures/bede_growth_condition_sig_gene_heatmap.pdf', width = 10, height = 10)
sig_gene_overlap <- sapply(sig_genes, function(x){sapply(sig_genes, function(y){length(intersect(x,y))/length(x)})})
width <- dim(sig_gene_overlap)[1]
den <- as.dendrogram(hclust(dist(sig_gene_overlap)))
heatmap.2(t(sig_gene_overlap), Rowv = den, Colv = den, revC = TRUE, trace = 'none', breaks = width + 1, margins = c(15,15),
          colsep = 1:width, rowsep = 1:width, sepwidth = c(0.025,0.025), sepcolor = 'black', col = cols(width))
dev.off()
