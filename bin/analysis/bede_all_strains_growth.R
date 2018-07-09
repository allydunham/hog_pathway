# Script performing analyses on bede's growth data across all strains
setwd('~/Projects/hog/')
library(rlang)
library(tidyverse)
library(magrittr)
library(ggpubr)

figure_root <- 'figures/bede_growth_all/'

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

ko_thresh <- 0.95
probs <- readRDS('data/Rdata/paff_all_genes.rds') %>%
  filter(strain %in% growth$strain) %>%
  mutate(ko = p_aff > 0.95) %>%
  left_join(., select(growth, strain, condition, score, qvalue), by='strain')

path <- readRDS('data/Rdata/hog_path_probs.rds') %>%
  filter(strain %in% growth$strain)

ko_growth <- readRDS('data/Rdata/gene_ko_growth.rds') %>%
  filter(gene %in% unique(probs$gene))

ko_growth_full <- readRDS('data/Rdata/gene_ko_growth_full.rds') %>%
  filter(gene %in% unique(probs$gene))


#### Analysis ####
sig_genes_strong <- sapply(unique(growth$condition), function(x){filter(ko_growth, condition == x, num_strains >= 2) %>% pull(gene)})
sig_genes_strong <- sig_genes_strong[sapply(sig_genes_strong, length) > 0]

# Generic function for gene ko growth box plots - taken from ko_gene_associations.R
plot_ko_growth_box <- function(condition, gene_names=NULL, gene_ids=NULL, prob_tbl=probs, ylab=NULL, xlab=NULL, nrow=NULL, ncol=NULL,
                               ylim=NULL, xvar='ko', yvar='qvalue'){
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
    select(strain, gene, !!yvar, !!xvar) %>%
    spread(key = 'gene', value = xvar) %>%
    rename_at(.vars = vars(!!gene_ids), .funs = function(x){gene_names[x]}) %>%
    gather(key = 'gene', value = !!xvar, -strain, -!!yvar) %>%
    mutate(!!xvar := factor(ifelse(!!sym(xvar), 'True', 'False'), levels = c('False', 'True'))) %>%
    mutate(gene = factor(gene, levels = gene_names))
  
  if (is.null(ylab)){
    ylab = paste0('S-Score in ', condition)
  }
  if (is.null(xlab)){
    xlab = paste0('P(Aff) > ', ko_thresh)
  }
  if (is.null(ylim)){
    r <- range(probs_filter[,yvar], na.rm = TRUE)
    err <- (r[2] - r[1]) * 0.3
    ylim = c(r[1] - err, r[2] + err)
  }
  
  p <- ggplot(probs_filter, aes_string(x=xvar, y=yvar)) +
    geom_boxplot(aes(fill=gene), varwidth = TRUE, alpha=0.4) +
    facet_wrap(~gene, nrow = nrow, ncol = ncol) +
    stat_compare_means(comparisons = list(c('False','True'))) +
    stat_summary(geom = 'text', fun.data = function(x){return(c(y = -0.03, label = length(x)))}) +
    xlab(xlab) +
    ylab(ylab) +
    ylim(ylim[1], ylim[2]) +
    guides(fill=FALSE)
  return(p)
}

assoc_plots <- lapply(names(sig_genes_strong), function(x){plot_ko_growth_box(x, gene_ids = sig_genes_strong[[x]])})


