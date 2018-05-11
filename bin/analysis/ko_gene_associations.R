# Script performing analyses on the liti growth data look for associations between genes 
# indicated as important in deletion studies
setwd('~/Projects/hog/')
library(tidyverse)
library(magrittr)
library(ggpubr)

#### Import Data ####
strains <- readRDS('data/Rdata/strain_meta.rds')
filtered_strains <- filter(strains, Ploidy == 2, Aneuploidies == 'euploid') %>% pull(`Standardized name`)
# Filter strains that are unusually distant from the reference genome
filtered_strains <-  setdiff(filtered_strains , c("AMH", "BAG", "BAH", "BAL", "CEG", "CEI"))
distance_to_ref <- structure(strains$`Total number of SNPs`, names=strains$`Standardized name`)

genes <- readRDS('data/Rdata/gene_meta_all.rds')

complexes <- readRDS('data/Rdata/complex_members.rds')

essential <- readRDS('data/Rdata/essential_genes.rds')
essential_genes <- filter(essential, essential == 'E') %>% pull(locus)
non_essential_genes <- filter(essential, essential == 'NE') %>% pull(locus)

impacts <- readRDS('data/Rdata/all_muts_impacts.rds')

growth <- readRDS('data/Rdata/growth_liti.rds') %>%
  filter(strain %in% filtered_strains)

allele_freqs <- readRDS('data/Rdata/allele_freqs.rds')

probs <- readRDS('data/Rdata/paff_all_genes.rds') %>%
  filter(strain %in% intersect(filtered_strains, growth$strain)) %>%
  left_join(. ,growth, by = 'strain') %>%
  gather(key='condition', value = 'growth', -strain, -gene, -p_aff) %>%
  mutate(gene_sig = FALSE)

# import Ko growth data showing which genes are significant (using 0.01 threshold)
ko_growth <- readRDS('data/Rdata/gene_ko_growth.rds') %>%
  filter(num_strains >= 2) %>% # Only interested in genes that show general ipact rather than in one strain
  filter(gene %in% unique(probs$gene))

ko_growth_full <- readRDS('data/Rdata/gene_ko_growth_full.rds') %>%
  filter(gene %in% unique(probs$gene))

#### Processing ####
ko_thresh <- 0.95

## Define conditions considered equivalent between bede and liti screens
equiv_conditions <- structure(c("NaCl 0.6M (48H)", "39ºC (48H)", "Caffeine 20mM (48H)", "Glycerol 2%  (48H)", "6-AU (48H)"),
                              names=c("ypdnacl15m", "ypd40", "ypdcafein40", "ypglycerol", "ypd6au"))

sig_genes <- lapply(equiv_conditions, function(x){filter(ko_growth, condition == x, num_strains >= 1) %>% pull(gene)})

for (con in names(sig_genes)){
  probs[probs$condition == con & probs$gene %in% sig_genes[[con]], 'gene_sig'] <- TRUE
}

#### Analysis ####
## By strain over all sig genes
strain_summary <- group_by(probs, condition) %>%
  group_by(strain, add=TRUE) %>%
  summarise(num_sig_kos = sum(gene_sig & p_aff > ko_thresh),
            growth = first(growth))

p_ypdnacl15m_sig_kos <- ggplot(filter(strain_summary, condition == 'ypdnacl15m'), aes(x=num_sig_kos, y=growth)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  ylab('Growth Relative to YPD Media') +
  xlab(paste0('Number of NaCl Stress Significant Genes with P(Aff) > ', ko_thresh))
fit <- lm(growth ~ num_sig_kos, data = filter(strain_summary, condition == 'ypdnacl15m'))
ggsave('figures/liti_growth/nacl15_growth_vs_number_sig_kos.pdf', p_ypdnacl15m_sig_kos, width = 12, height = 10)

## By individual gene
probs %<>% mutate(ko = p_aff > ko_thresh) %>%
  mutate(sig_ko = gene_sig & ko)

# Example with Pbs2 and Hog1 (Known to be sig in NaCl 1.5mM)
hog_id <- 'YLR113W'
pbs2_id <- 'YJL128C'
probs_nacl <- filter(probs, gene==!!hog_id | gene == !!pbs2_id, condition=='ypdnacl15m') %>%
  select(strain, gene, growth, ko) %>%
  spread(key = 'gene', value = 'ko') %>%
  rename(hog1 = !!hog_id, pbs2 = !!pbs2_id) %>%
  mutate(ko = ifelse(hog1, ifelse(pbs2, 'Both', 'Hog1'), ifelse(pbs2, 'Pbs2', 'Neither'))) %>%
  select(-hog1, -pbs2) %>%
  mutate(ypdnacl15m = growth) %>%
  mutate(ypd14 = !!growth$ypd14) %>%
  gather(key = 'condition', value = 'growth', -strain, -growth, -ko) %>%
  mutate(ko = factor(ko, levels = c('Neither', 'Hog1', 'Pbs2', 'Both'))) %>%
  mutate(condition = factor(condition, levels = c('ypdnacl15m', 'ypd14')))

p_osmotic_shock_ko_growth_box <- ggplot(probs_nacl, aes(x = ko, y = growth, colour = ko)) +
  facet_wrap(~ condition) +
  geom_boxplot(varwidth = TRUE) +
  xlab(paste0('Genes with P(aff) > ', ko_thresh)) +
  ylab('Growth Relative to YPD Media') + 
  stat_compare_means(method = 'wilcox.test',
                     comparisons = list(c('Hog1', 'Neither'), c('Pbs2', 'Neither'))) +
  stat_summary(geom = 'text', fun.data = function(x){return(c(y = -0.03, label = length(x)))}) +
  guides(colour = FALSE)
ggsave('figures/liti_growth/osmotic_shock_ko_growth_pbs2_hog1.pdf', p_osmotic_shock_ko_growth_box, width = 14, height = 10)

# Generic function for gene ko growth box plots
plot_ko_growth_box <- function(condition, gene_names=NULL, gene_ids=NULL){
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
  
  probs_filter <- filter(probs, gene %in% gene_ids, condition==!!condition) %>%
    select(strain, gene, growth, ko) %>%
    spread(key = 'gene', value = 'ko') %>%
    rename_at(.vars = vars(!!gene_ids), .funs = function(x){gene_names[x]}) %>%
    gather(key = 'gene', value = 'ko', -strain, -growth) %>%
    mutate(ko = factor(ifelse(ko, 'True', 'False'), levels = c('False', 'True'))) %>%
    mutate(gene = factor(gene, levels = gene_names))
  
  p <- ggplot(probs_filter, aes(x=ko, y=growth)) +
    geom_boxplot(aes(fill=gene), varwidth = TRUE, alpha=0.4) +
    facet_wrap(~gene) +
    stat_compare_means(comparisons = list(c('False','True'))) +
    stat_summary(geom = 'text', fun.data = function(x){return(c(y = -0.03, label = length(x)))}) +
    xlab(paste0('P(aff) > ', ko_thresh)) +
    ylab(paste0('Growth in ', condition, ' Relative to YPD Media')) +
    guides(fill=FALSE)
  return(p)
}

p_nacl15m_all_sig_kos_box <- plot_ko_growth_box('ypdnacl15m', gene_ids = sig_genes$ypdnacl15m)
ggsave('figures/liti_growth/osmotic_shock_ko_growth_all_genes.pdf', p_nacl15m_all_sig_kos_box, width = 20, height = 20)

p_heat_all_sig_kos_box <- plot_ko_growth_box('ypd40', gene_ids = sig_genes$ypd40)
ggsave('figures/liti_growth/high_temp_ko_growth_all_genes.pdf', p_heat_all_sig_kos_box, width = 30, height = 30)

p_cafein40_all_sig_kos_box <- plot_ko_growth_box('ypdcafein40', gene_ids = sig_genes$ypdcafein40)
ggsave('figures/liti_growth/cafein40_ko_growth_all_genes.pdf', p_cafein40_all_sig_kos_box, width = 30, height = 30)

p_glycerol_all_sig_kos_box <- plot_ko_growth_box('ypglycerol', gene_ids = sig_genes$ypglycerol)
ggsave('figures/liti_growth/glycerol_ko_growth_all_genes.pdf', p_glycerol_all_sig_kos_box, width = 30, height = 30)

p_6au_all_sig_kos_box <- plot_ko_growth_box('ypd6au', gene_ids = sig_genes$ypd6au)
ggsave('figures/liti_growth/6au_ko_growth_all_genes.pdf', p_6au_all_sig_kos_box, width = 30, height = 30)

## Analyse all genes
all_sig_genes <- unique(unname(unlist(sig_genes)))
gene_summary <- readRDS('data/Rdata/paff_all_genes.rds') %>%
  mutate(sig = gene %in% all_sig_genes) %>%
  mutate(sig_nacl = gene %in% sig_genes$ypdnacl15m) %>%
  mutate(ko = p_aff > ko_thresh) %>%
  group_by(gene) %>%
  summarise(count = sum(ko),
            sig_nacl = first(sig_nacl),
            sig=first(sig)) %>%
  mutate(essential = structure(essential$essential, names=essential$locus)[gene])

# Essential and sig KO genes are mutually exclusive as expected (S-Score would mean no diff between stressed and unstressed growth in essential gene)
table(gene_summary$sig_nacl, gene_summary$essential)

# Not a big difference in rate of KO between genes that are sig in some conditions
gene_summary %<>% mutate(sigess = if_else(essential == 'E', 'Essential', if_else(sig, 'Significant', 'Neither')))
p_ko_sig_rates <- ggplot(filter(gene_summary, !is.na(essential)), aes(x=count, colour=sigess)) +
  geom_density()

## In fact appear to be ko'd significantly more, albeit by a small amount
# t.test(count ~ sig, data = gene_summary)
# data:  count by sig
# t = -2.3647, df = 765.9, p-value = 0.01829
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -7.303760 -0.677893
# sample estimates:
#   mean in group FALSE  mean in group TRUE 
# 29.44828            33.43910 

## But significance is down to essential genes
# t.test(count ~ sig, data = filter(gene_summary, essential == 'NE'))
# data:  count by sig
# t = -1.0507, df = 755.77, p-value = 0.2937
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -5.537457  1.676464
# sample estimates:
#   mean in group FALSE  mean in group TRUE 
# 32.72335            34.65385

### Analyse genes ko rate against their significance
ko <- filter(ko_growth_full, condition == "NaCl 0.6M (48H)")

get_strain_sscore <- function(gene, strain, ko){
  if (gene %in% ko$gene){
    return(unlist(ko[ko$gene == gene, paste0(strain, '-score')]))
  } else {
    return(NA)
  }
}

gene_summary %<>% mutate(S288C = sapply(gene, get_strain_sscore, strain = 'S288C', ko = ko)) %>%
  mutate(UWOP = sapply(gene, get_strain_sscore, strain = 'UWOP', ko = ko)) %>%
  mutate(Y55 = sapply(gene, get_strain_sscore, strain = 'Y55', ko = ko)) %>%
  mutate(YPS = sapply(gene, get_strain_sscore, strain = 'YPS', ko = ko)) %>%
  gather(key = 'strain', value = 'sscore', -gene, -count, -sig, -essential, -sigess, -sig_nacl)

p_ko_count_significance <- ggplot(gene_summary, aes(x=sscore, y=count, colour=strain, shape=sig_nacl)) +
  geom_point() +
  xlab('S-Score in 0.6mM NaCl (48hr)') +
  ylab(paste0('Number of Strains with P(Aff) > ', ko_thresh)) +
  guides(colour=guide_legend(title = 'KO Strain for S-Score'), shape=guide_legend(title = 'Gene Significant in 2+ strains'))

ggsave('figures/liti_growth/strain_ko_count_vs_significance_nacl.pdf', p_ko_count_significance, width = 12, height = 10)

p_strain_sscore_density <- ggplot(gene_summary, aes(x=sscore, colour=strain)) +
  geom_density()
