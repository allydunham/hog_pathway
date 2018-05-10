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

ko_growth <- readRDS('data/Rdata/gene_ko_growth.rds')

#### Processing ####
ko_thresh <- 0.95

## Define conditions considered equivalent between bede and liti screens
equiv_conditions <- structure(c("NaCl 0.6M (48H)"),
                              names=c("ypdnacl15m"))

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
probs %<>% mutate(sig_ko = gene_sig & p_aff > ko_thresh)
hog_id <- 'YLR113W'
pbs2_id <- 'YJL128C'
p_hog_ko_growth_box <- ggplot(filter(probs, gene==!!hog_id , condition=='ypdnacl15m'), aes(x=sig_ko, y=growth)) +
  geom_boxplot(varwidth = TRUE) +
  xlab(paste0('Hog1 P(aff) > ', ko_thresh)) +
  ylab('Growth Relative to YPD Media') + 
  stat_compare_means(method = 't.test', hjust = -0.5) +
  stat_summary(geom = 'text', fun.data = function(x){return(c(y = 0, label = length(x)))})

p_pbs_ko_growth_box <- ggplot(filter(probs, gene==!!pbs2_id, condition=='ypdnacl15m'), aes(x=sig_ko, y=growth)) +
  geom_boxplot(varwidth = TRUE) +
  xlab(paste0('Pbs2 P(aff) > ', ko_thresh)) +
  ylab('') + 
  stat_compare_means(method = 't.test', hjust = -0.5) +
  stat_summary(geom = 'text', fun.data = function(x){return(c(y = 0, label = length(x)))})

p_osmotic_shock_ko_growth_box <- ggarrange(p_hog_ko_growth_box, p_pbs_ko_growth_box)
ggsave('figures/liti_growth/osmotic_shock_ko_growth.pdf', p_osmotic_shock_ko_growth_box, width = 14, height = 10)

all_sig_genes <- unique(unname(unlist(sig_genes)))
gene_summary <- readRDS('data/Rdata/paff_all_genes.rds') %>%
  mutate(sig = gene %in% all_sig_genes) %>%
  mutate(ko = p_aff > ko_thresh) %>%
  group_by(gene) %>%
  summarise(count = sum(ko),
            sig = first(sig)) %>%
  mutate(essential = structure(essential$essential, names=essential$locus)[gene])

# Essential and sig KO genes are mutually exclusive as expected (S-Score would mean no diff between stressed and unstressed growth in essential gene)
table(gene_summary$sig, gene_summary$essential)

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