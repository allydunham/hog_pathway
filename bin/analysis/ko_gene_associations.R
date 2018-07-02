# Script performing analyses on the liti growth data look for associations between genes 
# indicated as important in deletion studies
setwd('~/Projects/hog/')
library(tidyverse)
library(magrittr)
library(ggpubr)

figure_root <- 'figures/liti_growth/'

#### Import Data ####
strains <- readRDS('data/Rdata/strain_meta.rds')
distance_to_ref <- structure(strains$`Total number of SNPs`, names=strains$`Standardized name`)

filtered_strains_dip <- filter(strains, Ploidy == 2, Aneuploidies == 'euploid', Zygosity == 'homozygous') %>% pull(`Standardized name`)
filtered_strains_hap <- filter(strains, Ploidy == 1, Aneuploidies == 'euploid') %>% pull(`Standardized name`)
filtered_strains_both <- filter(strains, Ploidy <= 2, Aneuploidies == 'euploid') %>% pull(`Standardized name`)

# Filter strains that are unusually distant from the reference genome
filtered_strains_dip <-  setdiff(filtered_strains_dip , c("AMH", "BAG", "BAH", "BAL", "CEG", "CEI"))
filtered_strains_hap <-  setdiff(filtered_strains_hap , c("AMH", "BAG", "BAH", "BAL", "CEG", "CEI"))
filtered_strains_both <-  setdiff(filtered_strains_both , c("AMH", "BAG", "BAH", "BAL", "CEG", "CEI"))

genes <- readRDS('data/Rdata/gene_meta_all.rds')

complexes <- readRDS('data/Rdata/complex_members.rds')

essential <- readRDS('data/Rdata/essential_genes.rds')
essential_genes <- filter(essential, essential == 'E') %>% pull(locus)
non_essential_genes <- filter(essential, essential == 'NE') %>% pull(locus)

impacts <- readRDS('data/Rdata/all_muts_impacts.rds')

growth <- readRDS('data/Rdata/growth_liti.rds') %>%
  filter(strain %in% filtered_strains_both)

allele_freqs <- readRDS('data/Rdata/allele_freqs.rds')

probs <- readRDS('data/Rdata/paff_all_genes.rds') %>%
  filter(strain %in% growth$strain) %>%
  filter(strain %in% filtered_strains_both) %>%
  left_join(. ,growth, by = 'strain') %>%
  gather(key='condition', value = 'growth', -strain, -gene, -p_aff) %>%
  mutate(gene_sig = FALSE)

probs_norm <- readRDS('data/Rdata/norm_paff.rds') %>%
  rename(p_aff = norm_p_aff) %>%
  filter(strain %in% growth$strain) %>%
  filter(strain %in% filtered_strains_both) %>%
  left_join(. ,growth, by = 'strain') %>%
  gather(key='condition', value = 'growth', -strain, -gene, -p_aff) %>%
  mutate(gene_sig = FALSE)

probs_mat <- readRDS('data/Rdata/paff_all_genes_mat.rds') %>%
  filter(strain %in% growth$strain) %>%
filter(strain %in% filtered_strains_both)

probs_norm_mat <- readRDS('data/Rdata/norm_paff_mat.rds') %>%
  filter(strain %in% growth$strain) %>%
  filter(strain %in% filtered_strains_both)

# import Ko growth data showing which genes are significant (using 0.01 threshold)
ko_growth <- readRDS('data/Rdata/gene_ko_growth.rds') %>%
  filter(gene %in% unique(probs$gene))

ko_growth_full <- readRDS('data/Rdata/gene_ko_growth_full.rds') %>%
  filter(gene %in% unique(probs$gene))

#### Processing ####
ko_thresh <- 0.95
ko_norm_thresh <- mean(probs_norm$p_aff) + qnorm(0.95) * sd(probs_norm$p_aff)

## Define conditions considered equivalent between bede and liti screens
equiv_conditions <- structure(c("NaCl 0.6M (48H)", "39ºC (48H)", "Caffeine 20mM (48H)", "Glycerol 2%  (48H)", "6-AU (48H)"),
                              names=c("ypdnacl15m", "ypd40", "ypdcafein40", "ypglycerol", "ypd6au"))

sig_genes <- lapply(equiv_conditions, function(x){filter(ko_growth, condition == x, num_strains >= 1) %>% pull(gene)})
sig_genes_strong <- lapply(equiv_conditions, function(x){filter(ko_growth, condition == x, num_strains >= 2) %>% pull(gene)})
sig_genes_v_strong <- lapply(equiv_conditions, function(x){filter(ko_growth, condition == x, num_strains >= 3) %>% pull(gene)})

for (con in names(sig_genes)){
  probs[probs$condition == con & probs$gene %in% sig_genes[[con]], 'gene_sig'] <- TRUE
  probs_norm[probs_norm$condition == con & probs_norm$gene %in% sig_genes[[con]], 'gene_sig'] <- TRUE
}

# Sig genes plot


#### Analysis ####
## By strain over all sig genes
strain_summary <- group_by(probs, condition) %>%
  group_by(strain, add=TRUE) %>%
  summarise(num_sig_kos = sum(gene_sig & p_aff > ko_thresh),
            growth = first(growth))

fit <- lm(growth ~ num_sig_kos, data = filter(strain_summary, condition == 'ypdnacl15m'))
p_ypdnacl15m_sig_kos <- ggplot(filter(strain_summary, condition == 'ypdnacl15m'), aes(x=num_sig_kos, y=growth)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  ylab('Growth Relative to YPD Media') +
  xlab(paste0('Number of NaCl Stress Significant Genes with P(Aff) > ', ko_thresh)) #+
  annotate('text', x = 25, y = 0.425, label=paste0('y = ',
                                                 signif(fit$coefficients[2], 4),
                                                 'x + ',
                                                 signif(fit$coefficients[1], 4),
                                                 ', Adj. R-Squared = ',
                                                 signif(summary(fit)$adj.r.squared, 2)))

ggsave(paste0(figure_root, 'nacl15_growth_vs_number_sig_kos.pdf'), p_ypdnacl15m_sig_kos, width = 7, height = 5)

## By individual gene
probs %<>% mutate(ko = p_aff > ko_thresh)
probs_norm %<>% mutate(ko = p_aff > ko_norm_thresh)

# Example with Pbs2 and Hog1 (Known to be sig in NaCl 1.5mM)
hog_id <- 'YLR113W'
pbs2_id <- 'YJL128C'
kos_nacl <- filter(probs, gene==!!hog_id | gene == !!pbs2_id, condition=='ypdnacl15m') %>%
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

levels(kos_nacl$condition) <- c('1.5M NaCl', '14°C')

p_osmotic_shock_ko_growth_box <- ggplot(kos_nacl, aes(x = ko, y = growth, colour = ko)) +
  facet_wrap(~ condition, labeller = label_bquote(cols=.(condition))) +
  geom_boxplot(varwidth = TRUE) +
  xlab(paste0('Genes with P(Aff) > ', ko_thresh)) +
  ylab('Growth Relative to YPD Media') + 
  stat_compare_means(method = 'wilcox.test',
                     comparisons = list(c('Hog1', 'Neither'), c('Pbs2', 'Neither'))) +
  stat_summary(geom = 'text', fun.data = function(x){return(c(y = -0.03, label = length(x)))}) +
  guides(colour = FALSE)
ggsave(paste0(figure_root, 'osmotic_shock_ko_growth_pbs2_hog1.pdf'), p_osmotic_shock_ko_growth_box, width = 7, height = 5)


# Investigate full distribution of P(Aff)'s vs growth
probs_nacl <- filter(probs, gene==!!hog_id | gene == !!pbs2_id, condition=='ypdnacl15m') %>%
  mutate(gene = structure(c('Hog1', 'Pbs2'), names=c(hog_id, pbs2_id))[gene])

probs_norm_nacl <- filter(probs_norm, gene==!!hog_id | gene == !!pbs2_id, condition=='ypdnacl15m') %>%
  mutate(gene = structure(c('Hog1', 'Pbs2'), names=c(hog_id, pbs2_id))[gene])

p_hog_pbs_p_aff_growth <- ggplot(probs_nacl, aes(x=p_aff, y=growth, colour=gene)) +
  geom_point() +
  facet_wrap(~ gene) +
  guides(colour=FALSE) + 
  xlab('P(Aff)') + 
  ylab('Growth in 1.5mM NaCl relative to YPD Media')

p_hog_pbs_p_aff_growth_norm <- ggplot(probs_norm_nacl, aes(x=p_aff, y=growth, colour=gene)) +
  geom_point() +
  facet_wrap(~ gene) + 
  xlab('Normalised P(Aff)') + 
  ylab('Growth in 1.5mM NaCl relative to YPD Media')

p_hog_pbs_p_aff_growth_both <- ggarrange(p_hog_pbs_p_aff_growth, p_hog_pbs_p_aff_growth_norm)
ggsave(paste0(figure_root, 'p_aff_distribution_vs_growth.pdf'), p_hog_pbs_p_aff_growth_both, width = 22, height = 10)

# Generic function for gene ko growth box plots
plot_ko_growth_box <- function(condition, gene_names=NULL, gene_ids=NULL, prob_tbl=probs, ylab=NULL, xlab=NULL){
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
    select(strain, gene, growth, ko) %>%
    spread(key = 'gene', value = 'ko') %>%
    rename_at(.vars = vars(!!gene_ids), .funs = function(x){gene_names[x]}) %>%
    gather(key = 'gene', value = 'ko', -strain, -growth) %>%
    mutate(ko = factor(ifelse(ko, 'True', 'False'), levels = c('False', 'True'))) %>%
    mutate(gene = factor(gene, levels = gene_names))
  
  if (is.null(ylab)){
    ylab = paste0('Growth in ', condition, ' Relative to YPD Media')
  }
  if (is.null(xlab)){
    xlab = paste0('P(Aff) > ', ko_thresh)
  }
  
  p <- ggplot(probs_filter, aes(x=ko, y=growth)) +
    geom_boxplot(aes(fill=gene), varwidth = TRUE, alpha=0.4) +
    facet_wrap(~gene) +
    stat_compare_means(comparisons = list(c('False','True'))) +
    stat_summary(geom = 'text', fun.data = function(x){return(c(y = -0.03, label = length(x)))}) +
    xlab(xlab) +
    ylab(ylab) +
    ylim(-0.05,0.6) +
    guides(fill=FALSE)
  return(p)
}

p_nacl15m_all_sig_kos_box <- plot_ko_growth_box('ypdnacl15m', gene_ids = sig_genes_strong$ypdnacl15m, ylab = 'Relative Growth in 1.5M NaCl')
ggsave(paste0(figure_root, 'osmotic_shock_ko_growth_all_genes.pdf'), p_nacl15m_all_sig_kos_box, width = 7, height = 10)

p_nacl1m_all_sig_kos_box <- plot_ko_growth_box('ypdnacl1m', gene_ids = sig_genes_strong$ypdnacl15m)
ggsave(paste0(figure_root, 'low_osmotic_shock_ko_growth_all_genes.pdf'), p_nacl1m_all_sig_kos_box, width = 20, height = 20)

p_nacl15m_all_sig_vstrong_kos_box <- plot_ko_growth_box('ypdnacl15m', gene_ids = sig_genes_v_strong$ypdnacl15m)
ggsave(paste0(figure_root, 'osmotic_shock_ko_growth_all_genes_strong_sig.pdf'), p_nacl15m_all_sig_vstrong_kos_box, width = 20, height = 20)

p_nacl15m_all_sig_norm_kos_box <- plot_ko_growth_box('ypdnacl15m', gene_ids = sig_genes_strong$ypdnacl15m, prob_tbl = probs_norm)
ggsave(paste0(figure_root, 'normalised_associations/osmotic_shock_ko_growth_all_genes.pdf'), p_nacl15m_all_sig_norm_kos_box, width = 20, height = 20)

p_kcl2m_all_sig_kos_box <- plot_ko_growth_box('ypdkcl2m', gene_ids = sig_genes_strong$ypdnacl15m)
ggsave(paste0(figure_root, 'kcl_osmotic_shock_ko_growth_all_genes.pdf'), p_kcl2m_all_sig_kos_box, width = 20, height = 20)

p_heat_all_sig_kos_box <- plot_ko_growth_box('ypd40', gene_ids = sig_genes_strong$ypd40)
ggsave(paste0(figure_root, 'high_temp_ko_growth_all_genes.pdf'), p_heat_all_sig_kos_box, width = 30, height = 30)

p_cafein40_all_sig_kos_box <- plot_ko_growth_box('ypdcafein40', gene_ids = sig_genes_strong$ypdcafein40)
ggsave(paste0(figure_root, 'cafein40_ko_growth_all_genes.pdf'), p_cafein40_all_sig_kos_box, width = 30, height = 30)

p_glycerol_all_sig_kos_box <- plot_ko_growth_box('ypglycerol', gene_ids = sig_genes_strong$ypglycerol)
ggsave(paste0(figure_root, 'glycerol_ko_growth_all_genes.pdf'), p_glycerol_all_sig_kos_box, width = 30, height = 30)

p_6au_all_sig_kos_box <- plot_ko_growth_box('ypd6au', gene_ids = sig_genes_strong$ypd6au)
ggsave(paste0(figure_root, '6au_ko_growth_all_genes.pdf'), p_6au_all_sig_kos_box, width = 30, height = 30)

# Compare significance to effect size
probs_sig <- filter(probs,
                    condition %in% names(equiv_conditions),
                    gene %in% unique(unlist(sig_genes_strong)))

probs_sig_summary <- probs_sig %>%
  group_by(cg = paste(condition, gene)) %>%
  filter(any(ko)) %>%
  do(w = wilcox.test(growth ~ ko, data = .), condition = first(.$condition), gene = first(.$gene)) %>%
  summarise(condition,
            gene,
            pVal = w$p.value)

ko_growth_sig <- filter(ko_growth, condition %in% equiv_conditions, gene %in% unique(unlist(sig_genes_strong))) %>%
  mutate(maxSscore = pmax(`S288C-score`, `UWOP-score`, `Y55-score`, `YPS-score`, na.rm = TRUE)) %>%
  mutate(condition2 = unname(structure(names(equiv_conditions), names=equiv_conditions)[condition]))

gene_sscore_hash <- structure(ko_growth_sig$maxSscore, names=paste0(ko_growth_sig$condition2, ko_growth_sig$gene))

probs_sig_summary %<>% mutate(maxSscore = unname(gene_sscore_hash[paste0(condition, gene)])) %>%
  filter(!is.na(maxSscore))

p_sscore_p_val <- ggplot(probs_sig_summary, aes(x=maxSscore, y=pVal)) +
  geom_point()
# No Apparant relationship

## Test interactions between NaCl ko genes
probs_mat_nacl <- select(probs_mat, one_of(sig_genes_strong$ypdnacl15m))

# Filter to genes that have real impact
ko_thresh <- 0.5
real_effects <- c('STE50', 'PBS2', 'RXT3', 'PPZ1', 'ALG12')
real_effects_id <- structure(genes$id, names=genes$name)[real_effects]

probs_mat_nacl %<>% select(one_of(real_effects_id))
# sum(rowSums(probs_mat_nacl > 0.95) > 1)
# 11 strains have more than one likely KO in these genes

probs_mat_nacl %<>% mutate(kos = factor(rowSums(. > ko_thresh), levels = c('0', '1', '2', '3'))) %>%
  mutate(strain = growth$strain,
                           nacl15m = growth$ypdnacl15m,
                           nacl1m = growth$ypdnacl1m)
  
p_nacl_interaction <- ggplot(probs_mat_nacl, aes(x=kos, y=nacl15m)) +
  geom_boxplot(varwidth = TRUE) +
  stat_compare_means(comparisons = list(c('2','3'), c('1','2'), c('0','1'), c('1','3'), c('0','2'), c('0','3'))) +
  stat_summary(geom = 'text', fun.data = function(x){return(c(y = -0.03, label = length(x)))}) +
  xlab(paste0('Number of Osmotic Shock Genes with P(aff) > ', ko_thresh)) +
  ylab('Growth in 1.5mM NaCl Relative to YPD Media')
ggsave(paste0(figure_root, 'nacl15_ko_interactions.pdf'), p_nacl_interaction, width = 12, height = 10)

# Same with Cafeine
ko_thresh <- 0.95
probs_mat_caf <- select(probs_mat, one_of(sig_genes_strong$ypdcafein40)) %>%
  mutate(kos = factor(rowSums(. > ko_thresh))) %>%
  mutate(strain = growth$strain, cafein40 = growth$ypdcafein40)

p_cafe_interaction <- ggplot(probs_mat_caf, aes(x=kos, y=cafein40)) +
  geom_boxplot(varwidth = TRUE) +
  #stat_compare_means(comparisons = list(c('2','3'), c('1','2'), c('0','1'), c('1','3'), c('0','2'), c('0','3'))) +
  stat_summary(geom = 'text', fun.data = function(x){return(c(y = -0.03, label = length(x)))}) +
  xlab(paste0('Number of Caffeine Genes with P(aff) > ', ko_thresh)) +
  ylab('Growth in 40mM Caffeine Relative to YPD Media')
ggsave(paste0(figure_root, 'caffeine_ko_interactions.pdf'), p_cafe_interaction, width = 12, height = 10)



## Linear models for growth in condition
fit <- lm(nacl15m ~ ., data = select(probs_mat_nacl, -strain, -nacl1m))
fit_summary <- summary(fit)
top_genes <- as_tibble(fit_summary$coefficients, rownames = 'gene_id') %>%
  set_names(c('gene_id', 'coef', 'std_err', 't_val', 'p_val')) %>%
  mutate(p_adj = p.adjust(p_val)) %>%
  #filter(p_adj < 0.05) %>%
  mutate(gene_name = structure(genes$name, names=genes$id)[gene_id])

## Analyse all genes
all_sig_genes <- unique(unname(unlist(sig_genes_strong)))
gene_summary <- readRDS('data/Rdata/paff_all_genes.rds') %>%
  mutate(sig = gene %in% ko_growth$gene) %>%
  mutate(sig_nacl = gene %in% sig_genes_strong$ypdnacl15m) %>%
  mutate(ko = p_aff > ko_thresh) %>%
  left_join(., growth, by='strain') %>%
  group_by(gene) %>%
  summarise(count = sum(ko),
            growth_ko = mean(ypdnacl15m[ko], na.rm = TRUE),
            growth_not = mean(ypdnacl15m[!ko], na.rm = TRUE),
            sig_nacl = first(sig_nacl),
            sig=first(sig)) %>%
  mutate(essential = structure(essential$essential, names=essential$locus)[gene]) %>%
  mutate(sigess = if_else(essential == 'E', 'Essential', if_else(sig, 'Significant', 'Neither')))
  

# Essential and sig KO genes are mutually exclusive as expected (S-Score would mean no diff between stressed and unstressed growth in essential gene)
table(gene_summary$sig_nacl, gene_summary$essential)

# Not a big difference in rate of KO between genes that are sig in some conditions
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
  gather(key = 'strain', value = 'sscore', -gene, -count, -sig, -essential, -sigess, -sig_nacl, -growth_ko, -growth_not)

p_ko_count_significance <- ggplot(gene_summary, aes(x=sscore, y=count, colour=strain, shape=sig_nacl)) +
  geom_point() +
  xlab('S-Score in 0.6mM NaCl (48hr)') +
  ylab(paste0('Number of Strains with P(Aff) > ', ko_thresh)) +
  guides(colour=guide_legend(title = 'KO Strain for S-Score'), shape=guide_legend(title = 'Gene Significant in 2+ strains'))

ggsave(paste0(figure_root, 'strain_ko_count_vs_significance_nacl.pdf'), p_ko_count_significance, width = 12, height = 10)

p_strain_sscore_density <- ggplot(gene_summary, aes(x=sscore, colour=strain)) +
  geom_density()