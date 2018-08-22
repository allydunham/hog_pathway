# Script performing analyses on the growth data from the liti paper with respect to normalised P(Aff) Scores
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

gene_gene_cor <- readRDS('data/Rdata/all_gene_correlations_matrix.rds')

worst_paff <- readRDS('data/Rdata/norm_worst.rds') %>%
  filter(strain %in% intersect(filtered_strains, growth$strain))

counts <- readRDS('data/Rdata/norm_counts.rds') %>%
  filter(strain %in% growth$strain)

gene_scores <- readRDS('data/Rdata/norm_paff.rds') %>%
  filter(strain %in% filtered_strains, strain %in% growth$strain) %>%
  mutate(distance_to_ref = unname(distance_to_ref[strain])) %>% 
  mutate(length=unname(structure(genes$stop - genes$start, names=genes$id)[gene])) %>%
  mutate(count = counts$norm_count) %>%
  mutate(worst_p_aff = worst_paff$norm_worst_p_aff) %>%
  mutate(essential = NA)
gene_scores[gene_scores$gene %in% essential_genes, 'essential'] <- 'E'
gene_scores[gene_scores$gene %in% non_essential_genes, 'essential'] <- 'NE'

gene_scores <- left_join(gene_scores, growth, by = 'strain')

#### Analysis ####
## Growth PCA
pca <- prcomp(t(select(growth, -strain)), center = TRUE, scale. = TRUE)
pr_comps <- as_tibble(pca$x, rownames = 'condition')

# Some broad PCA clustering but nothing huge
p_growth_condition_pca <- ggplot(pr_comps, aes(x=PC1, y=PC2, label=condition)) + 
  geom_text(size=2.5)

## Unusual Gene scores
gene_scores %<>% mutate(norm_p_aff_high = norm_p_aff > mean(.$norm_p_aff) + 3 * sd(.$norm_p_aff)) %>%
  mutate(count_high = count > mean(.$count) + 3 * sd(.$count)) %>%
  mutate(worst_p_aff_high = worst_p_aff > mean(.$worst_p_aff) + 3 * sd(.$worst_p_aff))

essential_summary <- group_by(gene_scores, gene) %>%
  summarise(count_paff = sum(norm_p_aff_high),
            count_worst = sum(worst_p_aff_high),
            count_count = sum(count_high),
            essential = first(essential)) %>%
  filter(!is.na(essential))

p_essential_unusual_norm_paff <- ggplot(essential_summary, aes(x=essential, y=count_paff, colour=essential)) +
  geom_boxplot(varwidth = TRUE, notch = TRUE) +
  xlab('') +
  ylab('Number of Strains with Normalised P(Aff) > 3sd') +
  stat_compare_means(method = 'wilcox.test', comparisons = list(c('E', 'NE'))) +
  stat_summary(geom = 'text', fun.data = function(x){return(c(y = -25, label = length(x)))}) +
  stat_summary(geom ="text", fun.data = function(x){return(c(y = mean(x), label = signif(mean(x), digits = 3)))}, color="black") +
  guides(colour = FALSE)

p_essential_unusual_worst <- ggplot(essential_summary, aes(x=essential, y=count_worst, colour=essential)) +
  geom_boxplot(varwidth = TRUE, notch = TRUE) +
  xlab('') +
  ylab('Number of Strains with Normalised Worst P(Aff) > 3sd') +
  stat_compare_means(method = 'wilcox.test', comparisons = list(c('E', 'NE'))) +
  stat_summary(geom = 'text', fun.data = function(x){return(c(y = -25, label = length(x)))}) +
  stat_summary(geom ="text", fun.data = function(x){return(c(y = mean(x), label = signif(mean(x), digits = 3)))}, color="black") +
  guides(colour = FALSE)

p_essential_unusual_count <- ggplot(essential_summary, aes(x=essential, y=count_count, colour=essential)) + 
  geom_boxplot(varwidth = TRUE, notch = TRUE) +
  xlab('') +
  ylab('Number of Strains with Variant Count > 3sd') +
  stat_compare_means(method = 'wilcox.test', comparisons = list(c('E', 'NE'))) +
  stat_summary(geom = 'text', fun.data = function(x){return(c(y = -25, label = length(x)))}) +
  stat_summary(geom ="text", fun.data = function(x){return(c(y = mean(x), label = signif(mean(x), digits = 3)))}, color="black") +
  guides(colour = FALSE)

p_essential_unusual <- ggarrange(p_essential_unusual_norm_paff, p_essential_unusual_worst, p_essential_unusual_count)
ggsave('figures/paff_checks/unusual_score_counts_norm_essential_box.pdf', p_essential_unusual, width = 20, height = 20)

# find significant difference between essential and non-essential genes but not huge magnitude
t.test(count_count ~ essential, data = essential_summary)
t.test(count_worst ~ essential, data = essential_summary)
t.test(count_paff ~ essential, data = essential_summary)

## Strain Growth Analysis
strain_summary <- group_by(filter(gene_scores, !is.na(essential)), strain) %>%
  summarise(genes_high_paff = sum(norm_p_aff_high),
            genes_high_paff_ess = sum(norm_p_aff_high & essential == 'E'),
            genes_high_worst = sum(worst_p_aff_high),
            genes_high_worst_ess = sum(worst_p_aff_high & essential == 'E'),
            genes_high_count = sum(count_high),
            genes_high_count_ess = sum(count_high & essential == 'E')) %>%
  left_join(. , growth, by='strain') %>%
  rowwise() %>%
  mutate(growth_mean = mean(c(ypacetate:ypdformamide5)))

# No strong relationship between number of apparantly highly variant genes and mean stressed growth - need to look for more specific associations
p_growth_high_paff_mean <- ggplot(strain_summary, aes(x=genes_high_count, y=growth_mean)) +
  geom_point() +
  geom_smooth(method = 'lm')

strain_summary_melt <- gather(strain_summary, key = 'condition', value = 'growth', starts_with('yp')) %>%
  mutate(bin = cut(genes_high_count_ess, breaks = 10))

bin_summary <- group_by(strain_summary_melt, bin) %>%
  summarise(mean=mean(growth),
            sd=sd(growth))

p_growth_high_point <- ggplot(strain_summary_melt, aes(x=genes_high_count_ess, y=growth, colour=condition)) +
  geom_point() +
  xlab('Number of genes with >3 s.d. normalised variants in them') +
  ylab('Growth relative to YPD media')

p_growth_high_bin <- ggplot(strain_summary_melt, aes(x=bin, y=growth)) +
  geom_boxplot(notch = TRUE, varwidth = TRUE) +
  xlab('Number of genes with >3 s.d. normalised variants in them') +
  ylab('')

p_growth_high <- ggarrange(p_growth_high_point, p_growth_high_bin, widths = c(2,1))
ggsave('figures/liti_growth/growth_vs_unusual_variant_count_genes.pdf', p_growth_high, width = 20, height = 10)
## Does not appear to be a big decrease in growth in stressed media with any of the three normalised counts for gene impact

