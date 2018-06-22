# Script to analyse normalisation of P(Aff) scores
setwd('~/Projects/hog/')
library(moments)
library(tidyverse)
library(magrittr)
library(ggpubr)

#### Import Data ####
## Strains
strains <- readRDS('data/Rdata/strain_meta.rds')
distance_to_ref_full <- structure(strains$`Total number of SNPs`, names=strains$`Standardized name`)

filtered_strains_dip <- filter(strains, Ploidy == 2, Aneuploidies == 'euploid') %>% pull(`Standardized name`)
filtered_strains_hap <- filter(strains, Ploidy == 1, Aneuploidies == 'euploid') %>% pull(`Standardized name`)
filtered_strains_both <- filter(strains, Ploidy <= 2, Aneuploidies == 'euploid') %>% pull(`Standardized name`)

# Filter strains that are unusually distant from the reference genome
filtered_strains_dip <-  setdiff(filtered_strains_dip , c("AMH", "BAG", "BAH", "BAL", "CEG", "CEI"))
filtered_strains_hap <-  setdiff(filtered_strains_hap , c("AMH", "BAG", "BAH", "BAL", "CEG", "CEI"))
filtered_strains_both <-  setdiff(filtered_strains_both , c("AMH", "BAG", "BAH", "BAL", "CEG", "CEI"))

## Genes
genes <- readRDS('data/Rdata/gene_meta_all.rds') %>%
  mutate(length = stop - start)
complexes <- readRDS('data/Rdata/complex_members.rds')

essential <- readRDS('data/Rdata/essential_genes.rds')
essential_hash <- structure(essential$essential, names=essential$locus)
essential_genes <- filter(essential, essential == 'E') %>% pull(locus)
non_essential_genes <- filter(essential, essential == 'NE') %>% pull(locus)

## Genotypes
genotypes <- readRDS('data/Rdata/genotypes_all_genes.rds')
distance_to_ref <- colSums(select(genotypes, -mut_id))

## P(Aff)
probs <- readRDS('data/Rdata/paff_all_genes.rds') %>%
  filter(strain %in% filtered_strains_both) %>%
  mutate(length=unname(structure(genes$length, names=genes$id)[gene])) %>%
  mutate(ref_dist = unname(distance_to_ref_full[strain]))

worst_probs <- readRDS('data/Rdata/worst_probs_all.rds') %>%
  filter(strain %in% filtered_strains_both)

counts <- readRDS('data/Rdata/all_gene_mut_counts.rds') %>%
  filter(strain %in% filtered_strains_both) %>%
  gather(key = 'gene', value = 'count', -strain)

probs %<>% mutate(count = counts$count) %>%
  mutate(worst_paff = worst_probs$worst_p_aff) %>%
  mutate(essential = NA)

probs[probs$gene %in% essential_genes, 'essential'] <- 'E'
probs[probs$gene %in% non_essential_genes, 'essential'] <- 'NE'

#### Normalisation ####
### Of Overall P(Aff)
## Independently on SC288 distance
fit_ref_dist <- lm(p_aff ~ ref_dist, data = probs)
probs %<>% mutate(p_aff_ref_cor = unname(p_aff - (fit_ref_dist$coefficients[1] + fit_ref_dist$coefficients[2] * ref_dist)))

## Subsequently on gene length
fit_length <- lm(p_aff_ref_cor ~ length, data = probs)
probs %<>% mutate(p_aff_ref_length_cor = p_aff_ref_cor - (fit_length$coefficients[1] + fit_length$coefficients[2] * length))

## Normalise both together
fit_both <- lm(p_aff ~ length + ref_dist, data = probs)
probs %<>% mutate(p_aff_both_cor = p_aff - (fit_both$coefficients[1] + fit_both$coefficients[2] * length + fit_both$coefficients[3] * ref_dist))
saveRDS(fit_both$coefficients, 'data/Rdata/len_ref_dist_normalisation_lm_coefs.rds')

## Double Normalise Worst P(Aff)
fit_both <- lm(worst_paff ~ length + ref_dist, data = probs)
probs %<>% mutate(worst_paff_both_cor = worst_paff - (fit_both$coefficients[1] + fit_both$coefficients[2] * length + fit_both$coefficients[3] * ref_dist))
saveRDS(fit_both$coefficients, 'data/Rdata/worst_len_ref_dist_normalisation_lm_coefs.rds')

## Double Normalise Count
fit_both <- lm(count ~ length + ref_dist, data = probs)
probs %<>% mutate(count_both_cor = count - (fit_both$coefficients[1] + fit_both$coefficients[2] * length + fit_both$coefficients[3] * ref_dist))
saveRDS(fit_both$coefficients, 'data/Rdata/count_len_ref_dist_normalisation_lm_coefs.rds')

#### Analysis ####
strain_summary <- group_by(probs, strain) %>% 
  summarise(mean_paff = mean(p_aff),
            mean_ref_cor_paff = mean(p_aff_ref_cor),
            mean_both_cor_paff = mean(p_aff_both_cor),
            mean_worst_paff = mean(worst_paff),
            mean_worst_paff_cor = mean(worst_paff_both_cor),
            mean_count = mean(count),
            mean_count_cor = mean(count_both_cor)) %>%
  mutate(ref_dist = unname(distance_to_ref_full[strain]))

gene_summary_hap <- filter(probs, strain %in% filtered_strains_hap) %>%
  group_by(gene) %>%
  summarise(mean_paff = mean(p_aff),
            mean_paff_dist_cor = mean(p_aff_ref_cor),
            mean_paff_dist_length_cor = mean(p_aff_ref_length_cor),
            length = first(length),
            mean_paff_both_cor = mean(p_aff_both_cor),
            mean_worst_paff = mean(worst_paff),
            mean_worst_paff_cor = mean(worst_paff_both_cor),
            mean_count = mean(count),
            mean_count_cor = mean(count_both_cor)) %>%
  mutate(essential = unname(essential_hash[gene]))

gene_summary_dip <- filter(probs, strain %in% filtered_strains_dip) %>%
  group_by(gene) %>%
  summarise(mean_paff = mean(p_aff),
            mean_paff_dist_cor = mean(p_aff_ref_cor),
            mean_paff_dist_length_cor = mean(p_aff_ref_length_cor),
            length = first(length),
            mean_paff_both_cor = mean(p_aff_both_cor),
            mean_worst_paff = mean(worst_paff),
            mean_worst_paff_cor = mean(worst_paff_both_cor),
            mean_count = mean(count),
            mean_count_cor = mean(count_both_cor)) %>%
  mutate(essential = unname(essential_hash[gene]))

gene_summary_both <- group_by(probs, gene) %>%
  summarise(mean_paff = mean(p_aff),
            mean_paff_dist_cor = mean(p_aff_ref_cor),
            mean_paff_dist_length_cor = mean(p_aff_ref_length_cor),
            length = first(length),
            mean_paff_both_cor = mean(p_aff_both_cor),
            mean_worst_paff = mean(worst_paff),
            mean_worst_paff_cor = mean(worst_paff_both_cor),
            mean_count = mean(count),
            mean_count_cor = mean(count_both_cor)) %>%
  mutate(essential = unname(essential_hash[gene])) %>%
  mutate(essential = as.factor(essential))
levels(gene_summary_both$essential) <- c('Essential', 'Nonessential')

# Just ref dist normalisation
p_before_ref_correction <- ggplot(strain_summary, aes(x=ref_dist, y=mean_paff)) + 
  geom_point() +
  xlab("Differences from SC288") + 
  ylab("Mean P(Aff)") + 
  geom_smooth(method = 'lm')

p_after_ref_correction <- ggplot(strain_summary, aes(x=ref_dist, y=mean_ref_cor_paff)) + 
  geom_point() +
  xlab("Differences from SC288") + 
  ylab("Mean P(Aff)") + 
  geom_smooth(method = 'lm')

# Subsequent length normalisation
p_before_length_correction <- ggplot(gene_summary_both, aes(x=length)) + 
  geom_point(aes(y=mean_paff, colour = 'Pre Strain Correction')) +
  geom_point(aes(y=mean_paff_dist_cor, colour = 'Post Strain Correction')) +
  xlab("Gene Length") + 
  ylab("Mean P(Aff)")
# plot also shows ref dist normalisation doesn't affect gene length effect

p_after_length_correction <- ggplot(gene_summary_both, aes(y=mean_paff_dist_length_cor, x = length)) + 
  geom_point() +
  geom_smooth(method = 'lm') +
  xlab("Gene Length") + 
  ylab("Mean Length Corrected P(Aff)")

# Double normalisation
p_both_length <- ggplot(gene_summary, aes(y=mean_paff_both_cor, x = length)) + 
  geom_point() +
  geom_smooth(method = 'lm') +
  xlab("Gene Length") + 
  ylab("Mean Length Corrected P(Aff)")

p_both_ref <- ggplot(strain_summary, aes(x=ref_dist, y=mean_both_cor_paff)) + 
  geom_point() +
  xlab("Differences from SC288") + 
  ylab("Mean P(Aff)") + 
  geom_smooth(method = 'lm')

ggsave('figures/paff_checks/normalised_length.pdf', p_both_length, width = 12, height = 10)
ggsave('figures/paff_checks/normalised_ref_dist.pdf', p_both_ref, width = 12, height = 10)

# Behaves identically - as it must?

## Worst P(Aff)
p_worst_ref_before_correction <- ggplot(strain_summary, aes(x=ref_dist, y=mean_worst_paff)) + 
  geom_point() +
  xlab("Differences from SC288") + 
  ylab("Mean P(Aff)") + 
  geom_smooth(method = 'lm')

p_worst_length_before_correction <- ggplot(gene_summary, aes(y=mean_worst_paff, x = length)) + 
  geom_point() +
  geom_smooth(method = 'lm') +
  xlab("Gene Length") + 
  ylab("Mean Length Corrected P(Aff)")

p_worst_both_length <- ggplot(gene_summary, aes(y=mean_worst_paff_cor, x = length)) + 
  geom_point() +
  geom_smooth(method = 'lm') +
  xlab("Gene Length") + 
  ylab("Mean Length Corrected P(Aff)")

p_worst_both_ref <- ggplot(strain_summary, aes(x=ref_dist, y=mean_worst_paff_cor)) + 
  geom_point() +
  xlab("Differences from SC288") + 
  ylab("Mean P(Aff)") + 
  geom_smooth(method = 'lm')

## Count
p_count_ref_before_correction <- ggplot(strain_summary, aes(x=ref_dist, y=mean_count)) + 
  geom_point() +
  xlab("Differences from SC288") + 
  ylab("Mean P(Aff)") + 
  geom_smooth(method = 'lm')

p_count_length_before_correction <- ggplot(gene_summary, aes(y=mean_count, x = length)) + 
  geom_point() +
  geom_smooth(method = 'lm') +
  xlab("Gene Length") + 
  ylab("Mean Length Corrected P(Aff)")

p_count_both_length <- ggplot(gene_summary, aes(y=mean_count_cor, x = length)) + 
  geom_point() +
  geom_smooth(method = 'lm') +
  xlab("Gene Length") + 
  ylab("Mean Length Corrected P(Aff)")

p_count_both_ref <- ggplot(strain_summary, aes(x=ref_dist, y=mean_count_cor)) + 
  geom_point() +
  xlab("Differences from SC288") + 
  ylab("Mean P(Aff)") + 
  geom_smooth(method = 'lm')


## Analyse essentiality
p_corrected_essential <- ggplot(filter(gene_summary_both, !is.na(essential)), aes(x=essential, y=mean_paff_both_cor, colour=essential)) +
  geom_boxplot(notch = TRUE, varwidth = TRUE) + 
  ylim(-0.25, 0.35) +
  xlab('') +
  ylab('Mean Normalised P(Aff)') +
  stat_compare_means(method = 'wilcox.test', comparisons = list(c('Essential', 'Nonessential'))) +
#  stat_summary(geom = 'text', fun.data = function(x){return(c(y = -0.26, label = length(x)))}) +
  stat_summary(geom ="text", fun.data = function(x){return(c(y = mean(x), label = signif(mean(x), digits = 3)))}, color="black") +
  guides(colour = FALSE)
ggsave('figures/paff_checks/normalised_essential_box.pdf', p_corrected_essential, width = 12, height = 10)

p_corrected_essential_all <- ggplot(filter(probs, !is.na(essential)), aes(x=essential, y=p_aff_both_cor, colour=essential)) +
  geom_boxplot(notch = TRUE, varwidth = TRUE) + 
  xlab('Gene Essential?') +
  ylab('Normalised P(Aff)') +
  stat_compare_means(method = 'wilcox.test', comparisons = list(c('E', 'NE'))) +
  stat_summary(geom = 'text', fun.data = function(x){return(c(y = -0.75, label = length(x)))}) +
  stat_summary(geom ="text", fun.data = function(x){return(c(y = mean(x), label = signif(mean(x), digits = 3)))}, color="black") +
  guides(colour = FALSE) +
  ggtitle('Using euploid diploid and haploid strains')
ggsave('figures/paff_checks/normalised_essential_all_strains_box.jpg', p_corrected_essential_all, width = 12, height = 10)

p_corrected_worst_essential <- ggplot(filter(gene_summary_both, !is.na(essential)), aes(x=essential, y=mean_worst_paff_cor, colour=essential)) +
  geom_boxplot(notch = TRUE, varwidth = TRUE) + 
  xlab('Gene Essential?') +
  ylab('Mean Normalised Worst P(Aff)') +
  stat_compare_means(method = 'wilcox.test', comparisons = list(c('E', 'NE'))) +
  stat_summary(geom = 'text', fun.data = function(x){return(c(y = -0.26, label = length(x)))}) +
  stat_summary(geom ="text", fun.data = function(x){return(c(y = mean(x), label = signif(mean(x), digits = 3)))}, color="black") +
  guides(colour = FALSE)
ggsave('figures/paff_checks/normalised_worst_essential_box.pdf', p_corrected_worst_essential, width = 12, height = 10)
t.test(mean_worst_paff_cor ~ essential, data = gene_summary)

p_corrected_count_essential <- ggplot(filter(gene_summary_both, !is.na(essential)), aes(x=essential, y=mean_count_cor, colour=essential)) +
  geom_boxplot(notch = TRUE, varwidth = TRUE) + 
  ylim(-1.5, 2.6) +
  xlab('') +
  ylab('Mean Normalised Variant Count') +
  stat_compare_means(method = 'wilcox.test', comparisons = list(c('Essential', 'Nonessential'))) +
#  stat_summary(geom = 'text', fun.data = function(x){return(c(y = -1.5, label = length(x)))}) +
  stat_summary(geom ="text", fun.data = function(x){return(c(y = mean(x), label = signif(mean(x), digits = 3)))}, color="black") +
  guides(colour = FALSE)
ggsave('figures/paff_checks/normalised_count_essential_box.pdf', p_corrected_count_essential, width = 12, height = 10)
t.test(mean_count_cor ~ essential, data = gene_summary)

# Essential presentation figure (requires sift and soldX vs essential figures loaded from )
p <- ggarrange(p_corrected_count_essential, p_corrected_essential, p_essential_sift_per_base, p_essential_foldx_per_base, ncol = 2, nrow = 2, labels = 'auto')
ggsave('figures/paff_checks/essential_genes_boxes.pdf', p, width = 7, height = 7)

## Filter by haploid/diplod strains
p_corrected_essential_dip <- ggplot(filter(gene_summary_dip, !is.na(essential)), aes(x=essential, y=mean_paff_both_cor, colour=essential)) +
  geom_boxplot(notch = TRUE, varwidth = TRUE) + 
  xlab('Gene Essential?') +
  ylab('Mean Normalised P(Aff)') +
  stat_compare_means(method = 'wilcox.test', comparisons = list(c('E', 'NE'))) +
  stat_summary(geom = 'text', fun.data = function(x){return(c(y = -0.26, label = length(x)))}) +
  stat_summary(geom ="text", fun.data = function(x){return(c(y = mean(x), label = signif(mean(x), digits = 3)))}, color="black") +
  guides(colour = FALSE) +
  ggtitle('Using euploid diploid strains')
ggsave('figures/paff_checks/normalised_essential_box_dip.pdf', p_corrected_essential_dip, width = 12, height = 10)

p_corrected_essential_hap <- ggplot(filter(gene_summary_hap, !is.na(essential)), aes(x=essential, y=mean_paff_both_cor, colour=essential)) +
  geom_boxplot(notch = TRUE, varwidth = TRUE) + 
  xlab('Gene Essential?') +
  ylab('Mean Normalised P(Aff)') +
  stat_compare_means(method = 'wilcox.test', comparisons = list(c('E', 'NE'))) +
  stat_summary(geom = 'text', fun.data = function(x){return(c(y = -0.26, label = length(x)))}) +
  stat_summary(geom ="text", fun.data = function(x){return(c(y = mean(x), label = signif(mean(x), digits = 3)))}, color="black") +
  guides(colour = FALSE) +
  ggtitle('Using euploid haploid strains')
ggsave('figures/paff_checks/normalised_essential_box_hap.pdf', p_corrected_essential_hap, width = 12, height = 10)


p_corrected_essential_all_dip <- ggplot(filter(probs, !is.na(essential), strain %in% filtered_strains_dip),
                                        aes(x=essential, y=p_aff_both_cor, colour=essential)) +
  geom_boxplot(notch = TRUE, varwidth = TRUE) + 
  xlab('Gene Essential?') +
  ylab('Normalised P(Aff)') +
  stat_compare_means(method = 'wilcox.test', comparisons = list(c('E', 'NE'))) +
  stat_summary(geom = 'text', fun.data = function(x){return(c(y = -0.75, label = length(x)))}) +
  stat_summary(geom ="text", fun.data = function(x){return(c(y = mean(x), label = signif(mean(x), digits = 3)))}, color="black") +
  guides(colour = FALSE) +
  ggtitle('Using euploid diploid strains')
ggsave('figures/paff_checks/normalised_essential_all_dip_strains_box.jpg', p_corrected_essential_all_dip, width = 12, height = 10)

p_corrected_essential_all_hap <- ggplot(filter(probs, !is.na(essential), strain %in% filtered_strains_hap),
                                        aes(x=essential, y=p_aff_both_cor, colour=essential)) +
  geom_boxplot(notch = TRUE, varwidth = TRUE) + 
  xlab('Gene Essential?') +
  ylab('Normalised P(Aff)') +
  stat_compare_means(method = 'wilcox.test', comparisons = list(c('E', 'NE'))) +
  stat_summary(geom = 'text', fun.data = function(x){return(c(y = -0.75, label = length(x)))}) +
  stat_summary(geom ="text", fun.data = function(x){return(c(y = mean(x), label = signif(mean(x), digits = 3)))}, color="black") +
  guides(colour = FALSE) +
  ggtitle('Using euploid haploid strains')
ggsave('figures/paff_checks/normalised_essential_all_hap_strains_box.jpg', p_corrected_essential_all_hap, width = 12, height = 10)

