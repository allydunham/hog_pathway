# Script to analyse normalisation of P(Aff) scores
setwd('~/Projects/hog/')
library(moments)
library(tidyverse)
library(magrittr)
library(ggpubr)

#### Import Data ####
## Strains
strains <- readRDS('data/Rdata/strain_meta.rds')
filtered_strains <- filter(strains, Ploidy == 2, Aneuploidies == 'euploid') %>% pull(`Standardized name`)
distance_to_ref_full <- structure(strains$`Total number of SNPs`, names=strains$`Standardized name`)

# Filter strains that are unusually distant from the reference genome
filtered_strains <-  setdiff(filtered_strains , c("AMH", "BAG", "BAH", "BAL", "CEG", "CEI"))

## Genes
genes <- readRDS('data/Rdata/gene_meta_all.rds') %>%
  mutate(length = stop - start)
complexes <- readRDS('data/Rdata/complex_members.rds')
essential <- readRDS('data/Rdata/essential_genes.rds')
essential_hash <- structure(essential$essential, names=essential$locus)

## Genotypes
genotypes <- readRDS('data/Rdata/genotypes_all_genes.rds')
distance_to_ref <- colSums(select(genotypes, -mut_id))

## P(Aff)
probs <- readRDS('data/Rdata/paff_all_genes.rds') %>%
  filter(strain %in% filtered_strains) %>%
  mutate(length=unname(structure(genes$length, names=genes$id)[gene])) %>%
  mutate(ref_dist = unname(distance_to_ref_full[strain]))

#### Normalisation ####
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

#### Analysis ####
strain_summary <- group_by(probs, strain) %>% 
  summarise(mean_paff = mean(p_aff),
            mean_ref_cor_paff = mean(p_aff_ref_cor),
            mean_both_cor_paff = mean(p_aff_both_cor)) %>%
  mutate(ref_dist = unname(distance_to_ref_full[strain]))

gene_summary <- group_by(probs, gene) %>%
  summarise(mean_paff = mean(p_aff),
            mean_paff_dist_cor = mean(p_aff_ref_cor),
            mean_paff_dist_length_cor = mean(p_aff_ref_length_cor),
            length = first(length),
            mean_paff_both_cor = mean(p_aff_both_cor))

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
p_before_length_correction <- ggplot(gene_summary, aes(x=length)) + 
  geom_point(aes(y=mean_paff, colour = 'Pre Strain Correction')) +
  geom_point(aes(y=mean_paff_dist_cor, colour = 'Post Strain Correction')) +
  xlab("Gene Length") + 
  ylab("Mean P(Aff)")
# plot also shows ref dist normalisation doesn't affect gene length effect

p_after_length_correction <- ggplot(gene_summary, aes(y=mean_paff_dist_length_cor, x = length)) + 
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


## Analyse essentiality
gene_summary %<>% mutate(essential = unname(essential_hash[gene]))

p_corrected_essential <- ggplot(filter(gene_summary, !is.na(essential)), aes(x=essential, y=mean_paff_both_cor, colour=essential)) +
  geom_boxplot(notch = TRUE, varwidth = TRUE) + 
  xlab('Gene Essential?') +
  ylab('Mean Normalised P(Aff)')

ggsave('figures/paff_checks/normalised_essential_box.pdf', p_corrected_essential, width = 12, height = 10)

t.test(mean_paff_dist_length_cor ~ essential, data = gene_summary)


