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
  mutate(distance_to_ref = distance_to_ref[strain]) %>% 
  mutate(length=structure(genes$length, names=genes$id)[gene]) %>%
  mutate(ref_dist = distance_to_ref[strain])

#### Normalisation ####
## Independently on SC288 distance
fit_ref_dist <- lm(p_aff ~ ref_dist, data = probs)

probs %<>% mutate(p_aff_ref_cor = p_aff - (fit_ref_dist$coefficients[1] + fit_ref_dist$coefficients[2] * ref_dist))

strain_summary <- group_by(probs, strain) %>% 
  summarise(mean_paff = mean(p_aff), mean_cor_paff = mean(p_aff_ref_cor)) %>%
  mutate(ref_dist = distance_to_ref[strain])

p_before_ref_correction <- ggplot(strain_summary, aes(x=ref_dist, y=mean_paff)) + 
  geom_point() +
  xlab("Differences from SC288") + 
  ylab("Mean P(Aff)") + 
  geom_smooth(method = 'lm')

p_after_ref_correction <- ggplot(strain_summary, aes(x=ref_dist, y=mean_cor_paff)) + 
  geom_point() +
  xlab("Differences from SC288") + 
  ylab("Mean P(Aff)") + 
  geom_smooth(method = 'lm')

## Then normalise on gene length
fit_length <- lm(p_aff_ref_cor ~ length, data = probs)
probs %<>% mutate(p_aff_ref_length_cor = p_aff_ref_cor - (fit_length$coefficients[1] + fit_length$coefficients[2] * length))

gene_summary <- group_by(probs, gene) %>%
  summarise(mean_paff = mean(p_aff),
            mean_paff_dist_cor = mean(p_aff_ref_cor),
            mean_paff_dist_length_cor = mean(p_aff_ref_length_cor),
            length = first(length))

p_before_length_correction <- ggplot(gene_summary, aes(x=length)) + 
  geom_point(aes(y=mean_paff, colour = 'Pre Strain Correction')) +
  geom_point(aes(y=mean_paff_dist_cor, colour = 'Post Strain Correction')) +
  xlab("Gene Length") + 
  ylab("Mean P(Aff)")

p_after_length_correction <- ggplot(gene_summary, aes(y=mean_paff_dist_length_cor, x = length)) + 
  geom_point() +
  geom_smooth(method = 'lm') +
  xlab("Gene Length") + 
  ylab("Mean Length Corrected P(Aff)")

#### Analysis ####






