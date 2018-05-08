# Script performing analyses on the growth data from the liti paper with respect to normalised P(Aff) Scores
setwd('~/Projects/hog/')
library(tidyverse)
library(magrittr)

#### Import Data ####
strains <- readRDS('data/Rdata/strain_meta.rds')
filtered_strains <- filter(strains, Ploidy == 2, Aneuploidies == 'euploid') %>% pull(`Standardized name`)
# Filter strains that are unusually distant from the reference genome
filtered_strains <-  setdiff(filtered_strains , c("AMH", "BAG", "BAH", "BAL", "CEG", "CEI"))
distance_to_ref <- structure(strains$`Total number of SNPs`, names=strains$`Standardized name`)

genes <- readRDS('data/Rdata/gene_meta_all.rds')

complexes <- readRDS('data/Rdata/complex_members.rds')

essential <- readRDS('data/Rdata/essential_genes.rds')
essential_hash <- structure(essential$essential, names=essential$locus)

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
  mutate(essential = unname(essential_hash[gene]))


#### Analysis ####

