#!/usr/bin/env Rscript
# Script to normalise P(Aff) scores against gene length and strain distance from S288C
library(tidyverse)
library(magrittr)

args <- commandArgs(trailingOnly=TRUE)
probs_file <- args[1]

out_file <- paste0(probs_file, '.norm')

#### Import ####
model <- readRDS('data/Rdata/len_ref_dist_normalisation_lm_coefs.rds')

strains <- readRDS('data/Rdata/strain_meta.rds')
distance_to_ref <- structure(strains$`Total number of SNPs`, names=strains$`Standardized name`)

genes <- readRDS('data/Rdata/gene_meta_all.rds')
len <- structure(genes$stop - genes$start, names=genes$id)

probs <- read_tsv(probs_file, col_names = TRUE) %>%
  gather(key = 'gene', value = 'p_aff', -strain) %>%
  mutate(length = unname(len[gene])) %>%
  mutate(distance_to_ref = unname(distance_to_ref[strain]))

#### Normalise ####
## model is of the form: E[P(Aff)] = a + b * GeneLength + c * RefDist
probs %<>% mutate(p_aff = p_aff - (model[1] + model[2] * length + model[3] * distance_to_ref)) %>%
  select(strain, gene, p_aff) %>%
  spread(key = 'gene', value = 'p_aff')

#### Save ####
write_tsv(probs, out_file, col_names = TRUE)
