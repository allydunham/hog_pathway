#!/usr/bin/env Rscript
# Script to normalise P(Aff) scores against gene length and strain distance from S288C
# Uses a linear model based on P(Aff) calculated from FoldX and SIFT scores of all variants in a gene to factor out effect of gene length and strain diversity
library(tidyverse)
library(magrittr)

args <- commandArgs(trailingOnly=TRUE)
probs_file <- args[1]
model_file <- args[2]

out_file <- paste0(probs_file, '.norm')

#### Import ####
##Model expected to be a coefficients vector of the form seen on line 29
model <- readRDS(model_file)

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
