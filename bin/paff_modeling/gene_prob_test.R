# Script to check P(Aff) scores have been calculated correctly
setwd('~/Projects/hog/')
library(tidyverse)
library(magrittr)

gene_info <- readRDS('data/Rdata/gene_meta_hog.rds')

impacts <- read_tsv('data/hog-gene-variants.impact', col_names = TRUE, na = "NA",
                      col_types = cols(foldx_int_ddG = col_double(), foldx_int_ddG_sd = col_double()))

genotypes <- read_tsv('data/hog-gene-variants.low-freq-genotypes', col_names = TRUE)

probs <- read_tsv('data/hog-gene-variants.probs', col_names = TRUE)

probs_old <- readRDS('data/Rdata/paff_hog_genes.rds') %>%
  mutate(gene=structure(gene_info$id, names=gene_info$name)[gene]) %>%
  spread(key = 'gene', value = 'p_aff')

probs <- select(probs, names(probs_old))

pneut <- function(x){
  if (x[['type']] == 'nonsynonymous'){
    return(min(1,
               1/(1 + exp(0.21786182 * as.numeric(x[['foldx_ddG']]) + 0.07351653)),
               1/(1 + exp(-1.312424 * log(as.numeric(x[['sift_score']]) + 1.598027e-05) - 4.103955)),
               na.rm = TRUE)
           )
  } else if (x[['type']] %in% c('nonsense', 'frameshift')){
    if (as.numeric(x[['prop_aa']]) > 0.95){
      return(0.01)
    } else {
      return(0.95)
    }
  } else {
    return(1)
  }
}

impacts %<>% mutate(p_neut = apply(impacts, 1, pneut))

test_strain <- 'AAA'
muts <- unlist(genotypes[genotypes['AAA'] > 0, 'mut_id'])

imps <- filter(impacts, mut_id %in% muts) %>%
  group_by(gene) %>%
  summarise(p_aff = 1 - prod(p_neut))

## New method gives expected results but slightly different to old table

