#!/usr/bin/env Rscript
# Script to calculate gene-gene and gene-growth correlations
library(tidyverse)
library(magrittr)

args <- commandArgs(trailingOnly = TRUE)

meta_file <- args[1]
gene_file <- args[2]
growth_file <- args[3]
prob_file <- args[4]
complex_file <- args[5]

#### Import ####
## Strain information
meta <- read_tsv(meta_file, col_names = TRUE)

## Gene information
genes <- read_tsv(gene_file, col_names = TRUE)

sys_to_gene <- structure(genes$name, names=genes$id)

filtered_strains <- filter(meta, Ploidy == 2, Aneuploidies == 'euploid') %>% pull(`Standardized name`)

# Growth data
growth <- read_tsv(file = growth_file, col_names = TRUE) %>% 
  rename(strain=X1) %>%
  set_names(str_to_lower(names(.))) %>%
  filter(strain %in% filtered_strains)

## P(aff)
probs <- read_tsv(prob_file, col_names = TRUE) %>%
  filter(strain %in% growth$strain)

# Filter genes with no variation and save list
p_aff_sums <- colSums(select(probs, -strain))
no_prob_genes <- names(which(p_aff_sums == 0))
probs %<>% select(-one_of(no_prob_genes))

write_lines(no_prob_genes, 'cors_zero_pAff_genes')

# Complex Membership
complexes <- read_tsv(complex_file, col_names = TRUE)

#### Calculate values ####
## Gene/Gene Correlation
gene_gene_cor <- cor(select(probs, -strain)) %>%
  as_tibble(rownames='gene')

write_tsv(gene_gene_cor, 'cors_gene_gene_cor.tsv', col_names = TRUE)

same_complex <- function(x, y, comp=complexes){
  x_comps <- comp[comp$ORF == x,] %>% pull(Complex)
  y_comps <- comp[comp$ORF == y,] %>% pull(Complex)
  return(any(x_comps %in% y_comps))
}

gene_gene_cor[lower.tri(gene_gene_cor, diag = TRUE)] <- NA

cor_genes_melt <- gather(gene_gene_cor, key = 'gene2', value = 'cor', -gene) %>%
  drop_na(cor) %>%
  mutate(mag = abs(cor)) %>%
  arrange(desc(mag)) %>%
  mutate(name = sys_to_gene[gene], name2 = sys_to_gene[gene2]) %>%
  mutate(complex = map2_lgl(.$gene, .$gene2, same_complex))

write_tsv(cor_genes_melt, 'cors_gene_gene_cor_melt.tsv', col_names = TRUE)

## Gene/Growth Correlation
cor_growth <- cor(select(probs, -strain), select(growth, -strain)) %>%
  as_tibble(rownames='gene_id')

write_tsv(cor_growth, 'cors_growth_gene_cor.tsv', col_names = TRUE)

cor_growth_melt <- as_tibble(cor_growth, rownames = 'gene_id') %>%
  gather(key='condition', value = 'correlation', -gene_id) %>%
  group_by(condition) %>%
  mutate(gene_name = sys_to_gene[gene_id]) %>%
  mutate(mag = abs(correlation)) %>% 
  mutate(cor_sub_mean = correlation - meanGrowthCor) %>%
  mutate(cor_div_mean = correlation/meanGrowthCor) %>%
  arrange(desc(mag), .by_group = TRUE)

write_tsv(cor_growth_melt, 'cors_growth_gene_cor_melt.tsv', col_names = TRUE)
