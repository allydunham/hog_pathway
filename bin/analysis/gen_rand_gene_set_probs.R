### Script to generate a set of random gene sets and determine the expected number of sig genes in each strain/condition
library(rlang)
library(tidyverse)
library(magrittr)

sig_thresh <- 0.01
num_genes_per_set <- 25 

## Load growth data
ko_growth <- read_tsv('data/raw/ko_scores.txt', col_names = TRUE) %>%
  filter(!gene == 'WT') %>%
  filter(!duplicated(.[,c('strain', 'condition', 'gene')])) %>%
  select(-position) %>%
  mutate(condition = gsub('  ', ' ', condition)) %>% # Some conditions have double spaces in names
  mutate(name = if_else(is.na(name), gene, name))

genes <- unique(ko_growth$name)

## Determine expected number of significant genes for a set of size n in given strain/condition
random_gene_sets <- lapply(10:700, function(x){replicate(num_genes_per_set, sample(genes, x), simplify = FALSE)}) %>% unlist(recursive = FALSE)
names(random_gene_sets) <- paste0('randGroup', 1:length(random_gene_sets))
random_set_lengths <- sapply(random_gene_sets, length)

random_strain_con_num_sig <- data_frame(gene_set = rep(names(random_gene_sets), random_set_lengths),
                                        name = unname(unlist(random_gene_sets)),
                                        set_size = rep(random_set_lengths, random_set_lengths)) %>%
  left_join(., ko_growth, by='name') %>%
  mutate(abs_score = abs(score)) %>%
  group_by(strain, condition, gene_set) %>%
  summarise(set_size=first(set_size),
            n_sig = sum(qvalue < sig_thresh, na.rm = TRUE),
            prop_sig = n_sig/set_size) %>%
  group_by(set_size, add = TRUE) %>%
  summarise(mean_sig = mean(n_sig))

write_tsv(random_strain_con_num_sig, 'rand_gene_set_expected_sig_counts.tsv')
