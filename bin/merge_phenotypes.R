#!/us/bin/env Rscript
# Merge Phenotype data for Julia
library(tidyverse)

bede <- read_tsv('data/raw/yeasts_liti_fixed.tsv', col_names = TRUE, col_types = cols(strain = col_character())) %>%
  filter(subset == 'liti') %>%
  select(bede_strain_id = strain, strain = info, condition, score, qvalue)

liti <- read_tsv('data/raw/phenoMatrix_35ConditionsNormalizedByYPD.tab') %>% 
  rename(strain = X1) %>%
  pivot_longer(-strain, names_to = "condition", values_to = "score")

write_tsv(bede, "growth_bede.tsv")
write_tsv(liti, "growth_liti.tsv")
