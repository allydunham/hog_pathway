#!/us/bin/env Rscript
# Determine P(aff) / gene growth associations
library(tidyverse)
strains <- readRDS('data/Rdata/strain_meta.rds')
distance_to_ref <- structure(strains$`Total number of SNPs`, names=strains$`Standardized name`)

filtered_strains <- filter(strains, Ploidy <= 2, Aneuploidies == 'euploid') %>%
  pull(`Standardized name`) %>%
  setdiff(. , c("AMH", "BAG", "BAH", "BAL", "CEG", "CEI"))

genes <- readRDS('data/Rdata/gene_meta_all.rds')

growth <- read_tsv('data/raw/yeasts_liti_fixed.tsv', col_names = TRUE, col_types = cols(strain = col_character())) %>%
  filter(info %in% filtered_strains) %>% 
  filter(subset == 'liti') %>%
  rename(strain_id = strain) %>%
  rename(strain = info) %>%
  select(strain, condition, score, qvalue)

ko_thresh <- 0.5
probs <- readRDS('data/Rdata/paff_all_genes.rds') %>%
  filter(strain %in% growth$strain) %>%
  mutate(ko = p_aff > ko_thresh)

associations <- left_join(probs, growth, by = "strain") %>%
  group_by(gene, condition) %>%
  filter(!sum(ko) == 0) %>%
  summarise(p_value = wilcox.test(score[ko], score[!ko])$p.value,
            mean_active = mean(score[!ko], na.rm = TRUE),
            mean_ko = mean(score[ko], na.rm = TRUE),
            n = n(),
            n_ko = sum(ko),
            .groups = "drop") %>%
  mutate(n_active = n - n_ko,
         p_adj = p.adjust(p_value, method = "holm"))
write_tsv(associations, "data/gene_paff_growth_assocations.tsv")
