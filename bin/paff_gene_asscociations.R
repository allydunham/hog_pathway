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
  select(strain, condition, score)

ko_thresh <- 0.5
probs <- readRDS('data/Rdata/paff_all_genes.rds') %>%
  filter(strain %in% growth$strain) %>%
  mutate(ko = p_aff > ko_thresh)

test_gene <- function(tbl) {
  if (sum(tbl$ko) == 0) {
    return(tibble(condition=NA, n=nrow(tbl), n_ko=0, p_value=NA, mean_active = NA, mean_ko = NA))
  }
  left_join(tbl, growth, by = "strain") %>%
    group_by(condition) %>%
    summarise(p_value = ifelse(sum(ko) > 0, wilcox.test(score[ko], score[!ko])$p.value, NA),
              mean_active = mean(score[!ko], na.rm = TRUE),
              mean_ko = mean(score[ko], na.rm = TRUE),
              .groups = "drop") %>%
    mutate(n = nrow(tbl), n_ko = sum(tbl$ko)) %>%
    select(condition, n, n_ko, p_value, mean_active, mean_ko)
}

associations <- group_by(probs, gene) %>%
  group_modify(~test_gene(.x)) %>%
  mutate(n_active = n - n_ko,
         p_adj = p.adjust(p_value, method = "fdr"))
write_tsv(associations, "data/gene_paff_growth_assocations.tsv")
