# Script to analyse whole genome P(Aff) scores and their correlations
setwd('~/Projects/hog/')
library(tidyverse)
library(magrittr)
library(gplots)

#### Import ####
meta <- read_tsv('meta/strain_information.tsv', col_names = TRUE)

genes <- read_tsv('meta/sacc_gene_loci', col_names = TRUE)
sys_to_gene <- structure(genes$name, names=genes$id)


probs <- read_tsv('data/all-genes-no-missing.koprob', col_names = TRUE)

p_aff_sums <- colSums(select(probs, -strain))
no_prob_genes <- names(which(p_aff_sums == 0))

probs %<>% select(-one_of(no_prob_genes))

#### Analysis ####
cor <- cor(select(probs, -strain))

cor[lower.tri(cor, diag = TRUE)] <- NA

cor_melt <- as_tibble(cor) %>%
  add_column(gene = rownames(cor), .before = 1) %>%
  gather(key = 'gene2', value = 'cor', -gene) %>%
  drop_na(cor) %>%
  mutate(mag = abs(cor)) %>%
  arrange(desc(mag)) %>%
  mutate(name = sys_to_gene[gene], name2 = sys_to_gene[gene2])

# pdf('figures/all_genes_paff_heatmap.pdf', width = 50, height = 50)
# cols <- colorRampPalette(c("blue", "white","red"))(256)
# heatmap.2(cor, symm = TRUE, revC = TRUE, col=cols,
#           breaks=seq(-1,1,2/256), trace = "none")
# dev.off()
