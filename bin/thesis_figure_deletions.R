#!/us/bin/env Rscript
# Generate figure on deletion phenotypes for thesis
library(tidyverse)
library(ggpubr)
library(multipanelfigure)
library(ggtext)

theme_set(theme_pubclean() + theme(legend.position = 'right',
                                   plot.title = element_text(hjust = 0.5),
                                   plot.subtitle = element_text(hjust = 0.5),
                                   strip.background = element_blank(),
                                   legend.key = element_blank()))

blank_plot <- function(text = ''){
  ggplot(tibble(x=c(0, 1)), aes(x=x, y=x)) +
    geom_blank() +
    annotate(geom = 'text', x = 0.5, y = 0.5, label = text) +
    theme(panel.grid.major.y = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank())
}

#### Load Data ####
gene_meta <- readRDS('data/Rdata/gene_meta_all.rds') %>%
  mutate(name = ifelse(is.na(name), id, name)) %>%
  rename(gene = id)

ko_growth <- read_tsv('data/raw/ko_scores.txt', col_names = TRUE) %>%
  filter(!gene == 'WT') %>%
  filter(!duplicated(.[,c('strain', 'condition', 'gene')])) %>%
  select(-position) %>%
  mutate(condition = gsub('  ', ' ', condition)) %>% # Some conditions have double spaces in names
  mutate(name = if_else(is.na(name), gene, name))

# geno <- readRDS('data/Rdata/genotypes_all_genes.rds') %>%
#   select(UWOP = SACE_GAT, Y55 = ADA, YPS = AKN, S288C = SACE_GAV) %>%
#   mutate(across(everything(), ~ifelse(. > 0, 1, 0))) %>%
#   as.matrix() %>%
#   t()
# geno_dist <- dist(geno, method = "manhattan")

#### Panel - Strain Relationships - genetic and phenotypic (confusion) ####
# dendrgrams from genetic distance and KO profiles
growth_dists <- select(ko_growth, condition, strain, name, score) %>%
  pivot_wider(names_from = strain, values_from = score) %>%
  group_by(condition) %>%
  summarise(S288C_UWOP = sqrt(sum((S288C - UWOP)^2, na.rm = TRUE)),
            S288C_Y55 = sqrt(sum((S288C - Y55)^2, na.rm = TRUE)),
            S288C_YPS = sqrt(sum((S288C - YPS)^2, na.rm = TRUE)),
            Y55_UWOP = sqrt(sum((Y55 - UWOP)^2, na.rm = TRUE)),
            Y55_YPS = sqrt(sum((Y55 - YPS)^2, na.rm = TRUE)),
            YPS_UWOP = sqrt(sum((YPS - UWOP)^2, na.rm = TRUE))) %>%
  pivot_longer(-condition, names_to = "pair", values_to = "distance")

ggplot(growth_dists, aes(x = pair, y = distance)) +
  geom_boxplot()

p_strains <- blank_plot("Strain dendrograms\n(genetic/growth)")

#### Panel - Mal gene example heatmap ####
# MAL genes noticed first
p_mal <- blank_plot("MAL Gene heatmap")

#### Panel - Other example heatmap ####
# A second set of example genes showing switching behavior
p_example <- blank_plot("Example genes heatmap")

#### Panel - Condition correlation ####
# Conditions that correlate most between strains
p_condition_cor <- blank_plot("Condition correlation between strains")

#### Panel - Condition pair correlation ####
# Condition pair correlations within strains
p_pair_cor <- blank_plot("Condition pair correlations")

#### Panel - Relative Sensitivity per strain ####
# Which strain has most sig knockouts per condition
p_sensitivity <- blank_plot("Strain sensitivity to KO")

#### Figure Assembly ####
size <- theme(text = element_text(size = 8))
p1 <- p_strains + labs(tag = 'A') + size
p2 <- p_mal + labs(tag = 'B') + size
p3 <- p_example + labs(tag = 'C') + size
p4 <- p_condition_cor + labs(tag = 'D') + size
p5 <- p_pair_cor + labs(tag = 'E') + size
p6 <- p_sensitivity + labs(tag = 'E') + size

figure <- multi_panel_figure(width = 180, height = 120, columns = 3, rows = 2,
                             panel_label_type = 'none', row_spacing = 0, column_spacing = 0) %>%
  fill_panel(p1, row = 1, column = 1) %>%
  fill_panel(p2, row = 1, column = 2) %>%
  fill_panel(p3, row = 1, column = 3) %>%
  fill_panel(p4, row = 2, column = 1) %>%
  fill_panel(p5, row = 2, column = 2) %>%
  fill_panel(p6, row = 2, column = 3)

ggsave('figures/thesis_figure_deletions.pdf', figure, width = figure_width(figure), height = figure_height(figure), units = 'mm')
ggsave('figures/thesis_figure_deletions.tiff', figure, width = figure_width(figure), height = figure_height(figure), units = 'mm')
