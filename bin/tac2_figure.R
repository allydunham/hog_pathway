#!/usr/bin/env Rscript
# Script generating the yeast KO figure for TAC report 2

# Import packages
library(tidyverse)
library(rlang)
library(magrittr)
library(ggpubr)

#### Import General data ####
## KO Growth data
ko_growth <- read_tsv('data/raw/ko_scores.txt', col_names = TRUE) %>%
  filter(!gene == 'WT') %>%
  filter(!duplicated(.[,c('strain', 'condition', 'gene')])) %>%
  select(-position) %>%
  mutate(condition = gsub('  ', ' ', condition)) %>% # Some conditions have double spaces in names
  mutate(name = if_else(is.na(name), gene, name))

ko_growth_spread <- filter(ko_growth, qvalue < 0.01) %>%
  select(strain, condition, name, score) %>%
  spread(key = strain, value = score)

ko_growth_spread_all <- ko_growth %>%
  select(strain, condition, name, score) %>%
  spread(key = strain, value = score)
########

#### Panel 1 - genetic distance vs growth difference ####
growth_bede <- readRDS('data/Rdata/growth_bede.rds')
growth_bede <- as.matrix(select(growth_bede, -strain)) %>% set_rownames(growth_bede$strain)

growth_distance_bede <- dist(growth_bede) %>% 
  as.matrix() %>%
  as_tibble(rownames='strain1') %>%
  gather(key = 'strain2', value = 'growth_distance', -strain1)
genetic_distance_bede <- readRDS('data/Rdata/genetic_distance_matrix.rds') %>%
  as_tibble(rownames='strain1') %>%
  gather(key = 'strain2', value = 'genetic_distance', -strain1)

distances <- left_join(growth_distance_bede, genetic_distance_bede, by=c('strain1', 'strain2'))

p_distances <- ggplot(distances, aes(x=genetic_distance, y=growth_distance,
                                     group=cut_width(genetic_distance, 25000, boundary = 0))) +
  geom_point(shape=20, colour='cornflowerblue') +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  theme_pubclean() +
  ylab('S-Score Profile Distance') +
  xlab('Genetic Distance (Manhattan)')
########

#### Panel 2 - Summary of number of different genes ####
strain_ko_summary <- mutate(ko_growth, sig = (qvalue < 0.01) * sign(score)) %>%
  select(-score, -qvalue) %>%
  spread(key = strain, value = sig) %>%
  mutate(del = apply(select(., S288C, UWOP, Y55, YPS), 1, function(x){sum(x == -1, na.rm = TRUE)}),
         pos = apply(select(., S288C, UWOP, Y55, YPS), 1, function(x){sum(x == 1, na.rm = TRUE)}),
         neut = apply(select(., S288C, UWOP, Y55, YPS), 1, function(x){sum(x == 0, na.rm = TRUE)}),
         na = apply(select(., S288C, UWOP, Y55, YPS), 1, function(x){sum(is.na(x))})) %>%
  filter(na < 2) %>%
  summarise(tot = n(),
            `All Deleterious` = sum(del == 4 - na),
            `All Neutral` = sum(neut == 4 - na),
            `All beneficial` = sum(pos == 4 - na),
            `Deleterious & Beneficial Strains` = sum(del > 0 & pos > 0),
            `Beneficial & Neutral Strains` = sum(pos > 0 & neut > 0),
            `Deleterious & Neutral Strains` = sum(del > 0 & neut > 0),
            `Deleterious, Beneficial & Neutral Strains` = sum(del > 0 & neut > 0 & pos > 0)) %>%
  gather(key = 'metric', value = 'Count', -tot) %>%
  mutate(Percentage = Count / tot * 100)

p_ko_summary <- filter(strain_ko_summary, !metric == 'all_neut') %>%
  ggplot(aes(x = metric, y = Count, fill = metric)) +
  geom_col(position = position_dodge()) +
  scale_y_continuous(trans = 'pseudo_log', breaks = c(0, 1, 10, 1000)) +
  coord_flip() +
  theme(axis.title = element_blank(), axis.text.y = element_blank(),
        panel.background = element_blank(), axis.ticks.y = element_blank()) +
  guides(fill = FALSE)

p_ko_summary_tbl <- gather(strain_ko_summary, key = 'lab', value = 'value', Count, Percentage) %>%
  filter(!metric == 'all_neut') %>%
  ggplot(aes(x = lab, y = metric, label = round(value, digits = 4))) +
  geom_text() +
  scale_x_discrete(position = 'top') +
  theme(panel.background = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank())
########

#### Panel 3 - PCA of strain/condition responses? ####
########

#### Panel 4 - Examples of genes responsing very differently ####
########