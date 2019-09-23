#!/usr/bin/env Rscript
# Script generating the yeast KO figure for TAC report 2

# Import packages
library(tidyverse)
library(rlang)
library(magrittr)
library(ggpubr)

source('bin/general_functions.R')
source('bin/ko_growth_analysis/ko_analysis_functions.R')

#### Import General data ####
## KO Growth data
ko_growth <- read_tsv('data/raw/ko_scores.txt', col_names = TRUE) %>%
  filter(!gene == 'WT') %>%
  filter(!duplicated(.[,c('strain', 'condition', 'gene')])) %>%
  select(-position) %>%
  mutate(condition = gsub('  ', ' ', condition)) %>% # Some conditions have double spaces in names
  mutate(name = if_else(is.na(name), gene, name)) %>%
  mutate(sig = ifelse(qvalue < 0.01, '*', ''))

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

strain_ko_summary_condition <- mutate(ko_growth, sig = (qvalue < 0.01) * sign(score)) %>%
  select(-score, -qvalue) %>%
  spread(key = strain, value = sig) %>%
  mutate(del = apply(select(., S288C, UWOP, Y55, YPS), 1, function(x){sum(x == -1, na.rm = TRUE)}),
         pos = apply(select(., S288C, UWOP, Y55, YPS), 1, function(x){sum(x == 1, na.rm = TRUE)}),
         neut = apply(select(., S288C, UWOP, Y55, YPS), 1, function(x){sum(x == 0, na.rm = TRUE)}),
         na = apply(select(., S288C, UWOP, Y55, YPS), 1, function(x){sum(is.na(x))})) %>%
  filter(na < 2) %>%
  group_by(condition) %>%
  summarise(tot = n(),
            `All Deleterious` = sum(del == 4 - na),
            `All Neutral` = sum(neut == 4 - na),
            `All beneficial` = sum(pos == 4 - na),
            `Deleterious & Beneficial Strains` = sum(del > 0 & pos > 0),
            `Beneficial & Neutral Strains` = sum(pos > 0 & neut > 0),
            `Deleterious & Neutral Strains` = sum(del > 0 & neut > 0),
            `Deleterious, Beneficial & Neutral Strains` = sum(del > 0 & neut > 0 & pos > 0)) %>%
  gather(key = 'metric', value = 'Count', -tot, -condition) %>%
  mutate(Percentage = Count / tot * 100) %>%
  select(-Count) %>%
  spread(key = metric, value = Percentage)
########

#### Panel 3 - Relative condition sensitivity ####
condition_sensitivity <- ko_growth %>%
  filter(qvalue < 0.01) %>%
  select(strain, condition, name, score) %>%
  group_by(condition) %>%
  spread(key = strain, value = score) %>%
  summarise(S288C = sum(!is.na(S288C)),
            UWOP = sum(!is.na(UWOP)),
            Y55 = sum(!is.na(Y55)),
            YPS = sum(!is.na(YPS))) %>%
  mutate(condition = factor(condition, levels = condition[order(.$S288C/(.$S288C + .$UWOP + .$Y55 + .$YPS))], ordered = TRUE)) %>%
  gather(key = 'strain', value = 'sig_count', S288C:YPS) 

p_con_sens_bar <- ggplot(condition_sensitivity, aes(x=condition, y=sig_count, fill=strain)) +
  geom_col(position='fill') +
  theme_pubclean() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  ylab('Conditionally Significant Genes') + 
  xlab('Condition') +
  guides(fill = guide_legend(title = 'Strain'))
########

#### Panel 4 - Examples of genes responsing very differently ####
limits <- c(-11, 5)
  
clrs <- c('firebrick2', 'darkgoldenrod1', 'white', 'cornflowerblue', 'darkorchid')
vals <- scales::rescale(c(limits[1], limits[1]/2, 0, limits[2]/2, limits[2]))

plot_example_heatmaps <- function(genes, conditions, ko=ko_growth){
  tbl <- filter(ko, name %in% genes, condition %in% conditions)
  
  # Sort based on 288C
  mat <- filter(tbl, strain == 'S288C') %>%
    select(condition, name, score) %>%
    spread(key = condition, value = score) %>%
    tbl_to_matrix(row = 'name') %>%
    set_na(0)
  
  gene_dend <- hclust(dist(mat))
  con_dend <- hclust(dist(t(mat)))
  gene_order <- rownames(mat)[gene_dend$order]
  con_order <- colnames(mat)[con_dend$order]
  
  tbl <- mutate(tbl, name = factor(name, levels = gene_order),
                condition = factor(condition, levels = con_order))
  
  p <- ggplot(tbl, aes(x=condition, y=name, fill=score)) +
    facet_wrap(~strain, nrow = 1) +
    geom_raster() +
    geom_text(aes(label = sig)) +
    scale_fill_gradientn(colours = clrs, limits=limits, na.value = 'black', values = vals) +
    theme(panel.background = element_blank(), strip.background = element_blank(),
          axis.ticks = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    labs(x = '', y = '')
  return(p)
}

## Dichol genes
dichol_genes <- c('ALG9', 'GSY1', 'ALG12', 'ALG6', 'ALG8', 'ALG5', 'DIE2', 'ALG3', 'TPS1')
dichol_conditions <- c('39ºC (48H)', '6−AU + 39ºC (48H)', 'Caffeine 20mM (48H)',
                       'NaCl 0.6M (48H)', 'Acetic acid (48H)')

p_dichol_genes <- plot_example_heatmaps(dichol_genes, dichol_conditions)

## Mal genes
mal_genes <- c('MAL11', 'MAL12', 'MAL13', 'MAL31', 'MAL32', 'MAL33')
mal_condition <- 'Maltose 2% (48H)'

mal_tbl <- filter(ko_growth, condition == mal_condition, name %in% mal_genes) %>%
  mutate(name = factor(name, levels = mal_genes),
         strain = factor(strain, levels = c('S288C', 'YPS', 'Y55', 'UWOP')))

p_mal_genes <- ggplot(mal_tbl, aes(x=strain, y=name, fill=score)) +
  geom_raster() +
  geom_text(aes(label = sig)) +
  scale_fill_gradientn(colours = clrs, limits=limits, na.value = 'black', values = vals) +
  theme(panel.background = element_blank(), strip.background = element_blank(),
        axis.ticks = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        plot.title = element_text(hjust = 0.5)) +
  labs(x = '', y = '', title = mal_condition)

## NaCl genes
nacl_genes <- c('HOG1', 'PBS2', 'IXR1', 'RVS161', 'SNF4', 'VRP1', 'DUF1', 'GEM1')
nacl_conditions <- c('NaCl 0.4M (48H)', 'NaCl 0.6M (48H)', 'Maltose 2% (48H)', 'Glycerol 2% (48H)')

p_nacl_genes <- plot_example_heatmaps(nacl_genes, nacl_conditions)
########

#### Assemble plot ####

########
