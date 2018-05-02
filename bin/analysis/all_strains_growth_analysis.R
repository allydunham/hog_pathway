# Script performing analyses on the growth data from the liti paper across all strains
setwd('~/Projects/hog/')
library(tidyverse)
library(magrittr)
library(preprocessCore)

#### Import Data ####
## Meta Data
genes <- readRDS('data/Rdata/gene_meta_all.rds')

strains <- readRDS('data/Rdata/strain_meta.rds')
filtered_strains <- filter(strains, Ploidy == 2, Aneuploidies == 'euploid') %>% pull(`Standardized name`)

## Growth data
growth <- readRDS('data/Rdata/growth_liti.rds')
  
# initial growth distribution plot
p_con_growth_distribution <- ggplot(gather(growth, key = 'condition', value = 'growth', -strain), aes(x=condition, y=growth)) + 
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

p_strain_growth_distribution <- ggplot(gather(growth, key = 'condition', value = 'growth', -strain), aes(x=strain, y=growth)) + 
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Normalise growth between strains?
# Currently normalised relative to growth in unstressed media

# Filter strains
growth %<>% filter(strain %in% filtered_strains)

## Genetic Distance
genetic_distance <- readRDS('data/Rdata/genetic_distance.rds') %>%
  filter(strain %in% intersect(growth$strain, filtered_strains),
         strain2 %in% intersect(growth$strain, filtered_strains)) %>%
  rename(genetic_distance=distance)

genetic_distance_matrix <- readRDS('data/Rdata/genetic_distance_matrix.rds')

## Growth Distance
growth_dist <- function(x){
  name <- paste0(x, '_distance')
  t <- as.matrix(dist(growth[x]))
  t[lower.tri(t)] <- NA

  t <- as_tibble(t) %>%
    set_names(growth$strain) %>%
    add_column(strain = growth$strain, .before = 1) %>%
    gather(key = 'strain2', value = !!name, -strain) %>%
    drop_na(!!name) %>%
    select(!!name)

  return(t)
}

distance <- lapply(names(growth)[-1], growth_dist) %>% 
  bind_cols(genetic_distance, .) %>%
  gather(key = 'condition', value = 'growth', -strain, -strain2, -genetic_distance)

## Raw probs
probs_hog <- readRDS('data/Rdata/paff_hog_genes.rds') %>%
  spread(key = 'gene', value = 'p_aff')

## Pathway Data and collect relavent data
path <- readRDS('data/Rdata/hog_path_probs.rds') %>%
  mutate(total_paff = rowSums(select(probs_hog, -strain))) %>%
  mutate(count = rowSums(select(probs_hog, -strain) > 0.5)) %>%
  filter(strain %in% intersect(growth$strain, filtered_strains))

path <- bind_cols(path, select(growth, -strain)) %>% 
  gather(key="condition", value="growth", -strain, -hog_active, -hog_probability, -count, -total_paff)

## Mutation counts
counts <- readRDS('data/Rdata/hog_gene_mut_counts.rds') %>%
  filter(strain %in% intersect(growth$strain, filtered_strains))

#### Analysis ####
## NaCl Growth Against Hog Pathway
p_nacl_growth_vs_hog_prob <- ggplot(path, aes(x=hog_probability, y=growth, colour=condition)) + 
  geom_point() + 
  geom_smooth(method = 'lm')

ggsave('figures/liti_growth/hog_prob_growth.pdf', width = 12, height = 10, plot = p_nacl_growth_vs_hog_prob)

fit1 <- lm(growth ~ hog_probability, data = path, subset = path$condition == "ypdnacl1m")
fit15 <- lm(growth ~ hog_probability, data = path, subset = path$condition == "ypdnacl15m")

p_nacl_growth_vs_hog_sum <- ggplot(path, aes(x=count, y=growth, colour=condition)) + 
  geom_point() + 
  geom_smooth(method = 'lm')

ggsave('figures/liti_growth/hog_ko_count_growth.pdf', width = 12, height = 10, plot = p_nacl_growth_vs_hog_sum)

fit1 <- lm(growth ~ count, data = path, subset = path$condition == "ypdnacl1m")
fit15 <- lm(growth ~ count, data = path, subset = path$condition == "ypdnacl15m")

## Distance plots
p_growth_genetic_distance <- ggplot(filter(distance, condition %in% c("ypdnacl1m_distance", "ypdnacl15m_distance")),
                                    aes(x=genetic_distance, y=growth)) + 
  geom_point() + facet_wrap(~condition) + 
  xlab("Genetic Distance") + ylab("Growth")

ggsave('figures/liti_dist/genetic_dist_nacl_liti.pdf', plot = p_growth_genetic_distance, width = 16, height = 10)

fit <- lm(growth ~ genetic_distance, filter(distance, condition %in% c("ypdnacl15m_distance")))

## Binned SD analysis
bins = 40
distance %<>% mutate(bin=cut(distance$genetic_distance, bins, labels = FALSE))

gen_growth_binned <- as_tibble(matrix(data = 1:bins, nrow = bins, ncol = length(unique(distance$condition)))) %>%
  set_names(unique(distance$condition)) %>%
  gather(key = 'condition', value = 'bin')

gen_growth_binned$sd_growth_diff <- NA

for (con in unique(gen_growth_binned$condition)){
  for (bi in 1:bins){
    gen_growth_binned[gen_growth_binned$condition == con & gen_growth_binned$bin==bi, 'sd_growth_diff'] <- sd(
      filter(distance, condition==con, bin==bi)$growth
    )
  }
}

r <- range(distance$genetic_distance)
gen_growth_binned$midpoint <- r[1] + gen_growth_binned$bin * (r[2] - r[1])/(bins * 2)

p_binned_var <- ggplot(filter(gen_growth_binned, condition %in% c("ypdnacl1m_distance","ypdnacl15m_distance")),
                       aes(x=midpoint, y=sd_growth_diff, col=condition)) +
  geom_point() + 
  geom_smooth(method = 'lm', formula = y~x, level=0) +
  xlab('Genetic Distance') + ylab("Standard Deviation of Growth Differences")

fit <- lm(formula = sd_growth_diff ~ midpoint, data = filter(gen_growth_binned, condition %in% c("ypdnacl1m_distance")))

p_binned_sd_box <- ggplot(gen_growth_binned, aes(x=condition, y=sd_growth_diff)) + geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ylab("Standard Deviation of Growth Differences")

ggsave('figures/liti_dist/genetic_dist_sd_growth_box_liti.pdf', device = 'pdf', plot = p_binned_sd_box, width = 12, height = 10)
ggsave('figures/liti_dist/genetic_dist_sd_growth_liti.pdf', device = 'pdf', plot = p_binned_var, width = 12, height = 10)

# Using Median
gen_growth_median_binned <- as_tibble(matrix(data = 1:bins, nrow = bins, ncol = length(unique(distance$condition)))) %>%
  set_names(unique(distance$condition)) %>%
  gather(key = 'condition', value = 'bin')

gen_growth_median_binned$med_growth_diff <- NA

for (con in unique(gen_growth_median_binned$condition)){
  for (bi in 1:bins){
    gen_growth_median_binned[gen_growth_median_binned$condition == con & gen_growth_median_binned$bin==bi, 'med_growth_diff'] <- median(
      filter(distance, condition==con, bin==bi)$growth
    )
  }
}

r <- range(distance$genetic_distance)
gen_growth_median_binned$midpoint <- r[1] + gen_growth_median_binned$bin * (r[2] - r[1])/(bins * 2)

p_binned_med <- ggplot(filter(gen_growth_median_binned, condition %in% c("ypdnacl1m_distance","ypdnacl15m_distance")),
                       aes(x=midpoint, y=med_growth_diff, col=condition)) +
  geom_point() + 
  geom_smooth(method = 'lm', formula = y~x, level=0) +
  xlab('Genetic Distance') + ylab("Median of Growth Difference")

fit <- lm(formula = med_growth_diff ~ midpoint, data = filter(gen_growth_median_binned, condition %in% c("ypdnacl15m_distance")))

ggsave('figures/liti_dist/genetic_dist_median_growth_liti.pdf', device = 'pdf', plot = p_binned_med, width = 12, height = 10)


#### Use ratio for distance ####
growth_filtered <- select(growth, strain, ypdnacl15m, ypdnacl1m) %>%
  filter(!(ypdnacl15m == 0), !(ypdnacl1m == 0))

genetic_distance_filtered <- as_tibble(genetic_distance_matrix, rownames = 'strain') %>%
  gather(key = 'strain2', value = 'genetic_distance', -strain) %>%
  filter(strain %in% growth_filtered$strain, strain2 %in% growth_filtered$strain)

growth_ratio <- function(x){
  name <- paste0(x, '_ratio')
  gro <- structure(growth_filtered[[x]], names=growth_filtered$strain)

  t <- sapply(gro, function(i){
    sapply(gro, function(j){
      i/j # Gives col/row in matrix
    })
  })
  
  t <- as_tibble(t) %>%
    add_column(strain = growth_filtered$strain, .before = 1) %>%
    gather(key = 'strain2', value = !!name, -strain) %>%
    select(!!name)
  
  return(t)
}

ratios <- lapply(names(growth_filtered)[-1], growth_ratio) %>% 
  bind_cols(genetic_distance_filtered, .) %>%
  gather(key = 'condition', value = 'ratio', -strain, -strain2, -genetic_distance)

p_dist_growth_ratio <- ggplot(ratios, aes(x=genetic_distance, y=ratio, colour=condition)) + 
  geom_point()

## Analyse mutation counts
counts_melt <- gather(counts, key = 'gene', value = 'count', -strain) %>%
  mutate(ypdnacl1m=structure(growth$ypdnacl1m, names=growth$strain)[strain]) %>%
  mutate(ypdnacl15m=structure(growth$ypdnacl15m, names=growth$strain)[strain]) %>%
  mutate(ypdkcl2m=structure(growth$ypdkcl2m, names=growth$strain)[strain]) %>%
  mutate(ypdchx1=structure(growth$ypdchx1, names=growth$strain)[strain])
  
p_gene_mut_counts <- ggplot(counts_melt,
                            aes(x=count, y=ypdnacl15m)) + 
  geom_boxplot(aes(group=cut_width(count, 1)), varwidth = TRUE) +
  facet_wrap(~ gene)

  