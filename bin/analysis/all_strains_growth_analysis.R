# Script performing analyses on the growth data from the liti paper across all strains
setwd('~/Projects/hog/')
library(tidyverse)
library(magrittr)

#### Import Data ####
## Import new growth data
growth <- read_tsv(file = 'data/raw/phenoMatrix_35ConditionsNormalizedByYPD.tab', col_names = TRUE) %>% 
  rename(strain=X1) %>%
  set_names(str_to_lower(names(.)))

# initial growth distribution plot
p_con_growth_distribution <- ggplot(gather(growth, key = 'condition', value = 'growth', -strain), aes(x=condition, y=growth)) + 
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

p_strain_growth_distribution <- ggplot(gather(growth, key = 'condition', value = 'growth', -strain), aes(x=strain, y=growth)) + 
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


# Normalise growth between strains?
# Currently normalised relative to growth in unstressed media

## Calculate distance data
load('data/genetic_distance_old.Rdata')

# remove old growth records
rm(growth_distance_NaCl4, growth_distance_NaCl6, genetic_distance_growth)

genetic_distance[lower.tri(genetic_distance)] <- NA

# filter strains with no growth data
genetic_distance <- as_tibble(genetic_distance, rownames = 'strain') %>%
  select(strain, growth$strain) %>%
  filter(strain %in% growth$strain) %>%
  gather(key = 'strain2', value = 'genetic_distance', -strain) %>%
  drop_na(genetic_distance)

# Load growth distance
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

distance <- lapply(names(growth)[-1], growth_dist) %>% bind_cols(genetic_distance, .)
# saveRDS(distance, file = 'data/distances-growth-strains.Rdata')

# distance <- readRDS('data/distances-growth-strains.Rdata')
  
## Import raw probs
probs <- read_tsv(file = 'data/hog-gene-variants.probs.blosum', col_names = TRUE) %>% filter(strain %in% growth$strain)

## Import Pathway Data and collect relavent data
path <- read_tsv(file = 'data/hog-gene-variants.path.blosum', col_names = TRUE) %>% 
  filter(strain %in% growth$strain) %>% 
  mutate(hog_active, hog_active = as.factor(hog_active))

path$total_paff <- rowSums(select(probs, -strain))
path$count <- rowSums(select(probs, -strain) > 0.5)

path <- bind_cols(path, select(growth, -strain)) %>% 
  gather(key="condition", value="growth", -strain, -hog_active, -hog_probability, -count, -total_paff)

## Import meta data
meta <- read_tsv('meta/strain_information.tsv', col_names = TRUE)

filtered_strains <- filter(meta, Ploidy == 2, Aneuploidies == 'euploid') %>% pull(`Standardized name`)

distance <- filter(distance, strain %in% filtered_strains & strain2 %in% filtered_strains)
path <- filter(path, strain %in% filtered_strains, condition %in% c("ypdnacl1m", "ypdnacl15m"))

#### NaCl Growth Against Hog Pathway ####
p_nacl_growth_vs_hog_prob <- ggplot(path, aes(x=hog_probability, y=growth, colour=condition)) + 
  geom_point() + 
  geom_smooth(method = 'lm')

ggsave('figures/liti_growth_data/hog_prob_growth.pdf', width = 12, height = 10, plot = p_nacl_growth_vs_hog_prob)

fit1 <- lm(growth ~ hog_probability, data = path, subset = path$condition == "ypdnacl1m")
fit15 <- lm(growth ~ hog_probability, data = path, subset = path$condition == "ypdnacl15m")

p_nacl_growth_vs_hog_sum <- ggplot(path, aes(x=count, y=growth, colour=condition)) + 
  geom_point() + 
  geom_smooth(method = 'lm')

ggsave('figures/liti_growth_data/hog_ko_count_growth.pdf', width = 12, height = 10, plot = p_nacl_growth_vs_hog_sum)

fit1 <- lm(growth ~ count, data = path, subset = path$condition == "ypdnacl1m")
fit15 <- lm(growth ~ count, data = path, subset = path$condition == "ypdnacl15m")

#### Distance plots
p_growth_genetic_distance <- ggplot(distance, aes(x=genetic_distance, y=ypdnacl15m_distance)) + geom_point()


#### Use ratio for distance ####
load('data/genetic_distance_old.Rdata')

# remove old growth records
rm(growth_distance_NaCl4, growth_distance_NaCl6, genetic_distance_growth)

growth_filtered <- filter(growth, strain %in% filtered_strains) %>%
  select(strain, ypdnacl15m, ypdnacl1m)

# filter strains with no growth data
genetic_distance <- as_tibble(genetic_distance, rownames = 'strain') %>%
  select(strain, growth_filtered$strain) %>%
  filter(strain %in% growth_filtered$strain) %>%
  gather(key = 'strain2', value = 'genetic_distance', -strain)

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
  bind_cols(genetic_distance, .) %>%
  gather(key = 'condition', value = 'ratio', -strain, -strain2, -genetic_distance)

p_dist_growth_ratio <- ggplot(ratios, aes(x=genetic_distance, y=ratio, colour=condition)) + 
  geom_point()



  