# Script performing analyses on the growth data from the liti paper across all strains
setwd('~/Projects/hog/')
library(tidyverse)
library(magrittr)

#### Import Data ####
## Import new growth data
growth <- read_tsv(file = 'data/raw/phenoMatrix_35ConditionsNormalizedByYPD.tab', col_names = TRUE) %>% 
  rename(strain=X1)
names(growth) <- str_to_lower(names(growth))

# ## Load distance data
# load('data/genetic_distance.Rdata')
# 
# # remove old growth records
# rm(growth_distance_NaCl4, growth_distance_NaCl6, genetic_distance_growth)
# 
# # filter strains with no growth data
# genetic_distance <- as_tibble(genetic_distance, rownames = 'strain') %>% select(strain, growth$strain) %>% filter(strain %in% growth$strain)

## Import raw probs
probs <- read_tsv(file = 'data/hog-gene-variants.probs.blosum', col_names = TRUE) %>% filter(strain %in% growth$strain)

## Import Pathway Data and collect relavent data
path <- read_tsv(file = 'data/hog-gene-variants.path.blosum', col_names = TRUE) %>% 
  filter(strain %in% growth$strain) %>% 
  mutate(hog_active, hog_active = as.factor(hog_active))

path$growth_NaCl1mM <- growth$ypdnacl1m
path$growth_NaCl1.5mM <- growth$ypdnacl15m

path$total_paff <- rowSums(select(probs, -strain))
path$count <- rowSums(select(probs, -strain) > 0.5)

#### NaCl Growth Against Hog Pathway ####
p_nacl_growth_vs_hog_prob <- ggplot(path, aes(x=hog_probability)) + 
  geom_point(aes(y = growth_NaCl1mM, col='NaCl 1mM')) + geom_smooth(method = 'lm', aes(y = growth_NaCl1mM, col='NaCl 1mM')) + 
  geom_point(aes(y = growth_NaCl1.5mM, col='NaCl 1.5mM')) + geom_smooth(method = 'lm', aes(y = growth_NaCl1.5mM, col='NaCl 1.5mM'))

p_nacl_growth_vs_hog_prob_viol <- ggplot(path, aes(x=hog_active, y=growth_NaCl1mM)) + geom_violin()

p_nacl_growth_vs_hog_sum <- ggplot(path, aes(x=count, y=growth_NaCl1mM)) + geom_point() + geom_smooth(method = 'lm')
  
  
  
  
  