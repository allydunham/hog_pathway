setwd("~/Projects/hog/")
library(gplots)
library(tidyverse)
library(magrittr)

# Import strain meta info
meta <- read_tsv('meta/strain_information.tsv', col_names = TRUE)
strain_to_sys <- structure(meta$`Standardized name`, names=meta$`Isolate name`)
sys_to_strain <- structure(names(strain_to_sys), names = strain_to_sys)

filtered_strains <- filter(meta, Ploidy == 2, Aneuploidies == 'euploid') %>% pull(`Standardized name`)

# Import growth data
growth <- read_tsv("data/raw/bede_2017_parsed.tsv", col_names = TRUE) %>%
  mutate(strain=strain_to_sys[strain]) %>%
  drop_na(strain) %>%
#  filter(strain %in% filtered_strains) %>%
  arrange(strain)

# Import genotype data
genotypes <- read_tsv("data/hog-gene-variants.all-genotypes", col_names = TRUE) %>%
  gather(strain, val, -mut_id) %>%
  spread(mut_id, val) %>%
  filter(strain %in% growth$strain) %>% # & strain %in% filtered_strains)
  arrange(strain)

# Calculate genetic distance between strains with growth data
genetic_distance <- as.matrix(dist(select(genotypes, -strain), method = "man"))
genetic_distance[lower.tri(genetic_distance)] <- NA

genetic_distance <- as_tibble(genetic_distance) %>%
  set_names(genotypes$strain) %>%
  add_column(strain=genotypes$strain, .before = 1) %>%
  gather(key = 'strain2', value = 'genetic_distance', -strain) %>%
  drop_na(genetic_distance)
#  filter(strain %in% filtered_strains & strain2 %in% filtered_strains)
  
# Calculate growth distances and create distance table
growth_dist <- function(x){
  t <- as.matrix(dist(growth[x]))
  t[lower.tri(t)] <- NA
  
  t <- as_tibble(t) %>%
    set_names(growth$strain) %>%
    add_column(strain = growth$strain, .before = 1) %>%
    gather(key = 'strain2', value = !!x, -strain) %>%
    drop_na(!!x) %>%
    select(!!x)
  
  return(t)
}

distance <- lapply(names(growth)[-1], growth_dist) %>% bind_cols(genetic_distance, .)

## Load whole genome genetic distance
load('data/genetic_distance_old.Rdata')
#genetic_distance <- genetic_distance[rownames(genetic_distance) %in% filtered_strains,
#                                     colnames(genetic_distance) %in% filtered_strains]

rownames(genetic_distance) <- sys_to_strain[rownames(genetic_distance)]

gen_dist <- genetic_distance
gen_dist[lower.tri(gen_dist)] <- NA

gen_dist <- as_tibble(gen_dist) %>%
  select(growth$strain) %>%
  add_column(strain=colnames(genetic_distance), .before = 1) %>%
  filter(strain %in% growth$strain) %>%
  gather(key='strain2', value = 'genetic_distance_full', -strain) %>%
  drop_na()

distance %<>% rename(genetic_distance_hog=genetic_distance) %>%
  add_column(genetic_distance_full = gen_dist$genetic_distance_full, .before = 'genetic_distance_hog') %>%
  gather(key = 'condition', value = 'sscore_diff', -strain, -strain2, -genetic_distance_hog, -genetic_distance_full)

# Analyse hog distance against growth distance
p_gen_growth_dist_hog <- ggplot(filter(distance, condition %in% c("sodium chloride 0.4mM", "sodium chloride 0.6mM")),
                                aes(x=genetic_distance_hog, y=sscore_diff, colour=condition)) +
  geom_point() +
  xlab("Genetic Distance (Manhatten)") +
  ylab("Difference in S-Score")
ggsave('figures/bede_dist/genetic_dist_nacl.pdf', width = 14, height = 10, plot = p_gen_growth_dist_hog)


# Analyse whole genome distance structure
pdf('figures/heatmaps/genetic_dist_heatmap_small.pdf', width = 20, height = 20)
cols <- colorRampPalette(c("white","red"))(256)
heatmap.2(genetic_distance, symm = TRUE, revC = TRUE, col=cols,
          breaks=seq(0,max(genetic_distance),max(genetic_distance)/256), trace = "none")
dev.off()

# Genetic distance vs growth across all conds
p_gen_growth_dist_full <- ggplot(distance,
                                aes(x=genetic_distance_full, y=sscore_diff, colour=condition)) +
  geom_point() +
  xlab("Genetic Distance (Manhatten)") +
  ylab("Difference in S-Score")
ggsave('figures/bede_dist/genetic_dist_all.pdf', width = 28, height = 10, plot = p_gen_growth_dist_full)

# Try to determine difference between distribution of distances
fit <- lm(sscore_diff ~ genetic_distance_full:condition + 0, data = distance)
condition_effects <- summary(fit)$coefficients

p_gen_growth_dist_full_smooth <- ggplot(distance,
                                 aes(x=genetic_distance_full, y=sscore_diff, colour=condition)) +
  geom_smooth(method = 'lm', formula = y~x+0) +
  xlab("Genetic Distance (Manhatten)") +
  ylab("Difference in S-Score")
ggsave('figures/bede_dist/genetic_dist_all_smooth.pdf', width = 28, height = 10, plot = p_gen_growth_dist_full_smooth)

p_gen_growth_dist_full_nacl <- ggplot(filter(distance, condition %in% c("sodium chloride 0.4mM", "sodium chloride 0.6mM")),
                                 aes(x=genetic_distance_full, y=sscore_diff, colour=condition)) +
  geom_point() +
  xlab("Genetic Distance (Manhatten)") +
  ylab("Difference in S-Score")
ggsave('figures/bede_dist/genetic_dist_full_nacl.pdf', width = 14, height = 10, plot = p_gen_growth_dist_full_nacl)


# Lm coefficients are not normally distributed
#shapiro.test(fit_df$slope) 
# Shapiro-Wilk normality test
# 
# data:  fit_df$coef
# W = 0.94616, p-value = 0.04324
shapiro.test(fit_df$int) 
# Shapiro-Wilk normality test
# 
# data:  fit_df$int
# W = 0.94483, p-value = 0.03868

# Model variance of scorediff
bins=20
distance %<>% mutate(bin=cut(distance$genetic_distance_full, bins, labels = FALSE))

gen_growth_binned <- as_tibble(matrix(data = 1:bins, nrow = bins, ncol = length(unique(distance$condition)))) %>%
  set_names(unique(distance$condition)) %>%
  gather(key = 'condition', value = 'bin')

for (con in unique(gen_growth_binned$condition)){
  for (bi in 1:bins){
    gen_growth_binned[gen_growth_binned$condition == con & gen_growth_binned$bin==bi, 'sd_sscore_diff'] <- sd(
      filter(distance, condition==con, bin==bi)$sscore_diff
    )
  }
}

r <- range(distance$genetic_distance_full)
gen_growth_binned$midpoint <- r[1] + gen_growth_binned$bin * (r[2] - r[1])/(bins * 2)

p_binned_var <- ggplot(filter(gen_growth_binned, condition %in% c("sodium chloride 0.4mM","sodium chloride 0.6mM")),
                       aes(x=midpoint, y=sd_sscore_diff, col=condition)) +
  geom_point() + 
  geom_smooth(method = 'lm', formula = y~x, level=0) +
  xlab('Genetic Distance') + ylab("Standard Deviation of S-Score Differences")

fit <- lm(formula = sd_sscore_diff ~ midpoint, data = filter(gen_growth_binned, condition %in% c("sodium chloride 0.4mM")))

p_binned_sd_box <- ggplot(gen_growth_binned, aes(x=condition, y=sd_sscore_diff)) + geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ylab("Standard Deviation of S-Score Differences")

ggsave('figures/bede_dist/genetic_dist_sd_growth_box.pdf', device = 'pdf', plot = p_binned_sd_box, width = 12, height = 10)
ggsave('figures/bede_dist/genetic_dist_sd_growth.pdf', device = 'pdf', plot = p_binned_var, width = 12, height = 10)

# Model scorediff median
gen_growth_median_binned <- as_tibble(matrix(data = 1:bins, nrow = bins, ncol = length(unique(distance$condition)))) %>%
  set_names(unique(distance$condition)) %>%
  gather(key = 'condition', value = 'bin')

for (con in unique(gen_growth_median_binned$condition)){
  for (bi in 1:bins){
    gen_growth_median_binned[gen_growth_median_binned$condition == con & gen_growth_median_binned$bin==bi, 'sd_sscore_diff'] <- median(
      filter(distance, condition==con, bin==bi)$sscore_diff
    )
  }
}

gen_growth_median_binned$midpoint <- r[1] + gen_growth_median_binned$bin * (r[2] - r[1])/(bins * 2)

p_binned_median <- ggplot(filter(gen_growth_median_binned, condition %in% c("sodium chloride 0.4mM","sodium chloride 0.6mM")),
                       aes(x=midpoint, y=sd_sscore_diff, col=condition)) +
  geom_point() + 
  geom_smooth(method = 'lm', formula = y~x, level=0) +
  xlab('Genetic Distance') + ylab("Median of S-Score Differences")

fit <- lm(formula = sd_sscore_diff ~ midpoint, data = filter(gen_growth_binned, condition %in% c("sodium chloride 0.4mM")))

ggsave('figures/bede_dist/genetic_dist_median_growth.pdf', device = 'pdf', plot = p_binned_median, width = 12, height = 10)

# # S-Score normalises to condition
# p_growth_box <- ggplot(growth_melt, aes(reorder(Condition, SScore, sd), SScore)) + geom_boxplot() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1))
# 
# p_growth_dist_box <- ggplot(gen_growth_dist_melt, aes(reorder(Condition, SScore, sd), SScoreDiff)) + geom_boxplot() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1))
# 
# # Growth differences not all from same population
# kruskal.test(SScoreDiff~Condition, data = gen_growth_dist_melt)
# 
# # Apply KS-Test over different condition pairs and correct for multiple testing
# ks_tests <- sapply(unique(gen_growth_dist_melt$Condition), function(x){
#   sapply(unique(gen_growth_dist_melt$Condition), function(y){
#     if (!x == y){
#       return(ks.test(gen_growth_dist_melt[gen_growth_dist_melt$Condition==x, 'SScoreDiff'],
#                      gen_growth_dist_melt[gen_growth_dist_melt$Condition==y, 'SScoreDiff'])$p.value)
#     } else {
#       return(1)
#     }
#   })
# }) 
# dimnames(ks_tests) <- list(unique(gen_growth_dist_melt$Condition),unique(gen_growth_dist_melt$Condition))
# 
# ks_tests[lower.tri(ks_tests, diag = TRUE)] <- NA
# 
# ks_tests_melt <- subset(melt(ks_tests, value.name = 'p.value', variable.names=c('Cond1', 'Cond2')), !is.na(p.value))
# ks_tests_melt$p.adj <- p.adjust(ks_tests_melt$p.value, method = 'bonferroni')
# 
# ks_tests_melt <- ks_tests_melt[order(ks_tests_melt$p.adj),]
# ks_tests_melt_sig <- subset(ks_tests_melt, p.adj < 0.05)
