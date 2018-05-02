# Script looking at ko probability of genes comparaed to growth rates (based on Bede's growth screen)
setwd('~/Projects/hog/')
library(tidyverse)
library(magrittr)
library(purrr)
library(gplots)

#### Import ####
# Meta info
hog_meta <- readRDS('data/Rdata/gene_meta_hog.rds')
strain_meta <- readRDS('data/Rdata/strain_meta.rds')
strain_to_sys <- structure(strain_meta$`Standardized name`, names=strain_meta$`Isolate name`)

# Growth data
growth <- readRDS('data/Rdata/growth_bede.rds') %>%
  select(strain, `sodium chloride 0.4mM`, `sodium chloride 0.6mM`) %>%
  gather(key = 'condition', value = 'sscore', -strain)

# KO probabilities
prob_aff_hog <- readRDS('data/Rdata/paff_hog_genes.rds')

prob_aff_not_hog <- readRDS('data/Rdata/paff_all_genes.rds') %>%
  filter(strain %in% growth$strain, !gene %in% hog_meta$id)


#### Analysis ####
## Gene P(Aff) distribution
p_paff_box <- ggplot(prob_aff_hog, aes(x=gene, y=p_aff)) + geom_boxplot() + ylab('P(Aff)') + xlab('Gene')
ggsave('figures/p_affected_genes_box_corrected.pdf', width = 12, height = 10, plot = p_paff_box)

## Gene Correlations ##
gene_cor <- cor(as.matrix(spread(prob_aff_hog, key = gene, value = p_aff) %>% select(-strain)))

pdf('figures/heatmaps/gene_ko_probs_heatmap_corrected.pdf', width = 12, height = 12)
cols <- colorRampPalette(c("red","white","blue"))(256)
heatmap.2(gene_cor, symm=TRUE, col=cols, breaks=seq(-1,1,2/256), trace = "none", revC = TRUE, cellnote = round(gene_cor, 2))
dev.off()

## p(aff) against growth
# Filter to strains with growth data
prob_aff_hog %<>% filter(strain %in% growth$strain)

## Test number of ko'd genes
thresh <- 0.5
prob_aff_hog %<>% mutate(thresh=p_aff>thresh)

growth$ko_count <- NA
growth$p_aff_hog_sum <- NA
growth$p_aff_non_hog_sum <- NA
growth$p_aff_hog_mean <- NA
growth$p_aff_non_hog_mean <- NA
for (strain in unique(growth$strain)){
  growth[growth$strain == strain, 'ko_count'] <- sum(filter(prob_aff_hog, strain == !!strain)$thresh)
  growth[growth$strain == strain, 'p_aff_hog_sum'] <- sum(filter(prob_aff_hog, strain == !!strain)$p_aff)
  growth[growth$strain == strain, 'p_aff_non_hog_sum'] <- sum(filter(prob_aff_not_hog, strain == !!strain)$p_aff)
  growth[growth$strain == strain, 'p_aff_hog_mean'] <- mean(filter(prob_aff_hog, strain == !!strain)$p_aff)
  growth[growth$strain == strain, 'p_aff_non_hog_mean'] <- mean(filter(prob_aff_not_hog, strain == !!strain)$p_aff)
}

# Normalise to the number of genes
growth %<>% mutate(p_aff_hog_sum_norm = p_aff_hog_sum/length(unique(prob_aff_hog$gene)),
                   p_aff_non_hog_sum_norm = p_aff_non_hog_sum/length(unique(prob_aff_not_hog$gene)))

p_ko_count <- ggplot(growth, aes(x=ko_count, y=sscore, col=condition)) +
                       geom_point() + 
                       geom_smooth(method='lm', formula = y~x) + 
                       ggtitle(paste0('Effect of number of likely KOs on growth (Threshold: p(Aff) = ', thresh,')')) + 
                       xlab('KO Count') + ylab('S Score')

ggsave('figures/bede_growth/ko-count-vs-growth.pdf', p_ko_count, width = 12, height = 10)

fit <- lm(sscore ~ ko_count, data = filter(growth, condition=='sodium chloride 0.6mM'))
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -3.1501 -0.6994  0.1727  0.7093  2.3411 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)   
# (Intercept)  0.12803    0.20454   0.626  0.53296   
# ko_count    -0.09226    0.03414  -2.702  0.00827 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 1.035 on 88 degrees of freedom
# Multiple R-squared:  0.07661,	Adjusted R-squared:  0.06612 
# F-statistic: 7.301 on 1 and 88 DF,  p-value: 0.008267

p_ko_sum <- ggplot(growth, aes(x=p_aff_hog_sum, y=sscore, col=condition)) + 
  geom_point() + 
  geom_smooth(method='lm', formula = y~x) + 
  xlab('Sum(P(Aff))') + ylab('S-Score')

ggsave('figures/bede_growth/sum_p_aff-vs-growth.pdf', p_ko_sum, width = 12, height = 10)

fit <- lm(sscore ~ p_aff_hog_sum, data = filter(growth, condition=='sodium chloride 0.6mM'))

p_ko_sum_hog <- ggplot(growth, aes(y=sscore)) + 
  geom_point(aes(x=p_aff_non_hog_mean, colour='Not Hog')) + 
  geom_smooth(aes(x=p_aff_non_hog_mean, colour='Not Hog'), method='lm', formula = y~x) + 
  geom_point(aes(x=p_aff_hog_mean, colour='Hog')) + 
  geom_smooth(aes(x=p_aff_hog_mean, colour='Hog'), method='lm', formula = y~x) +
  xlab('Sum(P(Aff))') + ylab('S-Score')


## Test association against high probability ko
high_conf_ko <- read_tsv('data/hog-gene-variants.conf-ko0.9', col_names = TRUE) %>%
  select(strain, hog_meta$id) %>%
  set_names(c('strain', hog_meta$name)) %>%
  filter(strain %in% growth$strain) %>%
  gather(key = 'gene', value = 'ko', -strain)

growth$high_conf_ko_count <- NA
for (strain in unique(growth$strain)){
  growth[growth$strain == strain, 'high_conf_ko_count'] <- sum(filter(high_conf_ko, strain == !!strain)$ko)
}

p_conf_ko_count <- ggplot(growth, aes(x=high_conf_ko_count, y=sscore, col=condition)) + 
  geom_point() + 
  geom_smooth(method='lm', formula = y~x) + 
  ggtitle('Effect of number of high confidence KOs on growth (Threshold: Proportion = 0.9)') + 
  xlab('KO Count') + ylab('S Score')

ggsave('figures/bede_growth/conf-ko-count-vs-growth.pdf', p_conf_ko_count, width = 12, height = 10)

## Test if presence of any high conf is important
growth %<>% mutate(any = high_conf_ko_count > 0)

p_conf_ko_box <- ggplot(growth, aes(x=any, y=sscore)) + 
  geom_boxplot() + 
  ggtitle('Distribution of S Score with HOG knockouts') + 
  xlab('KOs') + ylab('S Score')

ggsave('figures/bede_growth/conf-kos-vs-growth-box.pdf', p_conf_ko_box, width = 12, height = 10)

t.test(sscore~any, data = growth, alternative='l')
# Welch Two Sample t-test
# 
# data:  sscore by any
# t = 0.3358, df = 43.721, p-value = 0.7386
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.3240032  0.4535324
# sample estimates:
#   mean in group FALSE  mean in group TRUE 
# -0.2761068          -0.3408714


#### Analyse whole pathway ko prob ####
path_active <- read_tsv('data/hog-gene-variants.path.blosum', col_names = TRUE) %>%
  filter(strain %in% growth$strain) %>%
  left_join(., growth, by = 'strain') %>%
  mutate(hog_active=as.logical(hog_active))

count_categories <- function(x){
  return(c(y = -3.5, label=length(x)))
}

p_bin_path_growth <- ggplot(path_active, aes(x=hog_active, y=sscore)) + 
  geom_boxplot() + 
  facet_wrap(~condition) + 
  xlab('HOG Path Active') + ylab('S-Score')

ggsave('figures/bede_growth/bin-path-vs-growth.pdf', width = 12, height = 10, plot = p_bin_path_growth)

t.test(sscore ~ hog_active, data = path_active, alternative='less')

## Hog1 activity probability against growth
p_prob_path_growth <- ggplot(path_active, aes(x=hog_probability, y=sscore, colour=condition)) + 
  geom_point() + 
  geom_smooth(method = 'lm') + 
  xlab('Hog1 Activity Probability') + 
  ylab('S-Score')

ggsave('figures/bede_growth/path-probability-vs-growth.pdf', width = 12, height = 10, plot = p_prob_path_growth)

fit <- lm(sscore~hog_probability, data = path_active, subset = path_active$condition == 'sodium chloride 0.6mM')

## growth against growth
p_growth_cor <- ggplot(spread(path_active, key = 'condition', value = 'sscore'), aes(x=`sodium chloride 0.4mM`, y=`sodium chloride 0.6mM`)) +
  geom_point() +
  geom_smooth(method = 'lm') + 
  xlab('S-Score (NaCl 0.4mM)') + 
  ylab('S-Score (NaCl 0.6mM)')

ggsave('figures/bede_growth/growth_correlation_nacl.pdf', width = 12, height = 10, plot = p_growth_cor)

# Strong correlation, increase in .6mM is slower
fit <- lm(`sodium chloride 0.6mM` ~ `sodium chloride 0.4mM`, spread(path_active, key = 'condition', value = 'sscore'))
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.88429 -0.25206 -0.00333  0.24655  0.86260 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)             -0.12534    0.03964  -3.162  0.00215 ** 
#   `sodium chloride 0.4mM`  0.91881    0.03557  25.833  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.3678 on 88 degrees of freedom
# Multiple R-squared:  0.8835,	Adjusted R-squared:  0.8822 
# F-statistic: 667.3 on 1 and 88 DF,  p-value: < 2.2e-16

                                 
