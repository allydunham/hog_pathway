# Script to analyse whole genome P(Aff) scores and their correlations
setwd('~/Projects/hog/')
library(MASS)
library(mixtools)
library(gplots)
library(tidyverse)
library(magrittr)

#### Import ####
## Strain information
meta <- read_tsv('meta/strain_information.tsv', col_names = TRUE)


## Gene information
genes <- read_tsv('meta/sacc_gene_loci', col_names = TRUE)
sys_to_gene <- structure(genes$name, names=genes$id)
filtered_strains <- filter(meta, Ploidy == 2, Aneuploidies == 'euploid') %>% pull(`Standardized name`)

hog_genes <- read_table2('meta/hog-gene-loci', col_names = FALSE, comment = '#') %>%
  set_names(c('chrom', 'start', 'stop', 'id', 'name', 'strand'))


# Growth data
growth <- read_tsv(file = 'data/raw/phenoMatrix_35ConditionsNormalizedByYPD.tab', col_names = TRUE) %>% 
  rename(strain=X1) %>%
  set_names(str_to_lower(names(.)))

## P(aff)
probs <- read_tsv('data/all-genes-no-missing.koprob', col_names = TRUE) %>%
  filter(strain %in% growth$strain)

# Filter genes with no variation
p_aff_sums <- colSums(select(probs, -strain))
no_prob_genes <- names(which(p_aff_sums == 0))
probs %<>% select(-one_of(no_prob_genes))

#### Analysis ####
## Gene/Gene Correlation
cor_genes <- cor(select(probs, -strain))

cor_genes[lower.tri(cor_genes, diag = TRUE)] <- NA

cor_genes_melt <- as_tibble(cor_genes) %>%
  add_column(gene = rownames(cor_genes), .before = 1) %>%
  gather(key = 'gene2', value = 'cor', -gene) %>%
  drop_na(cor) %>%
  mutate(mag = abs(cor)) %>%
  arrange(desc(mag)) %>%
  mutate(name = sys_to_gene[gene], name2 = sys_to_gene[gene2])

# pdf('figures/all_genes_paff_heatmap.pdf', width = 50, height = 50)
# cols <- colorRampPalette(c("blue", "white","red"))(256)
# heatmap.2(cor_genes, symm = TRUE, revC = TRUE, col=cols,
#           breaks=seq(-1,1,2/256), trace = "none")
# dev.off()

## Gene/Growth Correlation
cor_growth <- cor(select(probs, -strain), select(growth, -strain))
meanGrowthCor <- rowMeans(cor_growth)

pdf('figures/all_genes_growth_heatmap.pdf', width = 50, height = 50)
cols <- colorRampPalette(c("blue", "white","red"))(256)
heatmap.2(cor_growth, col=cols, breaks=seq(-1,1,2/256), trace = "none", symkey = FALSE)
dev.off()

# normalise by subtraction
pdf('figures/all_genes_growth_heatmap_norm_sub.pdf', width = 50, height = 50)
cols <- colorRampPalette(c("blue", "white","red"))(256)
heatmap.2(cor_growth - meanGrowthCor, col=cols, breaks=seq(-1,1,2/256), trace = "none", symkey = FALSE)
dev.off()

# Normalise by division
pdf('figures/all_genes_growth_heatmap_norm_div.pdf', width = 50, height = 50)
cols <- colorRampPalette(c("blue", "white","red"))(256)
heatmap.2(cor_growth / meanGrowthCor, col=cols, trace = "none", symkey = FALSE)
dev.off()

# Hog Genes Only
pdf('figures/hog_genes_growth_heatmap.pdf', width = 50, height = 50)
cols <- colorRampPalette(c("blue", "white","red"))(256)
heatmap.2(cor_growth[rownames(cor_growth) %in% hog_genes$id,], breaks=seq(-1,1,2/256), 
          col=cols, trace = "none", symkey = FALSE, cexRow = 2.5, cexCol = 2.5, margins = c(24,15))
dev.off()

cor_growth_melt <- as_tibble(cor_growth, rownames = 'gene_id') %>%
  gather(key='condition', value = 'correlation', -gene_id) %>%
  group_by(condition) %>%
  mutate(gene_name = sys_to_gene[gene_id]) %>%
  mutate(mag = abs(correlation)) %>% 
  mutate(cor_sub_mean = correlation - meanGrowthCor) %>%
  mutate(cor_div_mean = correlation/meanGrowthCor) %>%
  arrange(desc(mag), .by_group = TRUE)
  

# Genes analysis
gene_p_aff_means <- colMeans(select(probs, -strain))

# prob_melt <- gather(probs, key = 'gene', value = 'p_aff', -strain)
# 
# prob_mat <- as.matrix(select(probs, -strain))
# 
# fit <- normalmixEM(probs$YAL038W, mu = c(0,0.6))
# hist(prob_mat[!prob_mat == 0], probability = TRUE)
# curve(dexp(x, rate = fit$estimate), col='red', add=TRUE)

gene_summary <- tibble(id = names(select(probs, -strain))) %>%
  mutate(mean = colMeans(select(probs, -strain))) %>%
  mutate(mean_non_zero = apply(prob_mat, 2, function(x){mean(x[!x == 0])})) %>%
  mutate(sum_zero = apply(prob_mat, 2, function(x){sum(x==0)})) %>%
  mutate(sum_non_zero = apply(prob_mat, 2, function(x){sum(!x==0)})) %>%
  left_join(genes, by = 'id') %>%
  rename(chrom=`#chrom`) %>%
  mutate(length=abs(stop - start))

p_length <- ggplot(gene_summary, aes(x=length, y=mean_non_zero)) +
  geom_point() +
  xlab('Gene Length') +
  ylab('Mean P(Aff)')
  
ggsave('figures/correlations/length_vs_mean_paff.pdf', width = 12, height = 10, plot = p_length)


